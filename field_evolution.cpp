#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "common.h"
#include "initial_conditions.h"
#include "field_evolution.h"
#include "boundary_conditions.h"

inline double step(double x) noexcept {
    return x >= 0.0;
}

inline double MC(double a, double b, double c) noexcept
{
    if ((a > 0.0) && (b > 0.0) && (c > 0.0))
        return std::min({a, b, c});

    if ((a < 0.0) && (b < 0.0) && (c < 0.0))
        return std::max({a, b, c});

    return 0.0;
}

/*
    Computes current times 4*pi/c*(lapse function) components along cell edges in reduced units i.e., 4*pi*L_0/(c*B_0)*exp(nu/2)*J
    Works for both magnetizable B != H and non-magnetizable medium

    Inputs: B: magnetic field in reduced units at cell centers
            T: redshifted temperature in reduced units
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process objects containing information about the simulation domain and the current process
    Output: J: redshifted current density in reduced units
*/
void Compute_J(VectorField & B, VectorField & J, ScalarField & T, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process)
{
    double Deltatheta = dm.Deltatheta;
    double m1_r, m2_r, m4_r; //magnetization-dependent prefactors for the curl(B), grad(B)xBhat and grad(T)xBhat terms respectively evaluated at location of J_r
    double m1_th, m2_th, m3_th, m4_th; //magnetization-dependent prefactors for the curl(B), grad(B)xBhat, grad(mu_e)xBhat and grad(T)xBhat terms respectively evaluated at location of J_theta
    double m1_phi, m2_phi, m3_phi, m4_phi; //magnetization-dependent prefactors for the curl(B), grad(B)xBhat, grad(mu_e)xBhat and grad(T)xBhat terms respectively evaluated at location of J_phi
    double gradrB, gradthB, gradrmu_e, gradrT, gradthT; //components of gradients of Bmag, mu_e and T
    double gradrB_2, gradthB_2, gradrT_2, gradthT_2; //components of gradients of Bmag, mu_e and T. Used for computing value at lower left cell corners (location of J_phi^{i,j}). Don't need gradrmu_e_2 because mu is assumed spherically symmetric
    double r_imh, Deltar_n1, n_t_av, n_r_av;

    if(mC.quantization == true){
        double J_B_r, J_B_theta, J_B_phi; //curl(B) components
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
            Deltar_n1 = 0.5*( dm.Deltar[i-1] + dm.Deltar[i] );
            n_t_av = 0.5*( dm.n_t[i-1] + dm.n_t[i] );
            n_r_av = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            gradrmu_e = (mC.mu_e[i]-mC.mu_e[i-1])/Deltar_n1/n_r_av;

            double MagFactor = 1.;
            for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                if(mC.PolarAxisQuantization == true && ( dm.theta[j] - dm.Deltatheta < 0 || dm.theta[j] + dm.Deltatheta > pi )){
                    MagFactor = 0.;
                }
                else MagFactor = 1.;

                m1_r = ( 1. - MagFactor*4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) );
                m1_th = ( 1. - MagFactor*4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i-1][j]/mC.Bmag[i-1][j] ) );
                m1_phi = ( 1. - MagFactor*4.*pi*0.25*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] + mC.M[i-1][j]/mC.Bmag[i-1][j] + mC.M[i-1][j-1]/mC.Bmag[i-1][j-1] ) );
                m2_r = -MagFactor*4.*pi*( 0.5*( mC.chi[i][j] + mC.chi[i][j-1] ) - 0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) )*dm.n_t[i];
                m2_th = -MagFactor*4.*pi*( 0.5*( mC.chi[i][j]*dm.n_t[i] + mC.chi[i-1][j]*dm.n_t[i-1] ) - 0.5*( mC.M[i][j]/mC.Bmag[i][j]*dm.n_t[i] + mC.M[i-1][j]/mC.Bmag[i-1][j]*dm.n_t[i-1] ) );
                m2_phi = -MagFactor*4.*pi*( 0.25*( ( mC.chi[i][j] + mC.chi[i][j-1] )*dm.n_t[i] + ( mC.chi[i-1][j] + mC.chi[i-1][j-1] )*dm.n_t[i-1] )
                                 - 0.25*( ( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] )*dm.n_t[i] + ( mC.M[i-1][j]/mC.Bmag[i-1][j] + mC.M[i-1][j-1]/mC.Bmag[i-1][j-1] )*dm.n_t[i-1] ) );
                m3_th = -MagFactor*4.*pi*0.5*( mC.M_mu[i][j]*dm.n_t[i] + mC.M_mu[i-1][j]*dm.n_t[i-1] );
                m3_phi = -MagFactor*4.*pi*0.25*( ( mC.M_mu[i][j] + mC.M_mu[i][j-1] )*dm.n_t[i] + ( mC.M_mu[i-1][j] + mC.M_mu[i-1][j-1] )*dm.n_t[i-1] );
                m4_r = -MagFactor*4.*pi*0.5*( mC.C_m[i][j] + mC.C_m[i][j-1] ); //no n_t factor since already included in redshifted T
                m4_th = -MagFactor*4.*pi*0.5*( mC.C_m[i][j] + mC.C_m[i-1][j] ); //no n_t factor since already included in redshifted T
                m4_phi = -MagFactor*4.*pi*0.25*( mC.C_m[i][j] + mC.C_m[i][j-1] + mC.C_m[i-1][j] + mC.C_m[i-1][j-1] ); //no n_t factor since already included in redshifted T

                gradrB = gradr(mC.Bmag,Deltar_n1,i,j)/n_r_av;
                gradthB = gradth(mC.Bmag,dm.r[i]*Deltatheta,i,j);
                gradrT = gradr(T,Deltar_n1,i,j)/n_r_av;
                gradthT = gradth(T,dm.r[i]*Deltatheta,i,j);

                gradrB_2 = gradr(mC.Bmag,Deltar_n1,i,j-1)/n_r_av;
                gradthB_2 = gradth(mC.Bmag,dm.r[i-1]*Deltatheta,i-1,j);
                gradrT_2 = gradr(T,Deltar_n1,i,j-1)/n_r_av;
                gradthT_2 = gradth(T,dm.r[i-1]*Deltatheta,i-1,j);

                if( dm.theta[j]-0.5*Deltatheta > 0. ){
                    J_B_r = dm.n_t[i]*( ( 1./tan(Deltatheta/2.) + 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j] - ( 1./tan(Deltatheta/2.) - 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j-1] )/(2.*dm.r[i]);
                }
                else{
                    //if theta = 0, pi, use
                    J_B_r = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
                }
                if( i < B.shape()[1]-N_GC-1 ){
                    J_B_theta = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*dm.r[i]*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                    J_B_phi = ( dm.n_t[i]*dm.r[i]*B[1][i][j] - dm.n_t[i-1]*dm.r[i-1]*B[1][i-1][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av
                                        + n_t_av*( B[0][i][j-1] - B[0][i][j] )/(0.5*(dm.r[i]+dm.r[i-1])*Deltatheta);
                }
                else{
                    r_imh = dm.r[i] + 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
                    J_B_theta = dm.n_t[i]*( dm.r[i-1]*B[2][i-1][j] - (dm.r[i]+dm.Deltar[i])*B[2][i][j] )/( r_imh*Deltar_n1*n_r_av );
                    J_B_phi = ( dm.n_t[i]*(dm.r[i]+dm.Deltar[i])*B[1][i][j] - dm.n_t[i-1]*dm.r[i-1]*B[1][i-1][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av
                                        + n_t_av*( B[0][i][j-1] - B[0][i][j] )/(0.5*(dm.r[i]+dm.r[i-1])*Deltatheta);
                }
                J[0][i][j] = m1_r*J_B_r + ( m2_r*gradthB + m4_r*gradthT )*mC.rBhat[2][i][j];
                mC.J_B[0][i][j] = J_B_r;
                J[1][i][j] = m1_th*J_B_theta - ( m2_th*gradrB + m3_th*gradrmu_e + m4_th*gradrT )*mC.thBhat[2][i][j];
                mC.J_B[1][i][j] = J_B_theta;
                J[2][i][j] = m1_phi*J_B_phi + m2_phi*( 0.5*(gradrB+gradrB_2)*mC.phBhat[1][i][j] - 0.5*(gradthB+gradthB_2)*mC.phBhat[0][i][j] ) + m3_phi*gradrmu_e*mC.phBhat[1][i][j]
                                                  + m4_phi*( 0.5*(gradrT+gradrT_2)*mC.phBhat[1][i][j] - 0.5*(gradthT+gradthT_2)*mC.phBhat[0][i][j] );
                mC.J_B[2][i][j] = J_B_phi;
            }
            //Computes J_r for theta = pi
            if(process.nbrright < 0){
                size_t j = B.shape()[2]-N_GC;
                m1_r = ( 1. - 4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) );
                J_B_r = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1 - cos(Deltatheta/2.) )*dm.r[i] );
                J[0][i][j] = m1_r*J_B_r;
                mC.J_B[0][i][j] = J_B_r;
            }
        }
    }
    else{
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
            Deltar_n1 = 0.5*( dm.Deltar[i-1]+dm.Deltar[i] );
            n_t_av = 0.5*( dm.n_t[i-1] + dm.n_t[i] );
            n_r_av = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                if( dm.theta[j]-0.5*Deltatheta > 0. ){
                    J[0][i][j] = dm.n_t[i]*( ( 1./tan(Deltatheta/2.) + 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j] - ( 1./tan(Deltatheta/2.) - 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j-1] )/(2.*dm.r[i]);
                }
                else{
                    //if theta = 0, pi, use linear approximation to derivative (vanishing B_phi at these angles != vanishing J_r). This computes Jr for theta = 0
                    J[0][i][j] = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
                }
                if( i < B.shape()[1]-N_GC-1 ){
                    J[1][i][j] = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*dm.r[i]*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                    J[2][i][j] = ( dm.n_t[i]*dm.r[i]*B[1][i][j] - dm.n_t[i-1]*dm.r[i-1]*B[1][i-1][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av
                                        + n_t_av*( B[0][i][j-1] - B[0][i][j] )/(0.5*(dm.r[i]+dm.r[i-1])*Deltatheta);
                }
                else{
                    r_imh = dm.r[i] + 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
                    J[1][i][j] = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*(dm.r[i]+dm.Deltar[i])*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                    J[2][i][j] = ( dm.n_t[i]*(dm.r[i]+dm.Deltar[i])*B[1][i][j] - dm.n_t[i-1]*dm.r[i-1]*B[1][i-1][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av
                                        + n_t_av*( B[0][i][j-1] - B[0][i][j] )/(0.5*(dm.r[i]+dm.r[i-1])*Deltatheta);
                }
            }

            //Computes J_r for theta = pi
            if(process.nbrright < 0){
                size_t j = B.shape()[2]-N_GC;
                J[0][i][j] = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
            }
        }
    }

    exchng2Vector(J, N_GC, process, dm.stridetype_Vec);

    return;
}

/*
    Computes poloidal (r and theta)-components of current times 4*pi/c*(lapse function) components along cell edges in reduced units i.e., 4*pi*L_0/(c*B_0)*exp(nu/2)*J
    Works for both magnetizable B != H and non-magnetizable medium

    Inputs: B: magnetic field in reduced units at cell centers
            T: redshifted temperature in reduced units
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process objects containing information about the simulation domain and the current process
    Output: J: redshifted current density in reduced units
*/
void Compute_J_Poloidal(VectorField & B, VectorField & J, ScalarField & T, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process)
{
    double Deltatheta = dm.Deltatheta;
    double m1_r, m2_r, m4_r; //magnetization-dependent prefactors for the curl(B), grad(B)xBhat and grad(T)xBhat terms respectively evaluated at location of J_r
    double m1_th, m2_th, m3_th, m4_th; //magnetization-dependent prefactors for the curl(B), grad(B)xBhat, grad(mu_e)xBhat and grad(T)xBhat terms respectively evaluated at location of J_th
    double gradrB, gradthB, gradrmu_e, gradrT, gradthT; //components of gradients of Bmag, mu_e and T

    double r_imh, Deltar_n1, n_r_av;

    if(mC.quantization == true){
        double J_B_r, J_B_theta; //curl(B) components
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
            Deltar_n1 = 0.5*( dm.Deltar[i-1] + dm.Deltar[i] );
            n_r_av = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            gradrmu_e = (mC.mu_e[i]-mC.mu_e[i-1])/Deltar_n1/n_r_av;

            double MagFactor = 1.;
            for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                if(mC.PolarAxisQuantization == true && ( dm.theta[j] - dm.Deltatheta < 0 || dm.theta[j] + dm.Deltatheta > pi )){
                    MagFactor = 0.;
                }
                else MagFactor = 1.;

                m1_r = ( 1. - MagFactor*4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) );
                m1_th = ( 1. - MagFactor*4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i-1][j]/mC.Bmag[i-1][j] )  );
                m2_r = -MagFactor*4.*pi*( 0.5*( mC.chi[i][j] + mC.chi[i][j-1] ) - 0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) )*dm.n_t[i];
                m2_th = -MagFactor*4.*pi*( 0.5*( mC.chi[i][j]*dm.n_t[i] + mC.chi[i-1][j]*dm.n_t[i-1] ) - 0.5*( mC.M[i][j]/mC.Bmag[i][j]*dm.n_t[i] + mC.M[i-1][j]/mC.Bmag[i-1][j]*dm.n_t[i-1] ) );;
                m3_th = -MagFactor*4.*pi*0.5*( mC.M_mu[i][j]*dm.n_t[i] + mC.M_mu[i-1][j]*dm.n_t[i-1]);
                m4_r = -MagFactor*4.*pi*0.5*( mC.C_m[i][j] + mC.C_m[i][j-1] ); //no n_t factor since already included in redshifted T
                m4_th = -MagFactor*4.*pi*0.5*( mC.C_m[i][j] + mC.C_m[i-1][j] ); //no n_t factor since already included in redshifted T

                gradrB = gradr(mC.Bmag,Deltar_n1,i,j)/n_r_av;
                gradthB = gradth(mC.Bmag,dm.r[i]*Deltatheta,i,j);
                gradrT = gradr(T,Deltar_n1,i,j)/n_r_av;
                gradthT = gradth(T,dm.r[i]*Deltatheta,i,j);

                if( dm.theta[j] - 0.5*Deltatheta > 0 ){
                    J_B_r = dm.n_t[i]*( ( 1./tan(Deltatheta/2.) + 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j] - ( 1./tan(Deltatheta/2.) - 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j-1] )/(2.*dm.r[i]);
                    // J_B_r = dm.n_t[i]*( 1./tan(Deltatheta/2.)*(B[2][i][j] - B[2][i][j-1]) + 1./tan(dm.theta[j]-0.5*Deltatheta)*(B[2][i][j] + B[2][i][j-1]) )/(2.*dm.r[i]);
                }
                else{
                    //if theta = 0, pi, set Jacobian-dependent term equal to zero, since Bphi vanishes at this boundary. Avoids 1/tan(0) error.
                    J_B_r = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
                }
                if( i < B.shape()[1]-N_GC-1 ){
                    J_B_theta = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*dm.r[i]*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                }
                else{
                    r_imh = dm.r[i] + 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
                    J_B_theta = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*(dm.r[i]+dm.Deltar[i])*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                }
                J[0][i][j] = m1_r*J_B_r + ( m2_r*gradthB + m4_r*gradthT )*mC.rBhat[2][i][j];
                mC.J_B[0][i][j] = J_B_r;
                J[1][i][j] = m1_th*J_B_theta - ( m2_th*gradrB + m3_th*gradrmu_e + m4_th*gradrT )*mC.thBhat[2][i][j];
                mC.J_B[1][i][j] = J_B_theta;
            }
            if(process.nbrright < 0){
                size_t j = B.shape()[2]-N_GC;
                m1_r = ( 1. - 4.*pi*0.5*( mC.M[i][j]/mC.Bmag[i][j] + mC.M[i][j-1]/mC.Bmag[i][j-1] ) );
                J_B_r = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1 - cos(Deltatheta/2.) )*dm.r[i] );
                J[0][i][j] = m1_r*J_B_r;
                mC.J_B[0][i][j] = J_B_r;
            }
        }
    }
    else{
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
            Deltar_n1 = 0.5*( dm.Deltar[i-1] + dm.Deltar[i] );
            n_r_av = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                if( dm.theta[j] - 0.5*Deltatheta > 0 ){
                    J[0][i][j] = dm.n_t[i]*( ( 1./tan(Deltatheta/2.) + 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j] - ( 1./tan(Deltatheta/2.) - 1./tan(dm.theta[j]-0.5*Deltatheta) )*B[2][i][j-1] )/(2.*dm.r[i]);
                    // J[0][i][j] = dm.n_t[i]*( 1./tan(Deltatheta/2.)*(B[2][i][j] - B[2][i][j-1]) + 1./tan(dm.theta[j]-0.5*Deltatheta)*(B[2][i][j] + B[2][i][j-1]) )/(2.*dm.r[i]);
                }
                else{
                    //if theta = 0, pi, use linear approximation to derivative (vanishing B_phi at these angles != vanishing J_r). This computes Jr for theta = 0
                    J[0][i][j] = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
                }
                if( i < B.shape()[1]-N_GC-1 ){
                    J[1][i][j] = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*dm.r[i]*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                }
                else{
                    r_imh = dm.r[i] + 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
                    J[1][i][j] = ( dm.n_t[i-1]*dm.r[i-1]*B[2][i-1][j] - dm.n_t[i]*(dm.r[i]+dm.Deltar[i])*B[2][i][j] )/( r_imh*Deltar_n1 + ( dm.Deltar[i]*dm.Deltar[i] - dm.Deltar[i-1]*dm.Deltar[i-1] )/8. )/n_r_av;
                }
            }
            //Computes J_r for theta = pi
            if(process.nbrright < 0){
                size_t j = B.shape()[2]-N_GC;
                J[0][i][j] = dm.n_t[i]*(B[2][i][j] - B[2][i][j-1])*sin(Deltatheta/2.)/( 2.*( 1. - cos(Deltatheta/2.) )*dm.r[i] );
            }
        }
    }

    exchng2Vector(J, N_GC, process, dm.stridetype_Vec);
    if( mC.quantization == true ){
        exchng2Vector(mC.J_B, N_GC, process, dm.stridetype_Vec);
    }

    return;
}

/*
    Computes the update for the toroidal field component. Works for both magnetizable B != H and nonmagnetizable medium.
    MC slope limiter-based reconstruction performed using cell centroids in spherical coordinates from Mignone J. Comp. Phys. 270, 784 (2014)

    Uses the Godunov scheme with MC limiter
    Inputs: B, J, T: magnetic field, redshifted current density and redshifted temperature in reduced units
            tC: TransCoeffs object containing Ohmic and Hall diffusivities eta_O and eta_H in reduced units
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process objects containing information about the simulation domain and the current process
    Output: Qphi: update for Bphi
*/
void B_torEvolve(ScalarField & Qphi, ScalarField & q_SH, VectorField & B, VectorField & J, ScalarField & T, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process)
{
    double Deltatheta = dm.Deltatheta;
    double Fr, Ftheta; //flux difference function in r and theta directions and cross term
    double r_imh, r_iph, deltar_cent;
    double Deltar_1, Deltar_2, Deltar_3, Deltar_4; //twice the cell-centroid spacings in radial direction. Include corrections for variable cell spacing.
    double Deltatheta_1, Deltatheta_2, Deltatheta_3, Deltatheta_4; //twice the cell-centroid spacings in theta direction.
    double beta_r_m1, beta_r, beta_r_p1, beta_theta_m1, beta_theta, beta_theta_p1;
    double Delta_xir_im1_m1, Delta_xir_im1_p1, Delta_xir_i_m1, Delta_xir_i_p1, Delta_xir_ip1_m1, Delta_xir_ip1_p1; //half cell-centroid spacings in radial direction
    double Delta_xitheta_jm1_m1, Delta_xitheta_jm1_p1, Delta_xitheta_j_m1, Delta_xitheta_j_p1, Delta_xitheta_jp1_m1, Delta_xitheta_jp1_p1; //half cell-centroid spacings in theta direction
    double sigma_im1, sigma_i, sigma_ip1, sigma_jm1, sigma_j, sigma_jp1;
    double B_im1, B_im, B_ip, B_ip1;
    double B_jm1, B_jm, B_jp, B_jp1;
    double s_r_f, s_r_b, s_theta_f, s_theta_b;
    double s_r_f_p, s_r_f_m, s_r_b_p, s_r_b_m, s_theta_f_p, s_theta_f_m, s_theta_b_p, s_theta_b_m;
    double eta_O_th_ip1, eta_O_th_im1, eta_O_r_jp1, eta_O_r_jm1;
    double eta_O_perp_th_ip1, eta_O_perp_th_im1, eta_O_perp_r_jp1, eta_O_perp_r_jm1;
    double eta_O_delta_th_ip1, eta_O_delta_th_im1, eta_O_delta_r_jp1, eta_O_delta_r_jm1;
    double vphi_r_b, vphi_r_f, vphi_theta_b, vphi_theta_f;
    double Delta2, eta_hyp, EdgeLimit;
    double lapBphi_im1, lapBphi, lapBphi_ip1, lapBphi_jm1, lapBphi_jp1;

    double E_r_HB_jmh, E_r_HB_jph, E_theta_HB_imh, E_theta_HB_iph; //Burgers part of Hall electric field
    double E_r_HB_mag_jmh = 0., E_r_HB_mag_jph = 0., E_theta_HB_mag_imh = 0., E_theta_HB_mag_iph = 0.; // "cross" part of Hall electric field
    double theta_m, theta_p;
    double Jr_av_phi_lb, Jr_av_phi_lt, Jr_av_phi_rb, Jr_av_phi_rt, Jtheta_av_phi_lb, Jtheta_av_phi_lt, Jtheta_av_phi_rb, Jtheta_av_phi_rt;
    double eta_H_corner_lb, eta_H_corner_lt, eta_H_corner_rb, eta_H_corner_rt;
    double vr_av_phi_lb, vr_av_phi_lt, vr_av_phi_rb, vr_av_phi_rt, vtheta_av_phi_lb, vtheta_av_phi_lt, vtheta_av_phi_rb, vtheta_av_phi_rt, vphi_phi_lb, vphi_phi_lt, vphi_phi_rb, vphi_phi_rt;
    double Br_up_phi_lb, Br_up_phi_lt, Br_up_phi_rb, Br_up_phi_rt, Btheta_up_phi_lb, Btheta_up_phi_lt, Btheta_up_phi_rb, Btheta_up_phi_rt;

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        r_imh = dm.r[i] - 0.5*dm.Deltar[i]; //r_{i-1/2}, location of cell boundary between cells i and i-1
        r_iph = dm.r[i] + 0.5*dm.Deltar[i]; //r_{i+1/2}, location of cell boundary between cells i and i+1
        deltar_cent = 0.25*dm.Deltar[i-1] + 0.5*dm.Deltar[i] + 0.25*dm.Deltar[i+1];

        Deltar_1 = dm.Deltar[i-1] + dm.Deltar[i-2] + 1./3.*( dm.r[i-1]*dm.Deltar[i-1]*dm.Deltar[i-1]/(dm.r[i-1]*dm.r[i-1] + dm.Deltar[i-1]*dm.Deltar[i-1]/12.)
                                                             - dm.r[i-2]*dm.Deltar[i-2]*dm.Deltar[i-2]/(dm.r[i-2]*dm.r[i-2] + dm.Deltar[i-2]*dm.Deltar[i-2]/12.) );
        Deltar_2 = dm.Deltar[i] + dm.Deltar[i-1] + 1./3.*( dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i] + dm.Deltar[i]*dm.Deltar[i]/12.)
                                                           - dm.r[i-1]*dm.Deltar[i-1]*dm.Deltar[i-1]/(dm.r[i-1]*dm.r[i-1] + dm.Deltar[i-1]*dm.Deltar[i-1]/12.) );
        Deltar_3 = dm.Deltar[i+1] + dm.Deltar[i] + 1./3.*( dm.r[i+1]*dm.Deltar[i+1]*dm.Deltar[i+1]/(dm.r[i+1]*dm.r[i+1] + dm.Deltar[i+1]*dm.Deltar[i+1]/12.)
                                                           - dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i] + dm.Deltar[i]*dm.Deltar[i]/12.) );
        Deltar_4 = dm.Deltar[i+2] + dm.Deltar[i+1] + 1./3.*( dm.r[i+2]*dm.Deltar[i+2]*dm.Deltar[i+2]/(dm.r[i+2]*dm.r[i+2] + dm.Deltar[i+2]*dm.Deltar[i+2]/12.)
                                                             - dm.r[i+1]*dm.Deltar[i+1]*dm.Deltar[i+1]/(dm.r[i+1]*dm.r[i+1] + dm.Deltar[i+1]*dm.Deltar[i+1]/12.) );

        Delta_xir_im1_m1 = 0.5*dm.Deltar[i-1]*( 1. + 4.*dm.r[i-1]*dm.Deltar[i-1]/(12.*dm.r[i-1]*dm.r[i-1]+dm.Deltar[i-1]*dm.Deltar[i-1]) );
        Delta_xir_im1_p1 = 0.5*dm.Deltar[i-1]*( 1. - 4.*dm.r[i-1]*dm.Deltar[i-1]/(12.*dm.r[i-1]*dm.r[i-1]+dm.Deltar[i-1]*dm.Deltar[i-1]) );
        Delta_xir_i_m1 = 0.5*dm.Deltar[i]*( 1. + 4.*dm.r[i]*dm.Deltar[i]/(12.*dm.r[i]*dm.r[i]+dm.Deltar[i]*dm.Deltar[i]) );
        Delta_xir_i_p1 = 0.5*dm.Deltar[i]*( 1. - 4.*dm.r[i]*dm.Deltar[i]/(12.*dm.r[i]*dm.r[i]+dm.Deltar[i]*dm.Deltar[i]) );
        Delta_xir_ip1_m1 = 0.5*dm.Deltar[i+1]*( 1. + 4.*dm.r[i+1]*dm.Deltar[i+1]/(12.*dm.r[i+1]*dm.r[i+1]+dm.Deltar[i+1]*dm.Deltar[i+1]) );
        Delta_xir_ip1_p1 = 0.5*dm.Deltar[i+1]*( 1. - 4.*dm.r[i+1]*dm.Deltar[i+1]/(12.*dm.r[i+1]*dm.r[i+1]+dm.Deltar[i+1]*dm.Deltar[i+1]) );

        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++) {
            Fr = 0.;
            Ftheta = 0.;
            Delta2 = 1./( 1./(dm.Deltar[i]*dm.Deltar[i]) + 1./(dm.r[i]*dm.r[i]*Deltatheta*Deltatheta) );
            theta_m = dm.theta[j] - 0.5*dm.Deltatheta;
            theta_p = dm.theta[j] + 0.5*dm.Deltatheta;

            Deltatheta_1 = 2.*( Deltatheta + ( 1./tan(dm.theta[j-1]) - 1./tan(dm.theta[j-2]) )*dm.Deltatheta_dfactor);
            Deltatheta_2 = 2.*( Deltatheta + ( 1./tan(dm.theta[j]) - 1./tan(dm.theta[j-1]) )*dm.Deltatheta_dfactor );
            Deltatheta_3 = 2.*( Deltatheta + ( 1./tan(dm.theta[j+1]) - 1./tan(dm.theta[j]) )*dm.Deltatheta_dfactor );
            Deltatheta_4 = 2.*( Deltatheta + ( 1./tan(dm.theta[j+2]) - 1./tan(dm.theta[j+1]) )*dm.Deltatheta_dfactor );

            Delta_xitheta_jm1_m1 = 0.5*Deltatheta + 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
            Delta_xitheta_jm1_p1 = 0.5*Deltatheta - 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
            Delta_xitheta_j_m1 = 0.5*Deltatheta + 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
            Delta_xitheta_j_p1 = 0.5*Deltatheta - 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
            Delta_xitheta_jp1_m1 = 0.5*Deltatheta + 1./tan(dm.theta[j+1])*dm.Deltatheta_dfactor;
            Delta_xitheta_jp1_p1 = 0.5*Deltatheta - 1./tan(dm.theta[j+1])*dm.Deltatheta_dfactor;

            //Compute beta_r and beta_theta, parts of the toroidal B field advection velocities
            beta_r_p1 = dm.n_t[i+1]*dm.n_t[i+1]*dm.r[i+1]*( ( tC.psi_H[i+1][j+1] - tC.psi_H[i+1][j-1] )/( 2.*Deltatheta ) - 2.*tC.psi_H[i+1][j]/tan(dm.theta[j]) );
            beta_r = dm.n_t[i]*dm.n_t[i]*dm.r[i]*( ( tC.psi_H[i][j+1] - tC.psi_H[i][j-1] )/( 2.*Deltatheta ) - 2.*tC.psi_H[i][j]/tan(dm.theta[j]) );
            beta_r_m1 = dm.n_t[i-1]*dm.n_t[i-1]*dm.r[i-1]*( ( tC.psi_H[i-1][j+1] - tC.psi_H[i-1][j-1] )/( 2.*Deltatheta ) - 2.*tC.psi_H[i-1][j]/tan(dm.theta[j]) );

            beta_theta_p1 = -dm.n_t[i]*dm.n_t[i]/dm.n_r[i]*dm.r[i]*dm.r[i]*( tC.psi_H[i+1][j+1] - tC.psi_H[i-1][j+1])/( 2.*deltar_cent );
            beta_theta = -dm.n_t[i]*dm.n_t[i]/dm.n_r[i]*dm.r[i]*dm.r[i]*( tC.psi_H[i+1][j] - tC.psi_H[i-1][j] )/( 2.*deltar_cent ) ;
            beta_theta_m1 = -dm.n_t[i]*dm.n_t[i]/dm.n_r[i]*dm.r[i]*dm.r[i]*( tC.psi_H[i+1][j-1] - tC.psi_H[i-1][j-1] )/( 2.*deltar_cent );

            //Reconstruct values of toroidal field at the cell edges
            sigma_im1 = MC( ( B[2][i][j] - B[2][i-1][j] )/Deltar_2 + ( B[2][i-1][j] - B[2][i-2][j] )/Deltar_1,
                                ( B[2][i][j] - B[2][i-1][j] )/Delta_xir_im1_p1,
                                ( B[2][i-1][j] - B[2][i-2][j] )/Delta_xir_im1_m1 );
            sigma_i = MC( ( B[2][i+1][j] - B[2][i][j] )/Deltar_3 + ( B[2][i][j] - B[2][i-1][j] )/Deltar_2,
                                ( B[2][i+1][j] - B[2][i][j] )/Delta_xir_i_p1,
                                ( B[2][i][j] - B[2][i-1][j] )/Delta_xir_i_m1 );
            sigma_ip1 = MC( ( B[2][i+2][j] - B[2][i+1][j] )/Deltar_4 + ( B[2][i+1][j] - B[2][i][j] )/Deltar_3,
                                ( B[2][i+2][j] - B[2][i+1][j] )/Delta_xir_ip1_p1,
                                ( B[2][i+1][j] - B[2][i][j] )/Delta_xir_ip1_m1 );

            B_im1 = B[2][i-1][j] + Delta_xir_im1_p1*sigma_im1;
            B_im = B[2][i][j] - Delta_xir_i_m1*sigma_i;
            B_ip = B[2][i][j] + Delta_xir_i_p1*sigma_i;
            B_ip1 = B[2][i+1][j] - Delta_xir_ip1_m1*sigma_ip1;

            sigma_jm1 = MC( ( B[2][i][j] - B[2][i][j-1])/Deltatheta_2 + (B[2][i][j-1] - B[2][i][j-2])/Deltatheta_1,
                            ( B[2][i][j] - B[2][i][j-1] )/Delta_xitheta_jm1_p1,
                            ( B[2][i][j-1] - B[2][i][j-2] )/Delta_xitheta_jm1_m1 );
            sigma_j = MC( ( B[2][i][j+1] - B[2][i][j] )/Deltatheta_3 + ( B[2][i][j] - B[2][i][j-1] )/Deltatheta_2,
                          ( B[2][i][j+1] - B[2][i][j] )/Delta_xitheta_j_p1,
                          ( B[2][i][j] - B[2][i][j-1] )/Delta_xitheta_j_m1 );
            sigma_jp1 = MC( ( B[2][i][j+2] - B[2][i][j+1] )/Deltatheta_4 + ( B[2][i][j+1] - B[2][i][j] )/Deltatheta_3,
                            ( B[2][i][j+2] - B[2][i][j+1] )/Delta_xitheta_jp1_p1,
                            ( B[2][i][j+1] - B[2][i][j] )/Delta_xitheta_jp1_m1 );

            B_jm1 = B[2][i][j-1] + Delta_xitheta_jm1_p1*sigma_jm1;
            B_jm = B[2][i][j] - Delta_xitheta_j_m1*sigma_j;
            B_jp = B[2][i][j] + Delta_xitheta_j_p1*sigma_j;
            B_jp1 = B[2][i][j+1] - Delta_xitheta_jp1_m1*sigma_jp1;

            //Local propagation speeds at cell boundaries according to Kurganov, Noelle and Petrova, SIAM J. Sci. Comput. 23, 707-740 (2001) Eq. (3.14-3.15)
            // s_r_b_p = std::max( {beta_r_m1*B_im1, beta_r*B_im, 0.} );
            // s_r_b_m = std::min( {beta_r_m1*B_im1, beta_r*B_im, 0.} );
            // s_r_f_p = std::max( {beta_r*B_ip, beta_r_p1*B_ip1, 0.} );
            // s_r_f_m = std::min( {beta_r*B_ip, beta_r_p1*B_ip1, 0.} );
            //
            // s_theta_b_p = std::max( {beta_theta_m1*B_jm1, beta_theta*B_jm, 0.} );
            // s_theta_b_m = std::min( {beta_theta_m1*B_jm1, beta_theta*B_jm, 0.} );
            // s_theta_f_p = std::max( {beta_theta*B_jp, beta_theta_p1*B_jp1, 0.} );
            // s_theta_f_m = std::min( {beta_theta*B_jp, beta_theta_p1*B_jp1, 0.} );
            //
            // E_theta_HB_imh = ( s_r_b_p*beta_r_m1*B_im1*B_im1/2. - s_r_b_m*beta_r*B_im*B_im/2. + s_r_b_p*s_r_b_m*( B_im - B_im1 ) )/(s_r_b_p - s_r_b_m);
            // E_theta_HB_iph = ( s_r_f_p*beta_r*B_ip*B_ip/2. - s_r_f_m*beta_r_p1*B_ip1*B_ip1/2. + s_r_f_p*s_r_f_m*( B_ip1 - B_ip ) )/(s_r_f_p - s_r_f_m);
            // Fr += r_iph*E_theta_HB_iph - r_imh*E_theta_HB_imh;
            //
            // E_r_HB_jmh = -( s_theta_b_p*beta_theta_m1*B_jm1*B_jm1/2. - s_theta_b_m*beta_theta*B_jm*B_jm/2. + s_theta_b_p*s_theta_b_m*( B_jm - B_jm1 ) )/(s_theta_b_p - s_theta_b_m);
            // E_r_HB_jph = -( s_theta_f_p*beta_theta*B_jp*B_jp/2. - s_theta_f_m*beta_theta_p1*B_jp1*B_jp1/2. + s_theta_f_p*s_theta_f_m*( B_jp1 - B_jp ) )/(s_theta_f_p - s_theta_f_m);
            // Ftheta += -( E_r_HB_jph - E_r_HB_jmh);

            //KT central scheme
            // s_r_b = std::max( abs(beta_r_m1*B_im1), abs(beta_r*B_im) );
            // s_r_f = std::max( abs(beta_r*B_ip), abs(beta_r_p1*B_ip1) );
            //
            // s_theta_b = std::max( abs(beta_theta_m1*B_jm1), abs(beta_theta*B_jm) );
            // s_theta_f = std::max( abs(beta_theta*B_jp), abs(beta_theta_p1*B_jp1) );
            //
            // E_theta_HB_imh = 1./2.*( beta_r*B_im*B_im/2. + beta_r_m1*B_im1*B_im1/2. - s_r_b*( B_im - B_im1 ) );
            // E_theta_HB_iph = 1./2.*( beta_r_p1*B_ip1*B_ip1/2. + beta_r*B_ip*B_ip/2. - s_r_f*( B_ip1 - B_ip ) );
            // Fr += r_iph*E_theta_HB_iph - r_imh*E_theta_HB_imh;
            //
            // E_r_HB_jmh = -1./2.*( beta_theta*B_jm*B_jm/2. + beta_theta_m1*B_jm1*B_jm1/2. - s_theta_b*( B_jm - B_jm1 ) );
            // E_r_HB_jph = -1./2.*( beta_theta_p1*B_jp1*B_jp1/2. + beta_theta*B_jp*B_jp/2. - s_theta_f*( B_jp1 - B_jp ) );
            // Ftheta += -( E_r_HB_jph - E_r_HB_jmh);

            //Godunov method with MC reconstruction
            s_r_b_p = beta_r_m1*B_im1;
            s_r_b_m = beta_r*B_im;
            s_r_f_p = beta_r*B_ip;
            s_r_f_m = beta_r_p1*B_ip1;

            s_theta_b_p = beta_theta_m1*B_jm1;
            s_theta_b_m = beta_theta*B_jm;
            s_theta_f_p = beta_theta*B_jp;
            s_theta_f_m = beta_theta_p1*B_jp1;

            E_theta_HB_imh = Godunov(s_r_b_m, s_r_b_p, beta_r_m1, B_im1, beta_r, B_im);
            E_theta_HB_iph = Godunov(s_r_f_m, s_r_f_p, beta_r, B_ip, beta_r_p1, B_ip1);
            Fr += r_iph*E_theta_HB_iph - r_imh*E_theta_HB_imh;

            E_r_HB_jmh = -Godunov(s_theta_b_m, s_theta_b_p, beta_theta_m1, B_jm1, beta_theta, B_jm);
            E_r_HB_jph = -Godunov(s_theta_f_m, s_theta_f_p, beta_theta, B_jp, beta_theta_p1, B_jp1);
            Ftheta += -( E_r_HB_jph - E_r_HB_jmh);

            double Bphi_phi_lb, Bphi_phi_lt, Bphi_phi_rb, Bphi_phi_rt; //B_phi interpolated to the cell corners (location of J_phi)
            //Magnetization-dependent parts of B_phi^2 terms. Separated from B=H part of the B_phi evolution equation.
            if(mC.quantization == true){
                eta_H_corner_lb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j] + 1./tC.eta_H[i-1][j-1] ) );
                eta_H_corner_lt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j] + 1./tC.eta_H[i+1][j-1] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) );
                eta_H_corner_rb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j+1] + 1./tC.eta_H[i-1][j] ) );
                eta_H_corner_rt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j+1] + 1./tC.eta_H[i+1][j] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) );

                if( process.nbrleft < 0 && j == N_GC ){
                    Jtheta_av_phi_lb = 0.;
                    Jtheta_av_phi_lt = 0.;
                    Jtheta_av_phi_rb = 0.5*( sin(dm.theta[j+1])*( J[1][i][j+1] - mC.J_B[1][i][j+1] ) + sin(dm.theta[j])*( J[1][i][j] - mC.J_B[1][i][j] ) )/sin(theta_p);
                    Jtheta_av_phi_rt = 0.5*( sin(dm.theta[j+1])*( J[1][i+1][j+1] - mC.J_B[1][i+1][j+1] ) + sin(dm.theta[j])*( J[1][i+1][j] - mC.J_B[1][i+1][j] ) )/sin(theta_p);

                    Bphi_phi_lb = 0.;
                    Bphi_phi_lt = 0.;
                    Bphi_phi_rb = 0.25*( dm.r[i]*( sin(dm.theta[j+1])*B[2][i][j+1] + sin(dm.theta[j])*B[2][i][j] )
                                            + dm.r[i-1]*( sin(dm.theta[j+1])*B[2][i-1][j+1] + sin(dm.theta[j])*B[2][i-1][j] ) )/(r_imh*sin(theta_p));
                    Bphi_phi_rt = 0.25*( dm.r[i+1]*( sin(dm.theta[j+1])*B[2][i+1][j+1] + sin(dm.theta[j])*B[2][i+1][j] )
                                            + dm.r[i]*( sin(dm.theta[j+1])*B[2][i][j+1] + sin(dm.theta[j])*B[2][i][j] ) )/(r_iph*sin(theta_p));
                }
                else if( process.nbrright < 0 && j == B.shape()[2]-N_GC ) {
                    Jtheta_av_phi_lb = 0.5*( sin(dm.theta[j])*( J[1][i][j] - mC.J_B[1][i][j] ) + sin(dm.theta[j-1])*( J[1][i][j-1] - mC.J_B[1][i][j-1] ) )/sin(theta_m);
                    Jtheta_av_phi_lt = 0.5*( sin(dm.theta[j])*( J[1][i+1][j] - mC.J_B[1][i+1][j] ) + sin(dm.theta[j-1])*( J[1][i+1][j-1] - mC.J_B[1][i+1][j-1] ) )/sin(theta_m);
                    Jtheta_av_phi_rb = 0.;
                    Jtheta_av_phi_rt = 0.;

                    Bphi_phi_lb = 0.25*( dm.r[i]*( sin(dm.theta[j])*B[2][i][j] + sin(dm.theta[j-1])*B[2][i][j-1] )
                                            + dm.r[i-1]*( sin(dm.theta[j])*B[2][i-1][j] + sin(dm.theta[j-1])*B[2][i-1][j-1] ) )/(r_imh*sin(theta_m));
                    Bphi_phi_lt = 0.25*( dm.r[i+1]*( sin(dm.theta[j])*B[2][i+1][j] + sin(dm.theta[j-1])*B[2][i+1][j-1] )
                                            + dm.r[i]*( sin(dm.theta[j])*B[2][i][j] + sin(dm.theta[j-1])*B[2][i][j-1] ) )/(r_iph*sin(theta_m));
                    Bphi_phi_rb = 0.;
                    Bphi_phi_rt = 0.;
                }
                else{
                    Jtheta_av_phi_lb = 0.5*( sin(dm.theta[j])*( J[1][i][j] - mC.J_B[1][i][j] ) + sin(dm.theta[j-1])*( J[1][i][j-1] - mC.J_B[1][i][j-1] ) )/sin(theta_m);
                    Jtheta_av_phi_lt = 0.5*( sin(dm.theta[j])*( J[1][i+1][j] - mC.J_B[1][i+1][j] ) + sin(dm.theta[j-1])*( J[1][i+1][j-1] - mC.J_B[1][i+1][j-1] ) )/sin(theta_m);
                    Jtheta_av_phi_rb = 0.5*( sin(dm.theta[j+1])*( J[1][i][j+1] - mC.J_B[1][i][j+1] ) + sin(dm.theta[j])*( J[1][i][j] - mC.J_B[1][i][j] ) )/sin(theta_p);
                    Jtheta_av_phi_rt = 0.5*( sin(dm.theta[j+1])*( J[1][i+1][j+1] - mC.J_B[1][i+1][j+1] ) + sin(dm.theta[j])*( J[1][i+1][j] - mC.J_B[1][i+1][j] ) )/sin(theta_p);

                    Bphi_phi_lb = 0.25*( dm.r[i]*( sin(dm.theta[j])*B[2][i][j] + sin(dm.theta[j-1])*B[2][i][j-1] )
                                            + dm.r[i-1]*( sin(dm.theta[j])*B[2][i-1][j] + sin(dm.theta[j-1])*B[2][i-1][j-1] ) )/(r_imh*sin(theta_m));
                    Bphi_phi_lt = 0.25*( dm.r[i+1]*( sin(dm.theta[j])*B[2][i+1][j] + sin(dm.theta[j-1])*B[2][i+1][j-1] )
                                            + dm.r[i]*( sin(dm.theta[j])*B[2][i][j] + sin(dm.theta[j-1])*B[2][i][j-1] ) )/(r_iph*sin(theta_m));
                    Bphi_phi_rb = 0.25*( dm.r[i]*( sin(dm.theta[j+1])*B[2][i][j+1] + sin(dm.theta[j])*B[2][i][j] )
                                            + dm.r[i-1]*( sin(dm.theta[j+1])*B[2][i-1][j+1] + sin(dm.theta[j])*B[2][i-1][j] ) )/(r_imh*sin(theta_p));
                    Bphi_phi_rt = 0.25*( dm.r[i+1]*( sin(dm.theta[j+1])*B[2][i+1][j+1] + sin(dm.theta[j])*B[2][i+1][j] )
                                            + dm.r[i]*( sin(dm.theta[j+1])*B[2][i][j+1] + sin(dm.theta[j])*B[2][i][j] ) )/(r_iph*sin(theta_p));
                }

                vtheta_av_phi_lb = -eta_H_corner_lb*Jtheta_av_phi_lb;
                vtheta_av_phi_lt = -eta_H_corner_lt*Jtheta_av_phi_lt;
                vtheta_av_phi_rb = -eta_H_corner_rb*Jtheta_av_phi_rb;
                vtheta_av_phi_rt = -eta_H_corner_rt*Jtheta_av_phi_rt;

                Jr_av_phi_lb = 0.5*( dm.r[i]*( J[0][i][j] - mC.J_B[0][i][j] )+ dm.r[i-1]*( J[0][i-1][j] - mC.J_B[0][i-1][j] ) )/r_imh;
                Jr_av_phi_lt = 0.5*( dm.r[i+1]*( J[0][i+1][j] - mC.J_B[0][i+1][j] )+ dm.r[i]*( J[0][i][j] - mC.J_B[0][i][j] ) )/r_iph;
                Jr_av_phi_rb = 0.5*( dm.r[i]*( J[0][i][j+1] - mC.J_B[0][i][j+1] )+ dm.r[i-1]*( J[0][i-1][j+1] - mC.J_B[0][i-1][j+1] ) )/r_imh;
                Jr_av_phi_rt = 0.5*( dm.r[i+1]*( J[0][i+1][j+1] - mC.J_B[0][i+1][j+1] )+ dm.r[i]*( J[0][i][j+1] - mC.J_B[0][i][j+1] ) )/r_iph;

                vr_av_phi_lb = -eta_H_corner_lb*Jr_av_phi_lb;
                vr_av_phi_lt = -eta_H_corner_lt*Jr_av_phi_lt;
                vr_av_phi_rb = -eta_H_corner_rb*Jr_av_phi_rb;
                vr_av_phi_rt = -eta_H_corner_rt*Jr_av_phi_rt;

                E_r_HB_mag_jmh = 0.5*( -vtheta_av_phi_lb*Bphi_phi_lb + -vtheta_av_phi_lt*Bphi_phi_lt );
                E_r_HB_mag_jph = 0.5*( -vtheta_av_phi_rb*Bphi_phi_rb + -vtheta_av_phi_rt*Bphi_phi_rt );
                E_theta_HB_mag_imh = -0.5*( -vr_av_phi_lb*Bphi_phi_lb + -vr_av_phi_rb*Bphi_phi_rb );
                E_theta_HB_mag_iph = -0.5*( -vr_av_phi_lt*Bphi_phi_lt + -vr_av_phi_lt*Bphi_phi_lt );

                Fr += r_iph*E_theta_HB_mag_iph - r_imh*E_theta_HB_mag_imh;
                Ftheta += -( E_r_HB_mag_jph - E_r_HB_mag_jmh );
            }

            //Jtheta interpolation onto cell corners that lie on theta=0,pi boundaries
            if( process.nbrleft < 0 && j == N_GC ){
                Jtheta_av_phi_lb = 0.;
                Jtheta_av_phi_lt = 0.;
                Jtheta_av_phi_rb = 0.5*( sin(dm.theta[j+1])*J[1][i][j+1] + sin(dm.theta[j])*J[1][i][j] )/sin(theta_p);
                Jtheta_av_phi_rt = 0.5*( sin(dm.theta[j+1])*J[1][i+1][j+1] + sin(dm.theta[j])*J[1][i+1][j] )/sin(theta_p);
            }
            else if( process.nbrright < 0 && j == B.shape()[2]-N_GC ){
                Jtheta_av_phi_lb = 0.5*( sin(dm.theta[j])*J[1][i][j] + sin(dm.theta[j-1])*J[1][i][j-1] )/sin(theta_m);
                Jtheta_av_phi_lt = 0.5*( sin(dm.theta[j])*J[1][i+1][j] + sin(dm.theta[j-1])*J[1][i+1][j-1] )/sin(theta_m);
                Jtheta_av_phi_rb = 0.;
                Jtheta_av_phi_rt = 0.;
            }
            else{
                Jtheta_av_phi_lb = 0.5*( sin(dm.theta[j])*J[1][i][j] + sin(dm.theta[j-1])*J[1][i][j-1] )/sin(theta_m);
                Jtheta_av_phi_lt = 0.5*( sin(dm.theta[j])*J[1][i+1][j] + sin(dm.theta[j-1])*J[1][i+1][j-1] )/sin(theta_m);
                Jtheta_av_phi_rb = 0.5*( sin(dm.theta[j+1])*J[1][i][j+1] + sin(dm.theta[j])*J[1][i][j] )/sin(theta_p);
                Jtheta_av_phi_rt = 0.5*( sin(dm.theta[j+1])*J[1][i+1][j+1] + sin(dm.theta[j])*J[1][i+1][j] )/sin(theta_p);
            }

            //Quantities at location of J_phi: lb = {i-1/2, j-1/2}, lt = {i+1/2, j-1/2}, rb = {i-1/2, j+1/2}, rt = {i+1/2, j+1/2},
            eta_H_corner_lb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j] + 1./tC.eta_H[i-1][j-1] ) );
            eta_H_corner_lt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j] + 1./tC.eta_H[i+1][j-1] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) );
            eta_H_corner_rb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j+1] + 1./tC.eta_H[i-1][j] ) );
            eta_H_corner_rt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j+1] + 1./tC.eta_H[i+1][j] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) );
            Jr_av_phi_lb = 0.5*( dm.r[i]*J[0][i][j] + dm.r[i-1]*J[0][i-1][j] )/r_imh;
            Jr_av_phi_lt = 0.5*( dm.r[i+1]*J[0][i+1][j] + dm.r[i]*J[0][i][j] )/r_iph;
            Jr_av_phi_rb = 0.5*( dm.r[i]*J[0][i][j+1] + dm.r[i-1]*J[0][i-1][j+1] )/r_imh;
            Jr_av_phi_rt = 0.5*( dm.r[i+1]*J[0][i+1][j+1] + dm.r[i]*J[0][i][j+1] )/r_iph;
            vr_av_phi_lb = -eta_H_corner_lb*Jr_av_phi_lb;
            vr_av_phi_lt = -eta_H_corner_lt*Jr_av_phi_lt;
            vr_av_phi_rb = -eta_H_corner_rb*Jr_av_phi_rb;
            vr_av_phi_rt = -eta_H_corner_rt*Jr_av_phi_rt;
            vtheta_av_phi_lb = -eta_H_corner_lb*Jtheta_av_phi_lb;
            vtheta_av_phi_lt = -eta_H_corner_lt*Jtheta_av_phi_lt;
            vtheta_av_phi_rb = -eta_H_corner_rb*Jtheta_av_phi_rb;
            vtheta_av_phi_rt = -eta_H_corner_rt*Jtheta_av_phi_rt;
            vphi_phi_lb = -eta_H_corner_lb*J[2][i][j];
            vphi_phi_lt = -eta_H_corner_lt*J[2][i+1][j];
            vphi_phi_rb = -eta_H_corner_rb*J[2][i][j+1];
            vphi_phi_rt = -eta_H_corner_rt*J[2][i+1][j+1];

            Btheta_up_phi_lb = Reconstruct_Btheta(B, i, j, vr_av_phi_lb, Deltar_1, Deltar_2, Deltar_3, dm);
            Btheta_up_phi_lt = Reconstruct_Btheta(B, i+1, j, vr_av_phi_lt, Deltar_2, Deltar_3, Deltar_4, dm);
            Btheta_up_phi_rb = Reconstruct_Btheta(B, i, j+1, vr_av_phi_rb, Deltar_1, Deltar_2, Deltar_3, dm);
            Btheta_up_phi_rt = Reconstruct_Btheta(B, i+1, j+1, vr_av_phi_rt, Deltar_2, Deltar_3, Deltar_4, dm);

            Br_up_phi_lb = Reconstruct_BrEth(B, i, j, vtheta_av_phi_lb, Deltatheta_1, Deltatheta_2, Deltatheta_3, dm);
            Br_up_phi_lt = Reconstruct_BrEth(B, i+1, j, vtheta_av_phi_lt, Deltatheta_1, Deltatheta_2, Deltatheta_3, dm);
            Br_up_phi_rb = Reconstruct_BrEth(B, i, j+1, vtheta_av_phi_rb, Deltatheta_2, Deltatheta_3, Deltatheta_4, dm);
            Br_up_phi_rt = Reconstruct_BrEth(B, i+1, j+1, vtheta_av_phi_rt, Deltatheta_2, Deltatheta_3, Deltatheta_4, dm);

            double E_theta_b, E_theta_f, E_r_b, E_r_f;
            E_theta_b = 0.5*( -vphi_phi_lb*Br_up_phi_lb + -vphi_phi_rb*Br_up_phi_rb );
            E_theta_f = 0.5*( -vphi_phi_lt*Br_up_phi_lt + -vphi_phi_rt*Br_up_phi_rt );
            E_r_b = -0.5*( -vphi_phi_lb*Btheta_up_phi_lb + -vphi_phi_lt*Btheta_up_phi_lt );
            E_r_f = -0.5*( -vphi_phi_rb*Btheta_up_phi_rb + -vphi_phi_rt*Btheta_up_phi_rt );

            Fr += ( r_iph*E_theta_f - r_imh*E_theta_b );
            Ftheta += -( E_r_f - E_r_b );

            //Quantities at location of J_phi
            // eta_H_corner_lb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j] + 1./tC.eta_H[i-1][j-1] ) );
            // eta_H_corner_lt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j] + 1./tC.eta_H[i+1][j-1] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) );
            // eta_H_corner_rb = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j+1] + 1./tC.eta_H[i-1][j] ) );
            // eta_H_corner_rt = 2.*( dm.Deltar[i+1]+dm.Deltar[i] )/( dm.Deltar[i+1]*( 1./tC.eta_H[i+1][j+1] + 1./tC.eta_H[i+1][j] ) + dm.Deltar[i]*( 1./tC.eta_H[i][j+1] + 1./tC.eta_H[i][j] ) );
            // vphi_phi_lb = -eta_H_corner_lb*J[2][i][j];
            // vphi_phi_lt = -eta_H_corner_lt*J[2][i+1][j];
            // vphi_phi_rb = -eta_H_corner_rb*J[2][i][j+1];
            // vphi_phi_rt = -eta_H_corner_rt*J[2][i+1][j+1];
            // vphi_r_b = 0.5*( vphi_phi_lb + vphi_phi_rb );
            // vphi_r_f = 0.5*( vphi_phi_lt + vphi_phi_rt );
            // vphi_theta_b = 0.5*( vphi_phi_lb + vphi_phi_lt );
            // vphi_theta_f = 0.5*( vphi_phi_rb + vphi_phi_rt );
            // Fr += ( r_iph*(-vphi_r_f*B[0][i+1][j]) - r_imh*(-vphi_r_b*B[0][i][j]) );
            // Ftheta += ( (-vphi_theta_f*B[1][i][j+1]) - (-vphi_theta_b*B[1][i][j]) );

            mC.E_H_Pol[0][i][j] = E_r_HB_jmh + E_r_HB_mag_jph;//*EdgeLimit - (-vphi_theta_b*B[1][i][j])*( 1. - EdgeLimit );
            mC.E_H_Pol[1][i][j] = E_theta_HB_imh + E_theta_HB_mag_imh;//*EdgeLimit + (-vphi_r_b*B[0][i][j])*( 1. - EdgeLimit );
            mC.E_H_Pol[0][i][j+1] = E_r_HB_jph + E_r_HB_mag_jph;//*EdgeLimit - (-vphi_theta_f*B[1][i][j+1])*( 1. - EdgeLimit );
            mC.E_H_Pol[1][i+1][j] = E_theta_HB_iph + E_theta_HB_mag_iph;//*EdgeLimit + (-vphi_r_f*B[0][i+1][j])*( 1. - EdgeLimit );

            lapBphi = 1./dm.n_r[i]/(dm.r[i]*dm.r[i])*( 2.*pow(r_iph,2.)/(dm.n_r[i+1]+dm.n_r[i])*(B[2][i+1][j]*dm.n_t[i+1]-B[2][i][j]*dm.n_t[i])/(0.5*dm.Deltar[i+1]+0.5*dm.Deltar[i])
                                                        - 2.*pow(r_imh,2.)/(dm.n_r[i]+dm.n_r[i-1])*(B[2][i][j]*dm.n_t[i]-B[2][i-1][j]*dm.n_t[i-1])/(0.5*dm.Deltar[i]+0.5*dm.Deltar[i-1]) )/dm.Deltar[i]
                                                        + dm.n_t[i]/(dm.r[i]*dm.r[i]*Deltatheta)*( (B[2][i][j+1]-B[2][i][j])/Deltatheta - (B[2][i][j]-B[2][i][j-1])/Deltatheta + (B[2][i][j+1]-B[2][i][j-1])/(2.*tan(dm.theta[j])) );
            if (i == N_GC) {
                lapBphi_im1 = 1./dm.n_r[i-1]/pow(dm.r[i]-dm.Deltar[i],2.)*( 2.*pow(dm.r[i]-0.5*dm.Deltar[i-1],2.)/(dm.n_r[i]+dm.n_r[i-1])*(B[2][i][j]*dm.n_t[i]-B[2][i-1][j]*dm.n_t[i-1])/(0.5*dm.Deltar[i]+0.5*dm.Deltar[i-1])
                                                                - 2.*pow(dm.r[i]-1.5*dm.Deltar[i-1],2.)/(dm.n_r[i-1]+dm.n_r[i-2])*(B[2][i-1][j]*dm.n_t[i-1]-B[2][i-2][j]*dm.n_t[i-2])/(0.5*dm.Deltar[i-1]+0.5*dm.Deltar[i-2]) )/dm.Deltar[i-1]
                                                                + dm.n_t[i-1]/(pow(dm.r[i]-dm.Deltar[i],2.)*Deltatheta)*( (B[2][i-1][j+1]-B[2][i-1][j])/Deltatheta - (B[2][i-1][j]-B[2][i-1][j-1])/Deltatheta + (B[2][i-1][j+1]-B[2][i-1][j-1])/(2.*tan(dm.theta[j])) );
            }
            else {
                lapBphi_im1 = 1./dm.n_r[i-1]/(dm.r[i-1]*dm.r[i-1])*( 2.*pow(dm.r[i-1]+0.5*dm.Deltar[i-1],2.)/(dm.n_r[i]+dm.n_r[i-1])*(B[2][i][j]*dm.n_t[i]-B[2][i-1][j]*dm.n_t[i-1])/(0.5*dm.Deltar[i]+0.5*dm.Deltar[i-1])
                                                                - 2.*pow(dm.r[i-1]-0.5*dm.Deltar[i-1],2.)/(dm.n_r[i-1]+dm.n_r[i-2])*(B[2][i-1][j]*dm.n_t[i-1]-B[2][i-2][j]*dm.n_t[i-2])/(0.5*dm.Deltar[i-1]+0.5*dm.Deltar[i-2]) )/dm.Deltar[i-1]
                                                                + dm.n_t[i-1]/(dm.r[i-1]*dm.r[i-1]*Deltatheta)*( (B[2][i-1][j+1]-B[2][i-1][j])/Deltatheta - (B[2][i-1][j]-B[2][i-1][j-1])/Deltatheta + (B[2][i-1][j+1]-B[2][i-1][j-1])/(2.*tan(dm.theta[j])) );
            }
            if (i == B.shape()[1]-N_GC-2) {
                lapBphi_ip1 = 1./dm.n_r[i+1]/pow(dm.r[i]+dm.Deltar[i],2.)*( 2.*pow(dm.r[i]+1.5*dm.Deltar[i+1],2.)/(dm.n_r[i+2]+dm.n_r[i+1])*(B[2][i+2][j]*dm.n_t[i+2]-B[2][i+1][j]*dm.n_t[i+1])/(0.5*dm.Deltar[i+2]+0.5*dm.Deltar[i+1])
                                                               - 2.*pow(dm.r[i]+0.5*dm.Deltar[i+1],2.)/(dm.n_r[i+1]+dm.n_r[i])*(B[2][i+1][j]*dm.n_t[i+1]-B[2][i][j]*dm.n_t[i])/(0.5*dm.Deltar[i+1]+0.5*dm.Deltar[i]) )/dm.Deltar[i+1]
                                                                + dm.n_t[i+1]/(pow(dm.r[i]+dm.Deltar[i],2.)*Deltatheta)*( (B[2][i+1][j+1]-B[2][i+1][j])/Deltatheta - (B[2][i+1][j]-B[2][i+1][j-1])/Deltatheta + (B[2][i+1][j+1]-B[2][i+1][j-1])/(2.*tan(dm.theta[j])) );
            }
            else {
                lapBphi_ip1 = 1./dm.n_r[i+1]/(dm.r[i+1]*dm.r[i+1])*( 2.*pow(dm.r[i+1]+0.5*dm.Deltar[i+1],2.)/(dm.n_r[i+2]+dm.n_r[i+1])*(B[2][i+2][j]*dm.n_t[i+2]-B[2][i+1][j]*dm.n_t[i+1])/(0.5*dm.Deltar[i+2]+0.5*dm.Deltar[i+1])
                                                               - 2.*pow(dm.r[i+1]-0.5*dm.Deltar[i+1],2.)/(dm.n_r[i+1]+dm.n_r[i])*(B[2][i+1][j]*dm.n_t[i+1]-B[2][i][j]*dm.n_t[i])/(0.5*dm.Deltar[i+1]+0.5*dm.Deltar[i]) )/dm.Deltar[i+1]
                                                                + dm.n_t[i+1]/(dm.r[i+1]*dm.r[i+1]*Deltatheta)*( (B[2][i+1][j+1]-B[2][i+1][j])/Deltatheta - (B[2][i+1][j]-B[2][i+1][j-1])/Deltatheta + (B[2][i+1][j+1]-B[2][i+1][j-1])/(2.*tan(dm.theta[j])) );
            }
            lapBphi_jm1 = 1./dm.n_r[i]/(dm.r[i]*dm.r[i])*( 2.*pow(r_iph,2.)/(dm.n_r[i+1]+dm.n_r[i])*(B[2][i+1][j]*dm.n_t[i+1]-B[2][i][j]*dm.n_t[i])/(0.5*dm.Deltar[i+1]+0.5*dm.Deltar[i])
                                                         - 2.*pow(r_imh,2.)/(dm.n_r[i]+dm.n_r[i-1])*(B[2][i][j]*dm.n_t[i]-B[2][i-1][j]*dm.n_t[i-1])/(0.5*dm.Deltar[i]+0.5*dm.Deltar[i-1]) )/dm.Deltar[i]
                            + dm.n_t[i]/(dm.r[i]*dm.r[i]*Deltatheta)*( (B[2][i][j]-B[2][i][j-1])/Deltatheta - (B[2][i][j-1]-B[2][i][j-2])/Deltatheta + (B[2][i][j]-B[2][i][j-2])/(2.*tan(dm.theta[j-1])) );
            lapBphi_jp1 = 1./dm.n_r[i]/(dm.r[i]*dm.r[i])*( 2.*pow(r_iph,2.)/(dm.n_r[i+1]+dm.n_r[i])*(B[2][i+1][j]*dm.n_t[i+1]-B[2][i][j]*dm.n_t[i])/(0.5*dm.Deltar[i+1]+0.5*dm.Deltar[i])
                                                        - 2.*pow(r_imh,2.)/(dm.n_r[i]+dm.n_r[i-1])*(B[2][i][j]*dm.n_t[i]-B[2][i-1][j]*dm.n_t[i-1])/(0.5*dm.Deltar[i]+0.5*dm.Deltar[i-1]) )/dm.Deltar[i]
                            + dm.n_t[i]/(dm.r[i]*dm.r[i]*Deltatheta)*( (B[2][i][j+2]-B[2][i][j+1])/Deltatheta - (B[2][i][j+1]-B[2][i][j])/Deltatheta + (B[2][i][j+2]-B[2][i][j])/(2.*tan(dm.theta[j+1])) );

            //Hyperdiffusion term
            // eta_hyp = tC.C_hyp*Delta2*tC.eta_H[i][j]*OuterEdgeLimiter(dm.r[i],dm.r_max,dm.r_max/1000.);
            eta_hyp = tC.C_hyp*Delta2*tC.eta_H[i][j]*EdgeLimiter(dm.r[i],dm.r_max,dm.r_min,dm.r_max/100.)*PoleEdgeLimiter(dm.theta[j],pi/180.);
            Fr += eta_hyp*1./dm.n_r[i]*( ( (dm.r[i]+dm.Deltar[i])*lapBphi_ip1 - dm.r[i]*lapBphi)/(0.5*(dm.Deltar[i+1]+dm.Deltar[i]))
                                                                        - (dm.r[i]*lapBphi - (dm.r[i]-dm.Deltar[i])*lapBphi_im1)/(0.5*(dm.Deltar[i]+dm.Deltar[i-1])) );
            Ftheta += eta_hyp*( (sin(dm.theta[j]+Deltatheta)*lapBphi_jp1 - sin(dm.theta[j])*lapBphi)/(dm.r[i]*sin(dm.theta[j])*Deltatheta)
                                                    - (sin(dm.theta[j])*lapBphi - sin(dm.theta[j]-Deltatheta)*lapBphi_jm1)/(dm.r[i]*sin(dm.theta[j])*Deltatheta) );

            if(tC.conductivity_anisotropy == true){
                if(mC.quantization == true){
                    eta_O_delta_th_ip1 = (dm.Deltar[i+1]+dm.Deltar[i])/( dm.Deltar[i+1]/tC.eta_O_delta[i+1][j] + dm.Deltar[i]/tC.eta_O_delta[i][j] );
                    eta_O_delta_th_im1 = (dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]/tC.eta_O_delta[i][j] + dm.Deltar[i-1]/tC.eta_O_delta[i-1][j] );
                    eta_O_delta_r_jp1 = 2./( 1./tC.eta_O_delta[i][j+1] + 1./tC.eta_O_delta[i][j] );
                    eta_O_delta_r_jm1 = 2./( 1./tC.eta_O_delta[i][j] + 1./tC.eta_O_delta[i][j-1] );
                    eta_O_perp_th_ip1 = (dm.Deltar[i+1]+dm.Deltar[i])/( dm.Deltar[i+1]/tC.eta_O_perp[i+1][j] + dm.Deltar[i]/tC.eta_O_perp[i][j] );
                    eta_O_perp_th_im1 = (dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]/tC.eta_O_perp[i][j] + dm.Deltar[i-1]/tC.eta_O_perp[i-1][j] );
                    eta_O_perp_r_jp1 = 2./( 1./tC.eta_O_perp[i][j+1] + 1./tC.eta_O_perp[i][j] );
                    eta_O_perp_r_jm1 = 2./( 1./tC.eta_O_perp[i][j] + 1./tC.eta_O_perp[i][j-1] );

                    Qphi[i][j] = - Fr/(dm.n_r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*Deltatheta)
                                    - ( eta_O_delta_th_ip1*thdot(J,mC.thBhat,i+1,j)*mC.thBhat[1][i+1][j]*r_iph - eta_O_delta_th_im1*thdot(J,mC.thBhat,i,j)*mC.thBhat[1][i][j]*r_imh )/(dm.n_r[i]*dm.r[i]*dm.Deltar[i])
                                    - ( eta_O_perp_th_ip1*J[1][i+1][j]*r_iph - eta_O_perp_th_im1*J[1][i][j]*r_imh )/(dm.n_r[i]*dm.r[i]*dm.Deltar[i])
                                    + ( eta_O_delta_r_jp1*rdot(J,mC.rBhat,i,j+1)*mC.rBhat[0][i][j+1] - eta_O_delta_r_jm1*rdot(J,mC.rBhat,i,j)*mC.rBhat[0][i][j] )/(dm.r[i]*Deltatheta)
                                    + ( eta_O_perp_r_jp1*J[0][i][j+1] - eta_O_perp_r_jm1*J[0][i][j] )/(dm.r[i]*Deltatheta);
                }
                else{
                    eta_O_perp_th_ip1 = (dm.Deltar[i+1]+dm.Deltar[i])/( dm.Deltar[i+1]/tC.eta_O_perp[i+1][j] + dm.Deltar[i]/tC.eta_O_perp[i][j] );
                    eta_O_perp_th_im1 = (dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]/tC.eta_O_perp[i][j] + dm.Deltar[i-1]/tC.eta_O_perp[i-1][j] );
                    eta_O_perp_r_jp1 = 2./( 1./tC.eta_O_perp[i][j+1] + 1./tC.eta_O_perp[i][j] );
                    eta_O_perp_r_jm1 = 2./( 1./tC.eta_O_perp[i][j] + 1./tC.eta_O_perp[i][j-1] );

                    Qphi[i][j] = - Fr/(dm.n_r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*Deltatheta)
                                    - ( eta_O_perp_th_ip1*J[1][i+1][j]*r_iph - eta_O_perp_th_im1*J[1][i][j]*r_imh )/(dm.n_r[i]*dm.r[i]*dm.Deltar[i])
                                    + ( eta_O_perp_r_jp1*J[0][i][j+1] - eta_O_perp_r_jm1*J[0][i][j] )/(dm.r[i]*Deltatheta);
                }
            }
            else{
                eta_O_th_ip1 = (dm.Deltar[i+1]+dm.Deltar[i])/( dm.Deltar[i+1]/tC.eta_O[i+1][j] + dm.Deltar[i]/tC.eta_O[i][j] );
                eta_O_th_im1 = (dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]/tC.eta_O[i][j] + dm.Deltar[i-1]/tC.eta_O[i-1][j] );
                eta_O_r_jp1 = 2./( 1./tC.eta_O[i][j+1] + 1./tC.eta_O[i][j] );
                eta_O_r_jm1 = 2./( 1./tC.eta_O[i][j] + 1./tC.eta_O[i][j-1] );

                Qphi[i][j] = - Fr/(dm.n_r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*Deltatheta)
                                - ( eta_O_th_ip1*J[1][i+1][j]*r_iph - eta_O_th_im1*J[1][i][j]*r_imh )/(dm.n_r[i]*dm.r[i]*dm.Deltar[i])
                                + ( eta_O_r_jp1*J[0][i][j+1] - eta_O_r_jm1*J[0][i][j] )/(dm.r[i]*Deltatheta);

            }
            q_SH[i][j] = 0*( - 0.5*(beta_r + beta_r_m1)/(96.*pi)*0.5*(dm.n_t[i]+dm.n_t[i-1])*pow(B_im - B_im1,3.)*4.*pi*r_imh*r_imh*sin(dm.theta[j])*sin(0.5*Deltatheta)*step( -(beta_r + beta_r_m1)*(B_im-B_im1) )
                         - 0.5*(beta_r_p1 + beta_r)/(96.*pi)*0.5*(dm.n_t[i+1]+dm.n_t[i])*pow(B_ip1 - B_ip,3.)*4.*pi*r_iph*r_iph*sin(dm.theta[j])*sin(0.5*Deltatheta)*step( -(beta_r_p1 + beta_r)*(B_ip1-B_ip) )
                         - 0.5*(beta_theta + beta_theta_m1)/(96.*pi)*dm.n_t[i]*pow(B_jm - B_jm1,3.)*2.*pi*sin(dm.theta[j]-0.5*Deltatheta)*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*step( -(beta_theta + beta_theta_m1)*(B_jm-B_jm1) )
                         - 0.5*(beta_theta_p1 + beta_theta)/(96.*pi)*dm.n_t[i]*pow(B_jp1 - B_jp,3.)*2.*pi*sin(dm.theta[j]+0.5*Deltatheta)*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*step( -(beta_theta_p1 + beta_theta)*(B_jp1-B_jp) ) );

        }
    }

    //Exchange Qphi values in inter-process ghost cells and fill theta boundary ghost cells
    exchng2Scalar(Qphi, N_GC, process, dm.stridetype_Sca);

    return;
}

/*
        Godunov scheme with slope-limited reconstruction for computing Hall contribution to E_r and E_theta for B_tor evolution

        Inputs: u_L, u_R: left and right
                beta_L, B_L: wave propagation speed and reconstructed field at left side of cell interface
                beta_R, B_R: wave propagation speed and reconstructed field at right side of cell interface
        Output: Electric field (either E_r or E_theta) at an interface as determined by the Godunov method
 */
double Godunov(double u_L, double u_R, double beta_L, double B_L, double beta_R, double B_R){

    double F_G, s;
    double tol = 0.;

    if( u_L < u_R ){
        if( u_L > tol ){
            F_G = 0.5*u_L*B_L;
        }
        else if( u_R < -tol ){
            F_G = 0.5*u_R*B_R;
        }
        else{
            F_G = 0.;
        }
    }
    else{
        s = 0.5*(beta_L+beta_R)*(B_L + B_R)/2.;
        if( s > tol ){
            F_G = 0.5*u_L*B_L;
        }
        else if ( s < -tol ){
            F_G = 0.5*u_R*B_R;
        }
        else{
            F_G = 0.;
        }
    }

    return F_G;
}

/*
    Reconstruct B_theta from cell face centers to cell edge centers (the location of J_phi). Uses cell centroids in spherical coordinates from Mignone J. Comp. Phys. 270, 784 (2014)

    Inputs: B: magnetic field in reduced units at cell face centers
            i, j: cell indices labelling which B-field values are used for reconstruction. Reconstruction occurs at cell edge labelled by i-1/2, j-1/2
            vr_av_phi: average Hall advection velocity at cell edge to which B_theta is reconstructed
            Deltar_1, Deltar_2, Deltar_3: twice the cell-centroid spacings in radial direction. Include corrections for variable cell spacing.
            dm: Domain object containing information about the simulation domain
    Output: Btheta_up_phi: reconstructed magnetic field at location of J_phi
 */
double Reconstruct_Btheta(VectorField & B, size_t i, size_t j, double vr_av_phi, double Deltar_1, double Deltar_2, double Deltar_3, const Domain & dm)
{
    double Btheta_up_phi, alphatheta_im1, alphatheta_i;
    double Delta_xi_p1, Delta_xi_m1;
    double tol = 1e-9;

    if(vr_av_phi > tol){
        Delta_xi_p1 = 0.5*dm.Deltar[i-1]*( 1. - 4.*dm.r[i-1]*dm.Deltar[i-1]/(12.*dm.r[i-1]*dm.r[i-1]+dm.Deltar[i-1]*dm.Deltar[i-1]) );
        Delta_xi_m1 = 0.5*dm.Deltar[i-1]*( 1. + 4.*dm.r[i-1]*dm.Deltar[i-1]/(12.*dm.r[i-1]*dm.r[i-1]+dm.Deltar[i-1]*dm.Deltar[i-1]) );
        alphatheta_im1 = MC( ( B[1][i][j] - B[1][i-1][j] )/Deltar_2 + ( B[1][i-1][j] - B[1][i-2][j] )/Deltar_1,
                             ( B[1][i][j] - B[1][i-1][j] )/Delta_xi_p1, ( B[1][i-1][j] - B[1][i-2][j] )/Delta_xi_m1 );
        Btheta_up_phi = B[1][i-1][j] + Delta_xi_p1*alphatheta_im1;
    }
    else if(vr_av_phi < -tol){
        Delta_xi_p1 = 0.5*dm.Deltar[i]*( 1. - 4.*dm.r[i]*dm.Deltar[i]/(12.*dm.r[i]*dm.r[i]+dm.Deltar[i]*dm.Deltar[i]) );
        Delta_xi_m1 = 0.5*dm.Deltar[i]*( 1. + 4.*dm.r[i]*dm.Deltar[i]/(12.*dm.r[i]*dm.r[i]+dm.Deltar[i]*dm.Deltar[i]) );
        alphatheta_i = MC( ( B[1][i+1][j] - B[1][i][j] )/Deltar_3 + ( B[1][i][j] - B[1][i-1][j] )/Deltar_2,
                           ( B[1][i+1][j] - B[1][i][j] )/Delta_xi_p1, ( B[1][i][j] - B[1][i-1][j] )/Delta_xi_m1 );
        Btheta_up_phi = B[1][i][j] - Delta_xi_m1*alphatheta_i;
    }
    else{
        Btheta_up_phi = (B[1][i-1][j] + B[1][i][j])/2.;
    }

    return Btheta_up_phi;
}

/*
    Reconstruct B_r from cell face centers to cell edge centers (the location of J_phi). Uses cell centroids in spherical coordinates from Mignone J. Comp. Phys. 270, 784 (2014)

    Inputs: B: magnetic field in reduced units at cell face centers
            i, j: cell indices labelling which B-field values are used for reconstruction. Reconstruction occurs at cell edge labelled by i-1/2, j-1/2
            vtheta_av_phi: average Hall advection velocity at cell edge to which B_r is reconstructed
            Deltatheta_1, Deltatheta_2, Deltatheta_3: twice the cell-centroid spacings in theta direction.
    Output: Br_up_phi: reconstructed magnetic field at location of J_phi
 */
double Reconstruct_Br(VectorField & B, size_t i, size_t j, double vtheta_av_phi, double Deltatheta_1, double Deltatheta_2, double Deltatheta_3, const Domain & dm) {

    double Br_up_phi, alphar_jm1, alphar_j;
    double Br_up_phi_L, Br_up_phi_R;
    double Delta_xi_p1, Delta_xi_m1;
    double tol = 1e-9;

    if(vtheta_av_phi > tol){
        Delta_xi_m1 = 0.5*dm.Deltatheta + 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = 0.5*dm.Deltatheta - 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        alphar_jm1 = MC( (B[0][i][j]-B[0][i][j-1])/Deltatheta_2 + (B[0][i][j-1]-B[0][i][j-2])/Deltatheta_1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_p1, (B[0][i][j-1]-B[0][i][j-2])/Delta_xi_m1 );
        Br_up_phi = B[0][i][j-1] + Delta_xi_p1*alphar_jm1;
    }
    else if(vtheta_av_phi < -tol){
        Delta_xi_m1 = 0.5*dm.Deltatheta + 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = 0.5*dm.Deltatheta - 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        alphar_j = MC( (B[0][i][j+1]-B[0][i][j])/Deltatheta_3 + (B[0][i][j]-B[0][i][j-1])/Deltatheta_2, (B[0][i][j+1]-B[0][i][j])/Delta_xi_p1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_m1 );
        Br_up_phi = B[0][i][j] - Delta_xi_m1*alphar_j;
    }
    else{
        Delta_xi_m1 = 0.5*dm.Deltatheta + 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = 0.5*dm.Deltatheta - 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        alphar_jm1 = MC( (B[0][i][j]-B[0][i][j-1])/Deltatheta_2 + (B[0][i][j-1]-B[0][i][j-2])/Deltatheta_1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_p1, (B[0][i][j-1]-B[0][i][j-2])/Delta_xi_m1 );
        Br_up_phi_L = B[0][i][j-1] + Delta_xi_p1*alphar_jm1;

        Delta_xi_m1 = 0.5*dm.Deltatheta + 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = 0.5*dm.Deltatheta - 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        alphar_j = MC( (B[0][i][j+1]-B[0][i][j])/Deltatheta_3 + (B[0][i][j]-B[0][i][j-1])/Deltatheta_2, (B[0][i][j+1]-B[0][i][j])/Delta_xi_p1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_m1 );
        Br_up_phi_R = B[0][i][j] - Delta_xi_m1*alphar_j;

        Br_up_phi = (Br_up_phi_L + Br_up_phi_R)/2.;
    }

    return Br_up_phi;
}

/*
    Reconstruct B_r from cell face centers to cell edge centers (the location of J_phi). Uses cell centroids in spherical coordinates from Mignone J. Comp. Phys. 270, 784 (2014).
    Used in computing E_theta inside B_torEvolve only. Differs from Reconstruct_BrEth only in tol parameter

    Inputs: B: magnetic field in reduced units at cell face centers
            i, j: cell indices labelling which B-field values are used for reconstruction. Reconstruction occurs at cell edge labelled by i-1/2, j-1/2
            vtheta_av_phi: average Hall advection velocity at cell edge to which B_r is reconstructed
            Deltatheta_1, Deltatheta_2, Deltatheta_3: twice the cell-centroid spacings in theta direction.
    Output: Br_up_phi: reconstructed magnetic field at location of J_phi
 */
double Reconstruct_BrEth(VectorField & B, size_t i, size_t j, double vtheta_av_phi, double Deltatheta_1, double Deltatheta_2, double Deltatheta_3, const Domain & dm) {

    double Br_up_phi, alphar_jm1, alphar_j;
    double Br_up_phi_L, Br_up_phi_R;
    double Delta_xi_p1, Delta_xi_m1;
    static double tol = 1e-5;
    static double halfDeltatheta = 0.5*dm.Deltatheta;

    if(vtheta_av_phi > tol){
        Delta_xi_m1 = halfDeltatheta + 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = halfDeltatheta - 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        alphar_jm1 = MC( (B[0][i][j]-B[0][i][j-1])/Deltatheta_2 + (B[0][i][j-1]-B[0][i][j-2])/Deltatheta_1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_p1, (B[0][i][j-1]-B[0][i][j-2])/Delta_xi_m1 );
        Br_up_phi = B[0][i][j-1] + Delta_xi_p1*alphar_jm1;
    }
    else if(vtheta_av_phi < -tol){
        Delta_xi_m1 = halfDeltatheta + 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = halfDeltatheta - 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        alphar_j = MC( (B[0][i][j+1]-B[0][i][j])/Deltatheta_3 + (B[0][i][j]-B[0][i][j-1])/Deltatheta_2, (B[0][i][j+1]-B[0][i][j])/Delta_xi_p1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_m1 );
        Br_up_phi = B[0][i][j] - Delta_xi_m1*alphar_j;
    }
    else{
        Delta_xi_m1 = halfDeltatheta + 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = halfDeltatheta - 1./tan(dm.theta[j-1])*dm.Deltatheta_dfactor;
        alphar_jm1 = MC( (B[0][i][j]-B[0][i][j-1])/Deltatheta_2 + (B[0][i][j-1]-B[0][i][j-2])/Deltatheta_1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_p1, (B[0][i][j-1]-B[0][i][j-2])/Delta_xi_m1 );
        Br_up_phi_L = B[0][i][j-1] + Delta_xi_p1*alphar_jm1;

        Delta_xi_m1 = halfDeltatheta + 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        Delta_xi_p1 = halfDeltatheta - 1./tan(dm.theta[j])*dm.Deltatheta_dfactor;
        alphar_j = MC( (B[0][i][j+1]-B[0][i][j])/Deltatheta_3 + (B[0][i][j]-B[0][i][j-1])/Deltatheta_2, (B[0][i][j+1]-B[0][i][j])/Delta_xi_p1, (B[0][i][j]-B[0][i][j-1])/Delta_xi_m1 );
        Br_up_phi_R = B[0][i][j] - Delta_xi_m1*alphar_j;

        Br_up_phi = (Br_up_phi_L + Br_up_phi_R)/2.;
    }

    return Br_up_phi;

}

/*
    Computes c*lapse*electric field components along cell edges in reduced units i.e., c*t_0/(L_0*B_0)*E*exp(nu/2)
    Works for both isotropic and anisotropic conductivity
    Uses cell centroids in spherical coordinates from Mignone J. Comp. Phys. 270, 784 (2014)

    Inputs: B, J: magnetic field and redshifted current density in reduced units at cell face centers and cell face edges respectively
            Bn: updated intermediate values of B, in which only toroidal field is included
            N_GC: number of ghost cells
            tC: TransCoeffs object containing Ohmic and Hall diffusivities in reduced units
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            dm: Domain object containing information about the simulation domain
            process: Process objects containing information about the simulation domain and the current process
    Output: phE: electric field on cell face edge at location of J_phi (where E_phi lives, E_r and E_theta interpolated here)
            cE: conjugate electric field on cell face edge at location of J_phi (where E_phi lives, E_r and E_theta interpolated here)

*/
void Compute_E(VectorField & B, VectorField & Bn, VectorField & phE, VectorField & cE, VectorField & J, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process) {
    double Deltatheta = dm.Deltatheta;
    double vr_av_phi, vtheta_av_phi, vphi_av_phi; //averaged values of the velocity = -current density where J_phi lives.
    double Bphi_phi; //interpolated field values (B_phi interpolated to where J_phi lives).
    double Br_up_phi, Btheta_up_phi; //upwinded field values (B_r and B_theta upwinded to where J_phi lives).

    double eta_O_corner, eta_O_corner_perp, eta_O_corner_delta; //corner-averaged Ohmic diffusivities (where J_phi lives). Use the first variable for isotropic conductivity and the latter two for anisotropic conductivity.
    double eta_H_corner; //edge-averaged Hall diffusivities at corner where J_phi lives

    double Deltar_1, Deltar_2, Deltar_3; //twice the r-spacing between cell centroids of i-2 and i-1, i-1 and i, i and i+1 cells, respectively.
    double Deltatheta_1, Deltatheta_2, Deltatheta_3; //twice the theta-spacing between cell centroids of j-2 and j-1, j-1 and j, j and j+1 cells, respectively.
    double r_m, theta_m;
    double Jr_av_phi, Jtheta_av_phi; //interpolated currents: in theta, phi-directions where J_r lives; in r, phi-directions where J_th lives; in r, theta-directions where J_phi lives
    double Jr_av_c, Jtheta_av_c, vr_av_c, vtheta_av_c, vphi_av_c;
    double JdotBhat_phi, JdotBhat_c; //dot product of current density and Bhat evaluated where J_phi lives, and conjugate version of this
    //Compute electric field and conjugate electric field at cell edges using J and B
    //phE[0][i][j] = E_r^{i-1/2,j-1/2}, phE[1][i][j] = E_theta^{i-1/2,j-1/2}, phE[2][i][j] = E_phi^{i-1/2,j-1/2} (native)
    //cE[0][i][j] = E_r^{i-1/2,j-1/2}, phE[1][i][j] = E_theta^{i-1/2,j-1/2}, phE[2][i][j] = E_phi^{i-1/2,j-1/2} (native)

    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        r_m = dm.r[i] - 0.5*dm.Deltar[i];

        Deltar_1 = dm.Deltar[i-1] + dm.Deltar[i-2] + 1./3.*( dm.r[i-1]*dm.Deltar[i-1]*dm.Deltar[i-1]/(dm.r[i-1]*dm.r[i-1] + dm.Deltar[i-1]*dm.Deltar[i-1]/12.)
                                                             - dm.r[i-2]*dm.Deltar[i-2]*dm.Deltar[i-2]/(dm.r[i-2]*dm.r[i-2] + dm.Deltar[i-2]*dm.Deltar[i-2]/12.) );
        Deltar_2 = dm.Deltar[i] + dm.Deltar[i-1] + 1./3.*( dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i] + dm.Deltar[i]*dm.Deltar[i]/12.)
                                                           - dm.r[i-1]*dm.Deltar[i-1]*dm.Deltar[i-1]/(dm.r[i-1]*dm.r[i-1] + dm.Deltar[i-1]*dm.Deltar[i-1]/12.) );
        Deltar_3 = dm.Deltar[i+1] + dm.Deltar[i] + 1./3.*( dm.r[i+1]*dm.Deltar[i+1]*dm.Deltar[i+1]/(dm.r[i+1]*dm.r[i+1] + dm.Deltar[i+1]*dm.Deltar[i+1]/12.)
                                                           - dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i] + dm.Deltar[i]*dm.Deltar[i]/12.) );

        //Also computes E values at right boundary j = B.shape()[2]-N_GC;. This will only be retained for theta=pi, otherwise overwritten by ghost cell exchanges.
        for(size_t j=N_GC; j<B.shape()[2]-N_GC+1; j++) {
            theta_m = dm.theta[j] - 0.5*dm.Deltatheta;

            Deltatheta_1 = 2.*( Deltatheta + ( 1./tan(dm.theta[j-1]) - 1./tan(dm.theta[j-2]) )*dm.Deltatheta_dfactor );
            Deltatheta_2 = 2.*( Deltatheta + ( 1./tan(dm.theta[j]) - 1./tan(dm.theta[j-1]) )*dm.Deltatheta_dfactor );
            Deltatheta_3 = 2.*( Deltatheta + ( 1./tan(dm.theta[j+1]) - 1./tan(dm.theta[j]) )*dm.Deltatheta_dfactor );

            //Jtheta interpolation onto cell corners that lie on theta=0,pi boundaries
            if( process.nbrleft < 0 && j == N_GC ){
                Jtheta_av_phi = 0.;
                Jtheta_av_c = 0.;
                Bphi_phi = 0.;
            }
            else if( process.nbrright < 0 && j == B.shape()[2]-N_GC ){
                Jtheta_av_phi = 0.;
                Jtheta_av_c = 0.;
                Bphi_phi = 0.;
            }
            else{
                Jtheta_av_phi = 0.5*( sin(dm.theta[j])*J[1][i][j] + sin(dm.theta[j-1])*J[1][i][j-1] )/sin(theta_m);
                Jtheta_av_c = 0.5*( sin(dm.theta[j])*J[1][i][j] - sin(dm.theta[j-1])*J[1][i][j-1] )/sin(theta_m);
                Bphi_phi = 0.25*( dm.r[i]*( sin(dm.theta[j])*Bn[2][i][j] + sin(dm.theta[j-1])*Bn[2][i][j-1] )
                                    + dm.r[i-1]*( sin(dm.theta[j])*Bn[2][i-1][j] + sin(dm.theta[j-1])*Bn[2][i-1][j+1] ) )/(r_m*sin(theta_m));
            }

            //Quantities at location of J_phi
            eta_H_corner = 2.*( dm.Deltar[i]+dm.Deltar[i-1] )/( dm.Deltar[i]*( 1./tC.eta_H[i][j] + 1./tC.eta_H[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_H[i-1][j] + 1./tC.eta_H[i-1][j-1] ) );
            Jr_av_phi = 0.5*( dm.r[i]*J[0][i][j] + dm.r[i-1]*J[0][i-1][j] )/r_m;
            Jr_av_c = 0.5*( dm.r[i]*J[0][i][j] - dm.r[i-1]*J[0][i-1][j] )/r_m;
            vr_av_phi = -eta_H_corner*Jr_av_phi;
            vtheta_av_phi = -eta_H_corner*Jtheta_av_phi;
            vr_av_c = -eta_H_corner*Jr_av_c;
            vtheta_av_c = -eta_H_corner*Jtheta_av_c;
            vphi_av_phi = -eta_H_corner*J[2][i][j];
            vphi_av_c = 0.;

            Btheta_up_phi = Reconstruct_Btheta(B, i, j, vr_av_phi, Deltar_1, Deltar_2, Deltar_3, dm);
            Br_up_phi = Reconstruct_Br(B, i, j, vtheta_av_phi, Deltatheta_1, Deltatheta_2, Deltatheta_3, dm);

            mC.phE_H[0][i][j] = -( vtheta_av_phi*Bphi_phi - 0*vphi_av_phi*Btheta_up_phi );
            mC.phE_H[1][i][j] = -( 0*vphi_av_phi*Br_up_phi - vr_av_phi*Bphi_phi );

            if(tC.conductivity_anisotropy == true){
                if(mC.quantization == true){
                    eta_O_corner_perp = 2.*(dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]*( 1./tC.eta_O_perp[i][j] + 1./tC.eta_O_perp[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_O_perp[i-1][j] + 1./tC.eta_O_perp[i-1][j-1] ) );
                    eta_O_corner_delta = 2.*(dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]*( 1./tC.eta_O_delta[i][j] + 1./tC.eta_O_delta[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_O_delta[i-1][j] + 1./tC.eta_O_delta[i-1][j-1] ) );

                    JdotBhat_phi = Jr_av_phi*mC.phBhat[0][i][j] + Jtheta_av_phi*mC.phBhat[1][i][j] + J[2][i][j]*mC.phBhat[2][i][j];
                    phE[0][i][j] = eta_O_corner_delta*JdotBhat_phi*mC.phBhat[0][i][j] + eta_O_corner_perp*Jr_av_phi - ( vtheta_av_phi*Bphi_phi - vphi_av_phi*Btheta_up_phi );
                    phE[1][i][j] = eta_O_corner_delta*JdotBhat_phi*mC.phBhat[1][i][j] + eta_O_corner_perp*Jtheta_av_phi - ( vphi_av_phi*Br_up_phi - vr_av_phi*Bphi_phi );
                    phE[2][i][j] = eta_O_corner_delta*JdotBhat_phi*mC.phBhat[2][i][j] + eta_O_corner_perp*J[2][i][j] - ( vr_av_phi*Btheta_up_phi - vtheta_av_phi*Br_up_phi );

                    JdotBhat_c = Jr_av_c*mC.phBhat[0][i][j] + Jtheta_av_c*mC.phBhat[1][i][j];
                    cE[0][i][j] = eta_O_corner_delta*JdotBhat_c*mC.phBhat[0][i][j] + eta_O_corner_perp*Jr_av_c - 0*( vtheta_av_c*Bphi_phi - vphi_av_c*Btheta_up_phi );
                    cE[1][i][j] = eta_O_corner_delta*JdotBhat_c*mC.phBhat[1][i][j] + eta_O_corner_perp*Jtheta_av_c - 0*( vphi_av_c*Br_up_phi - vr_av_c*Bphi_phi );
                    cE[2][i][j] = eta_O_corner_delta*JdotBhat_c*mC.phBhat[2][i][j] - 0*( vr_av_c*Btheta_up_phi - vtheta_av_c*Br_up_phi );
                }
                else{
                    eta_O_corner_perp = 2.*(dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]*( 1./tC.eta_O_perp[i][j] + 1./tC.eta_O_perp[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_O_perp[i-1][j] + 1./tC.eta_O_perp[i-1][j-1] ) );

                    phE[0][i][j] = eta_O_corner_perp*Jr_av_phi - ( vtheta_av_phi*Bphi_phi - vphi_av_phi*Btheta_up_phi );
                    phE[1][i][j] = eta_O_corner_perp*Jtheta_av_phi - ( vphi_av_phi*Br_up_phi - vr_av_phi*Bphi_phi );
                    phE[2][i][j] = eta_O_corner_perp*J[2][i][j] - ( vr_av_phi*Btheta_up_phi - vtheta_av_phi*Br_up_phi );

                    cE[0][i][j] = eta_O_corner_perp*Jr_av_c - ( vtheta_av_c*Bphi_phi - vphi_av_c*Btheta_up_phi );
                    cE[1][i][j] = eta_O_corner_perp*Jtheta_av_c - ( vphi_av_c*Br_up_phi - vr_av_c*Bphi_phi );
                    cE[2][i][j] = - ( vr_av_c*Btheta_up_phi - vtheta_av_c*Br_up_phi );
                }
            }
            else {
                eta_O_corner = 2.*(dm.Deltar[i]+dm.Deltar[i-1])/( dm.Deltar[i]*( 1./tC.eta_O[i][j] + 1./tC.eta_O[i][j-1] ) + dm.Deltar[i-1]*( 1./tC.eta_O[i-1][j] + 1./tC.eta_O[i-1][j-1] ) );

                phE[0][i][j] = eta_O_corner*Jr_av_phi - ( vtheta_av_phi*Bphi_phi - vphi_av_phi*Btheta_up_phi );
                phE[1][i][j] = eta_O_corner*Jtheta_av_phi - ( vphi_av_phi*Br_up_phi - vr_av_phi*Bphi_phi );
                phE[2][i][j] = eta_O_corner*J[2][i][j] - ( vr_av_phi*Btheta_up_phi - vtheta_av_phi*Br_up_phi );

                cE[0][i][j] = eta_O_corner*Jr_av_c - ( vtheta_av_c*Bphi_phi - vphi_av_c*Btheta_up_phi );
                cE[1][i][j] = eta_O_corner*Jtheta_av_c - ( vphi_av_c*Br_up_phi - vr_av_c*Bphi_phi );
                cE[2][i][j] = - ( vr_av_c*Btheta_up_phi - vtheta_av_c*Br_up_phi );
            }
        }
    }

    return;
}

/*
    Computes electromagnetic force along closed path around cell face. Used to compute time derivative of magnetic field.
    Inputs: phE: redshifted electric field times c in reduced units evaluated at cell corners (where J_phi lives)
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
    Output: Qr, Qtheta: line integral of E divided by cell face area on bottom cell faces (to update Br) and left cell faces (to update Btheta) respectively
*/
void Compute_EMF(ScalarField & Qr, ScalarField & Qtheta, VectorField & phE, size_t N_GC, const Domain & dm, const Process & process)
{
    double Deltatheta = dm.Deltatheta;

    for(size_t i=N_GC; i<phE.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<phE.shape()[2]-N_GC; j++){

            if( i < phE.shape()[1]-N_GC-1 ){
                Qr[i][j] = -( ( 1./tan(dm.theta[j]) + 1./tan(Deltatheta/2.) )*phE[2][i][j+1] - ( -1./tan(dm.theta[j]) + 1./tan(Deltatheta/2.) )*phE[2][i][j] )/( 2.*(dm.r[i]-0.5*dm.Deltar[i]) );
                Qtheta[i][j] = -( ( 1. - 0.5*dm.Deltar[i]/dm.r[i] )*phE[2][i][j] - ( 1. + 0.5*dm.Deltar[i]/dm.r[i] )*phE[2][i+1][j] )/(dm.n_r[i]*dm.Deltar[i]);
            }
            //Because r[phE.shape()[1]-N_GC-1] = r[phE.shape()[1]-N_GC-2], change dm.r[i]-0.5*dm.Deltar[i] to dm.r[i]+0.5*dm.Deltar[i]
            else{
                Qr[i][j] = -( ( 1./tan(dm.theta[j]) + 1./tan(dm.Deltatheta/2.) )*phE[2][i][j+1] - ( -1./tan(dm.theta[j]) + 1./tan(dm.Deltatheta/2.) )*phE[2][i][j] )/( 2.*(dm.r[i]+0.5*dm.Deltar[i]) );
                // Qtheta[i][j] = -( ( 1. - 0.5*dm.Deltar[i]/ri )*phE[2][i][j] - ( 1. + 0.5*dm.Deltar[i]/ri )*phE[2][i+1][j] )/(dm.n_r[i]*dm.Deltar[i]); //this is never used
            }
        }
    }

    exchng2Scalar(Qr,N_GC,process,dm.stridetype_Sca);
    exchng2Scalar(Qtheta,N_GC,process,dm.stridetype_Sca);

    return;
}

/*
    Computes heat flux density in reduced units for both isotropic and anisotropic heat conductivities.

    Inputs: T: redshifted temperature in reduced units
            trC: TransCoeffs object containing thermal conductivity kappa
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process objects containing information about the simulation domain and the current process
    Output: q: heat flux density times lapse function squared in reduced units
*/
void Compute_q(ScalarField & T, VectorField & q, TransCoeffs & trC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process)
{
    double r_iph, r_imh; //radial coordinate at upper and lower cell face
    double theta_jph, theta_jmh; //theta coordinate at right and left cell face
    double n_r_imh, n_r_iph; //sqrt(g_rr)=n_r at r = r_{i-1/2}, r_{i+1/2} respectively
    double Deltatheta = dm.Deltatheta;
    double gradrT_imhjm1, gradrT_iphjm1, gradrT_imhj, gradrT_iphj;
    double gradrT_ijm1, gradrT_ij, gradrT_ijmh;
    double gradthT_im1jmh, gradthT_im1jph, gradthT_ijmh, gradthT_ijph;
    double gradthT_im1j, gradthT_ij, gradthT_imhj;

    if(trC.conductivity_anisotropy == true){
        double kappa_delta_imh, kappa_perp_imh, kappa_H_imh; //kappa values computed at lower radial cell boundaries (index i-1/2) using harmonic mean
        double kappa_delta_jmh, kappa_perp_jmh, kappa_H_jmh; //kappa values computed at left theta cell boundaries (index j-1/2) using harmonic mean

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i];
            r_iph = dm.r[i] + 0.5*dm.Deltar[i];
            n_r_imh = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            n_r_iph = 0.5*( dm.n_r[i] + dm.n_r[i+1] );
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){
                theta_jmh = dm.theta[j] - 0.5*Deltatheta;
                theta_jph = dm.theta[j] + 0.5*Deltatheta;

                gradrT_imhjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i,j-1)/n_r_imh;
                gradrT_iphjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i+1,j-1)/n_r_iph;
                gradrT_imhj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i,j)/n_r_imh;
                gradrT_iphj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i+1,j)/n_r_iph;

                gradrT_ijm1 = MC( 0.5*(gradrT_imhjm1*r_imh+gradrT_iphjm1*r_iph), 2.*gradrT_imhjm1*r_imh, 2.*gradrT_iphjm1*r_iph )/dm.r[i];
                gradrT_ij = MC( 0.5*(gradrT_imhj*r_imh+gradrT_iphj*r_iph), 2.*gradrT_imhj*r_imh, 2.*gradrT_iphj*r_iph )/dm.r[i];

                if( theta_jmh < 1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else if( theta_jmh > pi-1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else{
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1*sin(dm.theta[j-1])+gradrT_ij*sin(dm.theta[j])), 2.*gradrT_ijm1*sin(dm.theta[j-1]), 2.*gradrT_ij*sin(dm.theta[j]) )/sin(theta_jmh);
                }
                gradthT_im1jmh = gradth(T,Deltatheta*dm.r[i-1],i-1,j);
                gradthT_im1jph = gradth(T,Deltatheta*dm.r[i-1],i-1,j+1);
                gradthT_ijmh = gradth(T,Deltatheta*dm.r[i],i,j);
                gradthT_ijph = gradth(T,Deltatheta*dm.r[i],i,j+1);

                gradthT_im1j = MC( 0.5*(gradthT_im1jmh*sin(theta_jmh)+gradthT_im1jph*sin(theta_jph)), 2.*gradthT_im1jmh*sin(theta_jmh), 2.*gradthT_im1jph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_ij = MC( 0.5*(gradthT_ijmh*sin(theta_jmh)+gradthT_ijph*sin(theta_jph)), 2.*gradthT_ijmh*sin(theta_jmh), 2.*gradthT_ijph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_imhj = MC( 0.5*(gradthT_im1j*dm.r[i-1]+gradthT_ij*dm.r[i]), 2.*gradthT_im1j*dm.r[i-1], 2.*gradthT_ij*dm.r[i] )/r_imh;

                kappa_delta_jmh = 2./( 1./trC.kappa_delta[i][j-1] + 1./trC.kappa_delta[i][j] );
                kappa_perp_jmh = 2./( 1./trC.kappa_perp[i][j-1] + 1./trC.kappa_perp[i][j] );
                kappa_H_jmh = 2./( 1./trC.kappa_H[i][j-1] + 1./trC.kappa_H[i][j] );

                kappa_delta_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_delta[i-1][j] + dm.Deltar[i]/trC.kappa_delta[i][j] );
                kappa_perp_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_perp[i-1][j] + dm.Deltar[i]/trC.kappa_perp[i][j] );
                kappa_H_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_H[i-1][j] + dm.Deltar[i]/trC.kappa_H[i][j] );

                q[0][i][j] = -0.5*(dm.n_t[i-1]+dm.n_t[i])*( kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[0][i][j] + kappa_perp_imh )*gradrT_imhj
                            - 0.5*(dm.n_t[i-1]+dm.n_t[i])*( kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[1][i][j] - kappa_H_imh*mC.thBhat[2][i][j] )*
                                                            0.5/tan(Deltatheta/2.)*Deltatheta*gradthT_imhj;
                q[1][i][j] = -dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[1][i][j] + kappa_perp_jmh )*gradthT_ijmh
                            - dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[0][i][j] + kappa_H_jmh*mC.rBhat[2][i][j] )*gradrT_ijmh;
            }
        }
    }
    else{
        double kappa_imh, kappa_jmh;

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            n_r_imh = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){
                kappa_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa[i-1][j] + dm.Deltar[i]/trC.kappa[i][j] );
                q[0][i][j] = -0.5*(dm.n_t[i-1]+dm.n_t[i])*kappa_imh*gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i,j)/n_r_imh;
                if( j == 0 ){
                    q[1][i][j] = 0.; //sets heat flux in polar direction to zero at left edge of domain. If this edge is not at theta = 0, heat flux will be replaced by exchng2Scalar
                }
                else{
                    kappa_jmh = 2./( 1./trC.kappa[i][j-1] + 1./trC.kappa[i][j] );
                    q[1][i][j] = -dm.n_t[i]*kappa_jmh*gradth(T,dm.Deltatheta*dm.r[i],i,j);
                }
            }
        }
    }

    exchng2Vector(q, N_GC, process, dm.stridetype_Vec_q);

    return;
}

/*
    Computes the update for the temperature for both isotropic and anisotropic heat conductivities
    Inputs: T, q, J, rE, thE, phE: redshifted temperature, heat flux density, redshifted current density and redshifted electric field times c in reduced units (evaluated at locations of J_r, J_theta, J_phi)
            trC: TransCoeffs object containing thermal conductivity kappa
            thC: ThermCoeffs object containing specific heat capacity c_v and neutrino emissivity q_nu
            tparams: TParams object containing information about temperature evolution
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
    Output: QT: update for T
*/
void TEvolve(ScalarField & QT, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH, TransCoeffs & trC, ThermCoeffs & thC, TParams & tparams, size_t N_GC, const Domain & dm)
{
    double Fr, Ftheta; //flux difference function in r and theta directions and cross term
    double r_iph, r_imh; //radial coordinate at upper and lower cell face
    double theta_jph, theta_jmh; //theta coordinate at right and left cell face
    static double JH_pref = tparams.JouleHeating*B_0*B_0/(4.*pi*T_0*s_0); //prefactor to de-dimensionalize Joule heating term
    double JH; //Joule heating term without JH_pref.

    for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
        r_iph = dm.r[i] + 0.5*dm.Deltar[i];
        r_imh = dm.r[i] - 0.5*dm.Deltar[i];
        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            theta_jph = dm.theta[j] + 0.5*dm.Deltatheta;
            theta_jmh = dm.theta[j] - 0.5*dm.Deltatheta;

            Fr = r_iph*r_iph*q[0][i+1][j] - r_imh*r_imh*q[0][i][j];
            Ftheta = sin(theta_jph)*q[1][i][j+1] - sin(theta_jmh)*q[1][i][j];

            JH = JouleHeating(i,j,J,phE,cE,trC,dm);

            //Divide by overall factor of c_v[i][j]
            QT[i][j] = ( -Fr/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*2.*sin(dm.theta[j])*sin(dm.Deltatheta/2.))
                        + JH_pref*JH + 0*q_SH[i][j] + dm.n_t[i]*dm.n_t[i]*thC.q_nuCrust[i][j] )/thC.c_v[i][j];
        }
    }

    return;
}

/*
    Computes the Joule heating density for cell (i,j)
    Inputs: i, j: indices of cell for which to compute Joule heating rate
            J, phE, cE: redshifted current density, and redshifted electric field and conjugate electric field times c in reduced units (evaluated at locations of J_phi)
            dm: Domain object containing information about the simulation domain
    Output: Joule heating density of cell (i,j) in reduced units
*/
double JouleHeating(size_t i, size_t j, VectorField & J, VectorField & phE, VectorField & cE, TransCoeffs & tC, const Domain & dm)
{
    double JthEth_ph_ij, JthEth_ph_ip1j, JthEth_ph_ijp1, JthEth_ph_ip1jp1;
    double JthEth_c_ij, JthEth_c_ip1j, JthEth_c_ijp1, JthEth_c_ip1jp1;
    double JH_ph, JH_c; //Joule heating density contributions at location of currents in r, theta, phi-directions respectively

    double JH_ph_lb, JH_ph_rb, JH_ph_lt, JH_ph_rt;
    double JH_c_lb, JH_c_rb, JH_c_lt, JH_c_rt;

    double r_m = dm.r[i] - 0.5*dm.Deltar[i], r_p = dm.r[i] + 0.5*dm.Deltar[i];
    double theta_m = dm.theta[j] - 0.5*dm.Deltatheta, theta_p = dm.theta[j] + 0.5*dm.Deltatheta;

    double rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
    double Vij = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*dm.Deltatheta)*sin(dm.theta[j])*dm.n_r[i];
    double r_int_b = dm.n_r[i]*( dm.r[i]*dm.r[i]*dm.Deltar[i]/2.*rJacobian - dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/6. );
    double r_int_t = dm.n_r[i]*( dm.r[i]*dm.r[i]*dm.Deltar[i]/2.*rJacobian + dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/6. );
    double theta_int_l = cos( dm.theta[j] - 0.5*dm.Deltatheta ) - 2.*cos(dm.theta[j])*sin(0.5*dm.Deltatheta)/dm.Deltatheta;
    double theta_int_r = -cos( dm.theta[j] + 0.5*dm.Deltatheta ) + 2.*cos(dm.theta[j])*sin(0.5*dm.Deltatheta)/dm.Deltatheta;
    double lb_int = 2.*pi*r_int_b*theta_int_l, lt_int = 2.*pi*r_int_t*theta_int_l, rb_int = 2.*pi*r_int_b*theta_int_r, rt_int = 2.*pi*r_int_t*theta_int_r;

    if( dm.theta[j]-0.5*dm.Deltatheta > 0. ){
        JthEth_ph_ij = 0.5*( sin(dm.theta[j-1])*J[1][i][j-1] + sin(dm.theta[j])*J[1][i][j] )/sin(theta_m)*phE[1][i][j];
        JthEth_ph_ip1j= 0.5*( sin(dm.theta[j-1])*J[1][i+1][j-1] + sin(dm.theta[j])*J[1][i+1][j] )/sin(theta_m)*phE[1][i+1][j];
        JthEth_c_ij = 0.5*( sin(dm.theta[j])*J[1][i][j] - sin(dm.theta[j-1])*J[1][i][j-1] )/sin(theta_m)*cE[1][i][j];
        JthEth_c_ip1j = 0.5*( sin(dm.theta[j])*J[1][i+1][j] - sin(dm.theta[j-1])*J[1][i+1][j-1] )/sin(theta_m)*cE[1][i+1][j];
        // JthEth_ph_ij = 0.5*( J[1][i][j-1] + J[1][i][j] )*phE[1][i][j];
        // JthEth_ph_ip1j= 0.5*( J[1][i+1][j-1] + J[1][i+1][j] )*phE[1][i+1][j];
        // JthEth_c_ij = 0.5*( J[1][i][j] - J[1][i][j-1] )*cE[1][i][j];
        // JthEth_c_ip1j = 0.5*( J[1][i+1][j] - J[1][i+1][j-1] )*cE[1][i+1][j];
    }
    else{
        JthEth_ph_ij = 0.;
        JthEth_ph_ip1j = 0.;
        JthEth_c_ij = 0.;
        JthEth_c_ip1j = 0.;
    }

    if( dm.theta[j]+0.5*dm.Deltatheta + 1e-9 < pi ){
        JthEth_ph_ijp1 = 0.5*( sin(dm.theta[j])*J[1][i][j] + sin(dm.theta[j+1])*J[1][i][j+1] )/sin(theta_p)*phE[1][i][j+1];
        JthEth_ph_ip1jp1 = 0.5*( sin(dm.theta[j])*J[1][i+1][j] + sin(dm.theta[j+1])*J[1][i+1][j+1] )/sin(theta_p)*phE[1][i+1][j+1];
        JthEth_c_ijp1 = 0.5*( sin(dm.theta[j+1])*J[1][i][j+1] - sin(dm.theta[j])*J[1][i][j] )/sin(theta_p)*cE[1][i][j+1];
        JthEth_c_ip1jp1 = 0.5*( sin(dm.theta[j+1])*J[1][i+1][j+1] - sin(dm.theta[j])*J[1][i+1][j] )/sin(theta_p)*cE[1][i+1][j+1];
        // JthEth_ph_ijp1 = 0.5*( J[1][i][j] + J[1][i][j+1] )*phE[1][i][j+1];
        // JthEth_ph_ip1jp1 = 0.5*( J[1][i+1][j] + J[1][i+1][j+1] )*phE[1][i+1][j+1];
        // JthEth_c_ijp1 = 0.5*( J[1][i][j+1] - J[1][i][j] )*cE[1][i][j+1];
        // JthEth_c_ip1jp1 = 0.5*( J[1][i+1][j+1] - J[1][i+1][j] )*cE[1][i+1][j+1];
    }
    else{
        JthEth_ph_ijp1 = 0.;
        JthEth_ph_ip1jp1 = 0.;
        JthEth_c_ijp1 = 0.;
        JthEth_c_ip1jp1 = 0.;
    }

    JH_ph_lb = 0.5*( dm.r[i-1]*J[0][i-1][j] + dm.r[i]*J[0][i][j] )/r_m*phE[0][i][j] + JthEth_ph_ij + J[2][i][j]*phE[2][i][j];
    JH_ph_rb = 0.5*( dm.r[i-1]*J[0][i-1][j+1] + dm.r[i]*J[0][i][j+1] )/r_m*phE[0][i][j+1] + JthEth_ph_ijp1 + J[2][i][j+1]*phE[2][i][j+1];
    JH_ph_lt = 0.5*( dm.r[i]*J[0][i][j] + dm.r[i+1]*J[0][i+1][j] )/r_p*phE[0][i+1][j] + JthEth_ph_ip1j + J[2][i+1][j]*phE[2][i+1][j];
    JH_ph_rt = 0.5*( dm.r[i]*J[0][i][j+1] + dm.r[i+1]*J[0][i+1][j+1] )/r_p*phE[0][i+1][j+1] + JthEth_ph_ip1jp1 + J[2][i+1][j+1]*phE[2][i+1][j+1];

    // JH_ph = ( JH_ph_lb + JH_ph_rb + JH_ph_lt + JH_ph_rt )/4.;
    JH_ph = ( JH_ph_lb*lb_int + JH_ph_rb*rb_int + JH_ph_lt*lt_int + JH_ph_rt*rt_int )/Vij;

    JH_c_lb = 0.5*( dm.r[i]*J[0][i][j] - dm.r[i-1]*J[0][i-1][j] )/r_m*cE[0][i][j] + JthEth_c_ij;
    JH_c_rb = 0.5*( dm.r[i]*J[0][i][j+1] - dm.r[i-1]*J[0][i-1][j+1] )/r_m*cE[0][i][j+1] + JthEth_c_ijp1;
    JH_c_lt = 0.5*( dm.r[i+1]*J[0][i+1][j] - dm.r[i]*J[0][i][j] )/r_p*cE[0][i+1][j] + JthEth_c_ip1j;
    JH_c_rt = 0.5*( dm.r[i+1]*J[0][i+1][j+1] - dm.r[i]*J[0][i][j+1] )/r_p*cE[0][i+1][j+1] + JthEth_c_ip1jp1;

    // JH_c = ( JH_c_lb + JH_c_rb + JH_c_lt + JH_c_rt )/4.;
    JH_c = ( JH_c_lb*lb_int + JH_c_rb*rb_int + JH_c_lt*lt_int + JH_c_rt*rt_int )/Vij;

    return JH_ph + JH_c;
}

/*
    Computes the quasi-Joule heating density (rho.J.J_B) for cell (i,j)
    Inputs: i, j: indices of cell for which to compute Joule heating rate
            J, phE, cE: redshifted current density, and redshifted electric field and conjugate electric field times c in reduced units (evaluated at locations of J_phi)
            dm: Domain object containing information about the simulation domain
    Output: quasi-Joule heating density of cell (i,j) in reduced units
*/
double QuasiJouleHeating(size_t i, size_t j, VectorField & J_B, VectorField & phE, VectorField & cE, TransCoeffs & tC, MagCoeffs & mC, const Domain & dm)
{
    double JthEth_ph_ij, JthEth_ph_ip1j, JthEth_ph_ijp1, JthEth_ph_ip1jp1; //J_theta*E_theta evaluated at locations of J_phi
    double JthEthH_ph_ij, JthEthH_ph_ip1j, JthEthH_ph_ijp1, JthEthH_ph_ip1jp1; //J_theta*E_theta_H evaluated at locations of J_phi
    double JthEth_c_ij, JthEth_c_ip1j, JthEth_c_ijp1, JthEth_c_ip1jp1; //conjugate J_theta*E_theta evaluated at locations of J_phi
    double QuasiJH_ph, QuasiJH_c; //quasi-Joule heating density contributions at location of currents in r, theta, phi-directions respectively

    double JrErH_ph_ij, JrErH_ph_ip1j, JrErH_ph_ijp1, JrErH_ph_ip1jp1; //J_r*E_r_H evaluated at locations of J_phi
    double JH_ph_lb, JH_ph_rb, JH_ph_lt, JH_ph_rt;
    double JH_c_lb, JH_c_rb, JH_c_lt, JH_c_rt;
    double JrEr_l, JrEr_r; //difference between J_r*E_r_H at location of J_r and the sum of half its value at each of the neighbouring J_phi locations
    double JthEth_b, JthEth_t; //difference between J_theta*E_theta_H at location of J_theta and the sum of half its value at each of the neighbouring J_phi locations

    double r_m = dm.r[i] - 0.5*dm.Deltar[i], r_p = dm.r[i] + 0.5*dm.Deltar[i];
    double theta_m = dm.theta[j] - 0.5*dm.Deltatheta, theta_p = dm.theta[j] + 0.5*dm.Deltatheta;

    double rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
    double Vij = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*dm.Deltatheta)*sin(dm.theta[j])*dm.n_r[i];
    double r_int_b = dm.n_r[i]*( dm.r[i]*dm.r[i]*dm.Deltar[i]/2.*rJacobian - dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/6. );
    double r_int_t = dm.n_r[i]*( dm.r[i]*dm.r[i]*dm.Deltar[i]/2.*rJacobian + dm.r[i]*dm.Deltar[i]*dm.Deltar[i]/6. );
    double theta_int_l = cos( dm.theta[j] - 0.5*dm.Deltatheta ) - 2.*cos(dm.theta[j])*sin(0.5*dm.Deltatheta)/dm.Deltatheta;
    double theta_int_r = -cos( dm.theta[j] + 0.5*dm.Deltatheta ) + 2.*cos(dm.theta[j])*sin(0.5*dm.Deltatheta)/dm.Deltatheta;
    double lb_int = 2.*pi*r_int_b*theta_int_l, lt_int = 2.*pi*r_int_t*theta_int_l, rb_int = 2.*pi*r_int_b*theta_int_r, rt_int = 2.*pi*r_int_t*theta_int_r;

    if( dm.theta[j]-0.5*dm.Deltatheta > 0. ){
        JthEth_ph_ij = 0.5*( sin(dm.theta[j-1])*J_B[1][i][j-1] + sin(dm.theta[j])*J_B[1][i][j] )/sin(theta_m)*phE[1][i][j];
        JthEth_ph_ip1j= 0.5*( sin(dm.theta[j-1])*J_B[1][i+1][j-1] + sin(dm.theta[j])*J_B[1][i+1][j] )/sin(theta_m)*phE[1][i+1][j];
        JthEth_c_ij = 0.5*( sin(dm.theta[j])*J_B[1][i][j] - sin(dm.theta[j-1])*J_B[1][i][j-1] )/sin(theta_m)*cE[1][i][j];
        JthEth_c_ip1j = 0.5*( sin(dm.theta[j])*J_B[1][i+1][j] - sin(dm.theta[j-1])*J_B[1][i+1][j-1] )/sin(theta_m)*cE[1][i+1][j];
        JthEthH_ph_ij = 0.5*( sin(dm.theta[j-1])*J_B[1][i][j-1] + sin(dm.theta[j])*J_B[1][i][j] )/sin(theta_m)*mC.phE_H[1][i][j];
        JthEthH_ph_ip1j = 0.5*( sin(dm.theta[j-1])*J_B[1][i+1][j-1] + sin(dm.theta[j])*J_B[1][i+1][j] )/sin(theta_m)*mC.phE_H[1][i+1][j];
    }
    else{
        JthEth_ph_ij = 0.;
        JthEth_ph_ip1j = 0.;
        JthEth_c_ij = 0.;
        JthEth_c_ip1j = 0.;
        JthEthH_ph_ij = 0.;
        JthEthH_ph_ip1j = 0.;
    }

    if( dm.theta[j]+0.5*dm.Deltatheta + 1e-9 < pi ){
        JthEth_ph_ijp1 = 0.5*( sin(dm.theta[j])*J_B[1][i][j] + sin(dm.theta[j+1])*J_B[1][i][j+1] )/sin(theta_p)*phE[1][i][j+1];
        JthEth_ph_ip1jp1 = 0.5*( sin(dm.theta[j])*J_B[1][i+1][j] + sin(dm.theta[j+1])*J_B[1][i+1][j+1]) /sin(theta_p)*phE[1][i+1][j+1];
        JthEth_c_ijp1 = 0.5*( sin(dm.theta[j+1])*J_B[1][i][j+1] - sin(dm.theta[j])*J_B[1][i][j] )/sin(theta_p)*cE[1][i][j+1];
        JthEth_c_ip1jp1 = 0.5*( sin(dm.theta[j+1])*J_B[1][i+1][j+1] - sin(dm.theta[j])*J_B[1][i+1][j] )/sin(theta_p)*cE[1][i+1][j+1];
        JthEthH_ph_ijp1 = 0.5*( sin(dm.theta[j])*J_B[1][i][j] + sin(dm.theta[j+1])*J_B[1][i][j+1] )/sin(theta_p)*mC.phE_H[1][i][j+1];
        JthEthH_ph_ip1jp1 = 0.5*( sin(dm.theta[j])*J_B[1][i+1][j] + sin(dm.theta[j+1])*J_B[1][i+1][j+1] )/sin(theta_p)*mC.phE_H[1][i+1][j+1];
    }
    else{
        JthEth_ph_ijp1 = 0.;
        JthEth_ph_ip1jp1 = 0.;
        JthEth_c_ijp1 = 0.;
        JthEth_c_ip1jp1 = 0.;
        JthEthH_ph_ijp1 = 0.;
        JthEthH_ph_ip1jp1 = 0.;
    }

    JH_ph_lb = 0.5*( dm.r[i-1]*J_B[0][i-1][j] + dm.r[i]*J_B[0][i][j] )/r_m*phE[0][i][j] + JthEth_ph_ij + J_B[2][i][j]*phE[2][i][j];
    JH_ph_rb = 0.5*( dm.r[i-1]*J_B[0][i-1][j+1] + dm.r[i]*J_B[0][i][j+1] )/r_m*phE[0][i][j+1] + JthEth_ph_ijp1 + J_B[2][i][j+1]*phE[2][i][j+1];
    JH_ph_lt = 0.5*( dm.r[i]*J_B[0][i][j] + dm.r[i+1]*J_B[0][i+1][j] )/r_p*phE[0][i+1][j] + JthEth_ph_ip1j + J_B[2][i+1][j]*phE[2][i+1][j];
    JH_ph_rt = 0.5*( dm.r[i]*J_B[0][i][j+1] + dm.r[i+1]*J_B[0][i+1][j+1] )/r_p*phE[0][i+1][j+1] + JthEth_ph_ip1jp1 + J_B[2][i+1][j+1]*phE[2][i+1][j+1];

    JrErH_ph_ij = 0.5*( dm.r[i-1]*J_B[0][i-1][j] + dm.r[i]*J_B[0][i][j] )/r_m*mC.phE_H[0][i][j];
    JrErH_ph_ijp1 = 0.5*( dm.r[i-1]*J_B[0][i-1][j+1] + dm.r[i]*J_B[0][i][j+1] )/r_m*mC.phE_H[0][i][j+1];
    JrErH_ph_ip1j = 0.5*( dm.r[i]*J_B[0][i][j] + dm.r[i+1]*J_B[0][i+1][j] )/r_p*mC.phE_H[0][i+1][j];
    JrErH_ph_ip1jp1 = 0.5*( dm.r[i]*J_B[0][i][j+1] + dm.r[i+1]*J_B[0][i+1][j+1] )/r_p*mC.phE_H[0][i+1][j+1];

    JrEr_l = J_B[0][i][j]*mC.E_H_Pol[0][i][j] - 0.5*( JrErH_ph_ij + JrErH_ph_ip1j );
    JrEr_r = J_B[0][i][j+1]*mC.E_H_Pol[0][i][j+1] - 0.5*( JrErH_ph_ijp1 + JrErH_ph_ip1jp1 );

    JthEth_b = J_B[1][i][j]*mC.E_H_Pol[1][i][j] - 0.5*( JthEthH_ph_ij + JthEthH_ph_ijp1 );
    JthEth_t = J_B[1][i+1][j]*mC.E_H_Pol[1][i+1][j] - 0.5*( JthEthH_ph_ip1j + JthEthH_ph_ip1jp1 );

    // QuasiJH_ph = ( JH_ph_lb + JH_ph_rb + JH_ph_lt + JH_ph_rt )/4.;// + ( JrEr_l + JrEr_r + JthEth_b + JthEth_t )/2.;
    QuasiJH_ph = ( JH_ph_lb*lb_int + JH_ph_rb*rb_int + JH_ph_lt*lt_int + JH_ph_rt*rt_int )/Vij;// + ( JrEr_l + JrEr_r + JthEth_b + JthEth_t )/2.;

    JH_c_lb = 0.5*( dm.r[i]*J_B[0][i][j] - dm.r[i-1]*J_B[0][i-1][j] )/r_m*cE[0][i][j] + JthEth_c_ij;
    JH_c_rb = 0.5*( dm.r[i]*J_B[0][i][j+1] - dm.r[i-1]*J_B[0][i-1][j+1] )/r_m*cE[0][i][j+1] + JthEth_c_ijp1;
    JH_c_lt = 0.5*( dm.r[i+1]*J_B[0][i+1][j] - dm.r[i]*J_B[0][i][j] )/r_p*cE[0][i+1][j] + JthEth_c_ip1j;
    JH_c_rt = 0.5*( dm.r[i+1]*J_B[0][i+1][j+1] - dm.r[i]*J_B[0][i][j+1] )/r_p*cE[0][i+1][j+1] + JthEth_c_ip1jp1;

    // QuasiJH_c = ( JH_c_lb + JH_c_rb + JH_c_lt + JH_c_rt )/4.;
    QuasiJH_c = ( JH_c_lb*lb_int + JH_c_rb*rb_int + JH_c_lt*lt_int + JH_c_rt*rt_int )/Vij;

    return QuasiJH_ph + QuasiJH_c;
}

/*
    Computes update for the core temperature, which is assumed uniform.
    Inputs: T: redshifted temperature in reduced units
            Tcore: current core redshifted temperature in reduced units
            tC: TransCoeffs object containing thermal conductivity kappa
            N_GC: number of ghost cells
            dm, process: Domain and Process objects containing information about the simulation domain and the current process
            corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
    Output: QTore: update for Tcore
*/
double TCoreEvolve(ScalarField & T, double Tcore, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, const CoreThermalParams & corethermalparams)
{
    double C_v, Q_nu; //core volume-integrated specific heat capacity and neutrino emissivity in reduced units
    double kappa_crust, kappa_crust_H; //thermal conductivity for radial heat flow on crust side of crust-core transition.
    // For anisotropic conductivity, this must equal the perpendicular conductivity because Br vanishes so there is no heat flowing parallel to field lines from crust to core or vice versa

    double HeatFluxPartial = 0.; //Surface integral of heat flux at crust-core boundary divided by heat capacity of core, for current process only
    double r_inner = dm.r[N_GC]-0.5*dm.Deltar[N_GC]; //inner radius of simulation domain
    //Integrate heat flux over crust-core surface, a sphere of radius r[N_GC]-0.5*dm.Deltar[N_GC]. Take thermal conductivity as mean between values on either side of crust-core transition
    for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
        if( tC.conductivity_anisotropy == true ) {
            kappa_crust = 0.5*( tC.kappa_perp[N_GC][j] + tC.kappa_perp[N_GC-1][j] );
            kappa_crust_H = 0.5*( tC.kappa_H[N_GC][j] +  tC.kappa_H[N_GC-1][j] );
            HeatFluxPartial += -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust*2.*(T[N_GC][j]-Tcore)/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]))
                                +2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*dm.n_t[N_GC]*kappa_crust_H*mC.thBhat[2][N_GC][j]*0.5/tan(dm.Deltatheta/2.)*( 0.25*T[N_GC][j+1] - 0.25*T[N_GC][j-1] );
        }
        else {
            kappa_crust = 0.5*( tC.kappa[N_GC][j] + tC.kappa[N_GC-1][j] );//0.5*( tC.kappa[N_GC][j] + kappa_e + kappa_n/Tcore );
            HeatFluxPartial += -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust*2.*(T[N_GC][j]-Tcore)/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]));
        }
    }

    std::vector<double> HeatFluxVec(process.num_procs);
    double HeatFlux; //Surface integral of heat flux at crust-core boundary for entire domain
    double QTcore; //update for Tcore

    MPI_Allgather(&HeatFluxPartial, 1, MPI_DOUBLE, HeatFluxVec.data(), 1, MPI_DOUBLE, process.comm1D);

    if(process.world_rank == 0){
        HeatFlux = std::reduce(HeatFluxVec.begin(), HeatFluxVec.end()); //add HeatFluxPartials for all processes
        C_v = gsl_spline_eval( corethermalparams.C_vCore, Tcore, corethermalparams.C_vCore_acc );
        Q_nu = gsl_spline_eval( corethermalparams.Q_nuCore, Tcore, corethermalparams.Q_nuCore_acc );
        QTcore = Q_nu/C_v - HeatFlux/C_v; //add neutrino emissivity of entire core divided by heat capacity of entire core to HeatFlux to obtain update for Tcore
    }
    MPI_Bcast(&QTcore, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of QTcore to all processes

    return QTcore;
}

/*
    Computes IMEX-intermediate update for the core temperature, which is assumed uniform.
    Inputs: T: redshifted temperature in reduced units
            Tcore: current core redshifted temperature in reduced units
            tC: TransCoeffs object containing thermal conductivity kappa
            N_GC: number of ghost cells
            dm, process: Domain and Process objects containing information about the simulation domain and the current process
            corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
            a_ijIm: off-diagonal entry in Butcher tableau for this implicit RK step
    Output: Tcore_star: IMEX-intermediate redshifted core temperature in reduced units
*/
double TCore_starEvolve(ScalarField & T, double Tcore, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, const CoreThermalParams & corethermalparams, double a_ijIm)
{
    double C_v, Q_nu; //core volume-integrated specific heat capacity and neutrino emissivity in reduced units
    double kappa_crust, kappa_crust_H; //thermal conductivity for radial heat flow on crust side of crust-core transition.
    // For anisotropic conductivity, this must equal the perpendicular conductivity because Br vanishes so there is no heat flowing parallel to field lines from crust to core or vice versa

    double HeatFluxPartial = 0.; //Surface integral of heat flux at crust-core boundary divided by heat capacity of core, for current process only
    double r_inner = dm.r[N_GC]-0.5*dm.Deltar[N_GC]; //inner radius of simulation domain
    //Integrate heat flux over crust-core surface, a sphere of radius r[N_GC]-0.5*dm.Deltar[N_GC]. Take thermal conductivity as mean between values on either side of crust-core transition
    for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
        if( tC.conductivity_anisotropy == true ) {
            kappa_crust = 0.5*( tC.kappa_perp[N_GC][j] + tC.kappa_perp[N_GC-1][j] );
            kappa_crust_H = 0.5*( tC.kappa_H[N_GC][j] + tC.kappa_H[N_GC-1][j] );
            // std::cout << "j = " << j << ", kappa_crust = " << kappa_crust << ", kappa_par = " << kappa_crust + 0.5*( tC.kappa_delta[N_GC][j] +  tC.kappa_delta[N_GC-1][j] ) <<
            //             ", kappa_H = " << 0.5*( tC.kappa_H[N_GC][j] + tC.kappa_H[N_GC-1][j]) << ", Bmag = " << mC.Bmag[N_GC][j];
            // std::cout << ", Bhat_r = " << std::setprecision(8) << mC.thBhat[0][N_GC][j] << ", Bhat_th = " << mC.thBhat[1][N_GC][j] << ", Bhat_phi = " << mC.thBhat[2][N_GC][j] << std::endl;
            HeatFluxPartial += -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust*(T[N_GC][j]-T[N_GC-1][j])/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]))
                                +2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*dm.n_t[N_GC]*kappa_crust_H*mC.thBhat[2][N_GC][j]*0.5/tan(dm.Deltatheta/2.)*( 0.25*T[N_GC][j+1] - 0.25*T[N_GC][j-1] );
        }
        else {
            kappa_crust = 0.5*( tC.kappa[N_GC][j] + tC.kappa[N_GC-1][j] );
            HeatFluxPartial += -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust*(T[N_GC][j]-T[N_GC-1][j])/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]));
        }
    }

    std::vector<double> HeatFluxVec(process.num_procs);
    double HeatFlux; //Surface integral of heat flux at crust-core boundary for entire domain
    double QTcore; //update for Tcore

    MPI_Allgather(&HeatFluxPartial, 1, MPI_DOUBLE, HeatFluxVec.data(), 1, MPI_DOUBLE, process.comm1D);

    if(process.world_rank == 0){
        HeatFlux = std::reduce(HeatFluxVec.begin(), HeatFluxVec.end()); //add HeatFluxPartials for all processes
        C_v = gsl_spline_eval( corethermalparams.C_vCore, Tcore, corethermalparams.C_vCore_acc );
        Q_nu = gsl_spline_eval( corethermalparams.Q_nuCore, Tcore, corethermalparams.Q_nuCore_acc );
        QTcore = Q_nu/C_v - HeatFlux/C_v; //add neutrino emissivity of entire core divided by heat capacity of entire core to HeatFlux to obtain update for Tcore
    }
    MPI_Bcast(&QTcore, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of QTcore to all processes

    return Tcore + a_ijIm*dm.Deltat*QTcore;
}

/*
    Computes update for the core temperature, which is assumed uniform, using IMEX method
    Inputs: T: redshifted temperature in reduced units
            Tcore: current core redshifted temperature in reduced units
            Tcore_star: intermediate temperature for IMEX method
            tC: TransCoeffs object containing thermal conductivity kappa
            N_GC: number of ghost cells
            dm, process: Domain and Process objects containing information about the simulation domain and the current process
            corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
            a_ii: diagonal entry in Butcher tableau for this implicit RK step
    Output: Tcore_new: updated Tcore using implicit method
*/
double TCoreEvolveImplicit(ScalarField & T, double Tcore, double Tcore_star, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, CoreThermalParams & corethermalparams, double a_ii)
{
    double tempHFP, tempHF_H; //temporary variable to hold prefactor for kappa_perp part of heat flux and total kappa_H part of heat flux of each crust-core boundary cell
    double C_v, Q_nu; //core volume-integrated specific heat capacity and neutrino emissivity in reduced units
    double kappa_crust, kappa_crust_H; //thermal conductivities for radial heat flow on crust side of crust-core transition.
    // For anisotropic conductivity, there are contributions to the heat flow proportional to both kappa_perp and kappa_H

    double HeatFluxPartial = 0.; //Surface integral of heat flux at crust-core boundary, for current process only
    double HeatFluxPartial_pref = 0.; //Surface integral of heat flux proportional to kappa_perp divided by core temperature at crust-core boundary, for current process only
    double HeatFluxPartial_H = 0.; //Surface integral of thermal Hall part of heat flux across crust-core boundary, for current process only
    double r_inner = dm.r[N_GC]-0.5*dm.Deltar[N_GC]; //inner radius of simulation domain
    //Integrate heat flux over crust-core surface, a sphere of radius r[N_GC]-0.5*dm.Deltar[N_GC]. Take thermal conductivity as mean between values on either side of crust-core transition
    for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
        if( tC.conductivity_anisotropy == true ) {
            kappa_crust = 0.5*( tC.kappa_perp[N_GC][j] + tC.kappa_perp[N_GC-1][j] );
            kappa_crust_H = 0.5*( tC.kappa_H[N_GC][j] + tC.kappa_H[N_GC-1][j] );
            tempHFP = -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]));
            tempHF_H = 2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*dm.n_t[N_GC]*kappa_crust_H*mC.thBhat[2][N_GC][j]*0.5/tan(dm.Deltatheta/2.)*( 0.25*T[N_GC][j+1] - 0.25*T[N_GC][j-1] );
            HeatFluxPartial_pref += tempHFP;
            HeatFluxPartial_H += tempHF_H;
            // HeatFluxPartial += tempHFP*T[N_GC][j] + tempHF_H; //total heat flux is tempHFP*(T[N_GC][j]-Tcore) + tempHF_H
            HeatFluxPartial += tempHFP*( T[N_GC][j] - T[N_GC-1][j] ) + tempHF_H; //total heat flux is tempHFP*(T[N_GC][j]-Tcore) + tempHF_H

        }
        else {
            kappa_crust = 0.5*( tC.kappa[N_GC][j] + tC.kappa[N_GC-1][j] );
            // tempHFP = -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust*2./(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]));
            tempHFP = -2.*pi*sin(dm.theta[j])*dm.Deltatheta*r_inner*r_inner*dm.n_t[N_GC]*kappa_crust/(dm.n_r[N_GC]*0.5*(dm.Deltar[N_GC-1] + dm.Deltar[N_GC]));
            HeatFluxPartial_pref += tempHFP;
            // HeatFluxPartial += tempHFP*T[N_GC][j];
            // HeatFluxPartial += tempHFP*( T[N_GC][j] - Tcore );
            HeatFluxPartial += tempHFP*( T[N_GC][j] - T[N_GC-1][j] );
        }
    }

    std::vector<double> HeatFluxVec(process.num_procs);
    std::vector<double> HeatFlux_prefVec(process.num_procs);
    std::vector<double> HeatFlux_HVec(process.num_procs);
    double HeatFlux; //Surface integral of heat flux at crust-core boundary for entire domain.
    double HeatFlux_pref; //Prefactor for surface integral of heat flux at crust-core boundary for entire domain.
    double HeatFlux_H; //Surface integral of thermal Hall heat flux at crust-core boundary for entire domain.
    double dQ_nuCoredT, denom;
    double Tcore_new; //update for Tcore. Equal to the new T_core.

    MPI_Allgather(&HeatFluxPartial, 1, MPI_DOUBLE, HeatFluxVec.data(), 1, MPI_DOUBLE, process.comm1D);
    MPI_Allgather(&HeatFluxPartial_pref, 1, MPI_DOUBLE, HeatFlux_prefVec.data(), 1, MPI_DOUBLE, process.comm1D);
    MPI_Allgather(&HeatFluxPartial_H, 1, MPI_DOUBLE, HeatFlux_HVec.data(), 1, MPI_DOUBLE, process.comm1D);

    if(process.world_rank == 0){
        HeatFlux = std::reduce(HeatFluxVec.begin(), HeatFluxVec.end()); //add HeatFluxPartial for all processes
        HeatFlux_pref = std::reduce(HeatFlux_prefVec.begin(), HeatFlux_prefVec.end()); //add HeatFluxPartial_pref for all processes
        HeatFlux_H = std::reduce(HeatFlux_HVec.begin(), HeatFlux_HVec.end()); //add HeatFluxPartial_H for all processes

        C_v = gsl_spline_eval( corethermalparams.C_vCore, Tcore, corethermalparams.C_vCore_acc );
        Q_nu = gsl_spline_eval( corethermalparams.Q_nuCore, Tcore, corethermalparams.Q_nuCore_acc );

        //Compute finite-difference approximation to dQ_nu/dT at T=Tcore
        double DeltaT = 1e-4*Tcore;
        dQ_nuCoredT = ( gsl_spline_eval( corethermalparams.Q_nuCore, Tcore+0.5*DeltaT, corethermalparams.Q_nuCore_acc ) - gsl_spline_eval( corethermalparams.Q_nuCore, Tcore-0.5*DeltaT, corethermalparams.Q_nuCore_acc ) )/DeltaT;

        // denom = 1. - a_ii*dm.Deltat*( dQ_nuCoredT/C_v + HeatFlux_pref/C_v ); //denominator of implicit calculation of T_core
        denom = 1. - a_ii*dm.Deltat*( dQ_nuCoredT/C_v ); //denominator of implicit calculation of T_core
        //add neutrino emissivity of entire core divided by heat capacity of entire core to HeatFlux to obtain update for Tcore
        // Tcore_new = ( Tcore_star + a_ii*dm.Deltat*( Q_nu/C_v - dQ_nuCoredT*Tcore/C_v - HeatFlux/C_v ) )/denom;
        Tcore_new = ( Tcore_star + a_ii*dm.Deltat*( Q_nu/C_v - dQ_nuCoredT*Tcore/C_v - HeatFlux/C_v ) )/denom;
        // std::cout << "T_core = " << Tcore_new*T_0 << "K, Q_nu - dQ_nuCoredT*Tcore = " << Q_nu - dQ_nuCoredT*Tcore << ", -HeatFlux = " << -HeatFlux <<  std::endl;
    }

    MPI_Bcast(&Tcore_new, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of QTcore to all processes
    MPI_Bcast(&HeatFlux_pref, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of HeatFlux_pref to all processes

    corethermalparams.HeatFlux = HeatFlux; //store last computed value of HeatFlux for use in boundary condition imposed in implicit timestepper for T
    corethermalparams.HeatFlux_pref = HeatFlux_pref; //store last computed value of HeatFlux_pref for use in boundary condition imposed in implicit timestepper for T
    corethermalparams.HeatFlux_H = HeatFlux_H; //store last computed value of HeatFlux_H for use in boundary condition imposed in implicit timestepper for T

    return Tcore_new;
}

/*
    Computes an intermediate auxiliary temperature appearing on the RHS of the implicit evolution
    Inputs: T_init: initial value (at start of time step) of redshifted temperature in reduced units
            T: IMEX-intermediate redshifted temperature in reduced units. This is the intermediate temperature for the IMEX substeps, not for the ADI substeps within each IMEX substep
            q, J, phE, cE: heat flux, redshifted current density and redshifted electric field and conjugate electric field times c (evaluated at location of J_phi) in reduced units
            trC: TransCoeffs object containing thermal conductivity kappa
            thC: TransCoeffs object containing specific heat capacity c_v and neutrino emissivity q_nu
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            tparams: TParams object containing information about temperature evolution
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process object containing information about the current process
            a_ijIm: off-diagonal entries in Butcher tableaux for this implicit RK step
    Output: Tstar: IMEX-intermediate auxiliary redshifted temperature in reduced units
*/
void Compute_Tstar_IMEX(ScalarField & Tstar, ScalarField & T_init, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH,
    ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, TParams & tparams, size_t N_GC, const Domain & dm, const Process & process, double a_ijIm)
{
    static double JH_pref = tparams.JouleHeating*B_0*B_0/(4.*pi*T_0*s_0); //prefactor to de-dimensionalize Joule heating term
    double JH; //Joule heating term without JH_pref.
    double Fr, Ftheta; //flux difference function in r and theta directions and cross term
    double r_iph, r_imh; //radial coordinate at upper and lower cell face
    double theta_jph, theta_jmh; //theta coordinate at right and left cell face
    double n_r_imh, n_r_iph; //sqrt(g_rr)=n_r at r = r_{i-1/2}, r_{i+1/2} respectively
    double Deltatheta = dm.Deltatheta;
    double gradrT_imhjm1, gradrT_iphjm1, gradrT_imhj, gradrT_iphj;
    double gradrT_ijm1, gradrT_ij, gradrT_ijmh;
    double gradthT_im1jmh, gradthT_im1jph, gradthT_ijmh, gradthT_ijph;
    double gradthT_im1j, gradthT_ij, gradthT_imhj;

    if(trC.conductivity_anisotropy == true){
        double kappa_delta_imh, kappa_perp_imh, kappa_H_imh; //kappa values computed at lower radial cell boundaries (index i-1/2) using harmonic mean
        double kappa_delta_jmh, kappa_perp_jmh, kappa_H_jmh; //kappa values computed at left theta cell boundaries (index j-1/2) using harmonic mean

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i];
            r_iph = dm.r[i] + 0.5*dm.Deltar[i];
            n_r_imh = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            n_r_iph = 0.5*( dm.n_r[i] + dm.n_r[i+1] );
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){
                theta_jmh = dm.theta[j] - 0.5*Deltatheta;
                theta_jph = dm.theta[j] + 0.5*Deltatheta;

                gradrT_imhjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i,j-1)/n_r_imh;
                gradrT_iphjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i+1,j-1)/n_r_iph;
                gradrT_imhj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i,j)/n_r_imh;
                gradrT_iphj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i+1,j)/n_r_iph;

                gradrT_ijm1 = MC( 0.5*(gradrT_imhjm1*r_imh+gradrT_iphjm1*r_iph), 2.*gradrT_imhjm1*r_imh, 2.*gradrT_iphjm1*r_iph )/dm.r[i];
                gradrT_ij = MC( 0.5*(gradrT_imhj*r_imh+gradrT_iphj*r_iph), 2.*gradrT_imhj*r_imh, 2.*gradrT_iphj*r_iph )/dm.r[i];
                gradrT_ijmh = MC( 0.5*(gradrT_ijm1*sin(dm.theta[j-1])+gradrT_ij*sin(dm.theta[j])), 2.*gradrT_ijm1*sin(dm.theta[j-1]), 2.*gradrT_ij*sin(dm.theta[j]) )/sin(theta_jmh);

                if( theta_jmh < 1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else if( theta_jmh > pi-1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else{
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1*sin(dm.theta[j-1])+gradrT_ij*sin(dm.theta[j])), 2.*gradrT_ijm1*sin(dm.theta[j-1]), 2.*gradrT_ij*sin(dm.theta[j]) )/sin(theta_jmh);
                }
                gradthT_im1jmh = gradth(T,Deltatheta*dm.r[i-1],i-1,j);
                gradthT_im1jph = gradth(T,Deltatheta*dm.r[i-1],i-1,j+1);
                gradthT_ijmh = gradth(T,Deltatheta*dm.r[i],i,j);
                gradthT_ijph = gradth(T,Deltatheta*dm.r[i],i,j+1);

                gradthT_im1j = MC( 0.5*(gradthT_im1jmh*sin(theta_jmh)+gradthT_im1jph*sin(theta_jph)), 2.*gradthT_im1jmh*sin(theta_jmh), 2.*gradthT_im1jph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_ij = MC( 0.5*(gradthT_ijmh*sin(theta_jmh)+gradthT_ijph*sin(theta_jph)), 2.*gradthT_ijmh*sin(theta_jmh), 2.*gradthT_ijph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_imhj = MC( 0.5*(gradthT_im1j*dm.r[i-1]+gradthT_ij*dm.r[i]), 2.*gradthT_im1j*dm.r[i-1], 2.*gradthT_ij*dm.r[i] )/r_imh;

                kappa_delta_jmh = 2./( 1./trC.kappa_delta[i][j-1] + 1./trC.kappa_delta[i][j] );
                kappa_perp_jmh = 2./( 1./trC.kappa_perp[i][j-1] + 1./trC.kappa_perp[i][j] );
                kappa_H_jmh = 2./( 1./trC.kappa_H[i][j-1] + 1./trC.kappa_H[i][j] );

                kappa_delta_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_delta[i-1][j] + dm.Deltar[i]/trC.kappa_delta[i][j] );
                kappa_perp_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_perp[i-1][j] + dm.Deltar[i]/trC.kappa_perp[i][j] );
                kappa_H_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_H[i-1][j] + dm.Deltar[i]/trC.kappa_H[i][j] );

                q[0][i][j] = -0.5*(dm.n_t[i-1]+dm.n_t[i])*( kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[0][i][j] + kappa_perp_imh )*gradrT_imhj
                            - 0.5*(dm.n_t[i-1]+dm.n_t[i])*( kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[1][i][j] - kappa_H_imh*mC.thBhat[2][i][j] )*
                                                            0.5/tan(Deltatheta/2.)*Deltatheta*gradthT_imhj;
                q[1][i][j] = -dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[1][i][j] + kappa_perp_jmh )*gradthT_ijmh
                            - dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[0][i][j] + kappa_H_jmh*mC.rBhat[2][i][j] )*gradrT_ijmh;
            }
        }
    }
    else{
        double kappa_imh, kappa_jmh;

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            n_r_imh = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){

                kappa_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa[i-1][j] + dm.Deltar[i]/trC.kappa[i][j] );
                q[0][i][j] = -0.5*(dm.n_t[i-1]+dm.n_t[i])*kappa_imh*gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i,j)/n_r_imh;
                if( j == 0 ){
                    q[1][i][j] = 0.; //sets heat flux in polar direction to zero at left edge of domain. If this edge is not at theta = 0, heat flux will be replaced by exchng2Scalar
                }
                else{
                    kappa_jmh = 2./( 1./trC.kappa[i][j-1] + 1./trC.kappa[i][j] );
                    q[1][i][j] = -dm.n_t[i]*kappa_jmh*gradth(T,dm.Deltatheta*dm.r[i],i,j);
                }
            }
        }
    }

    exchng2Vector(q,N_GC,process,dm.stridetype_Vec_q);

    for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
        r_iph = dm.r[i] + 0.5*dm.Deltar[i];
        r_imh = dm.r[i] - 0.5*dm.Deltar[i];
        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            theta_jph = dm.theta[j] + 0.5*dm.Deltatheta;
            theta_jmh = dm.theta[j] - 0.5*dm.Deltatheta;

            Fr = r_iph*r_iph*q[0][i+1][j] - r_imh*r_imh*q[0][i][j];
            Ftheta = sin(theta_jph)*q[1][i][j+1] - sin(theta_jmh)*q[1][i][j];

            JH = JouleHeating(i,j,J,phE,cE,trC,dm);

            Tstar[i][j] = T_init[i][j] + dm.Deltat*a_ijIm*( -Fr/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*2.*sin(dm.theta[j])*sin(dm.Deltatheta/2.))
                                                            + JH_pref*JH + 0*q_SH[i][j] + dm.n_t[i]*dm.n_t[i]*thC.q_nuCrust[i][j] )/thC.c_v[i][j];
        }
    }

    exchng2Scalar(Tstar,N_GC,process,dm.stridetype_Sca);

    return;
}

/*
    Computes the explicit ADI-split contribution to the stiff part of the temperature evolution
    This includes the divergence of the theta-derivative contribution to q_r and the entirety of the q_theta contribution, the entirety of the Jacobian terms in this divergence,
    the Joule heating term, and the neutrino emissivity term, all evaluated at the temperature T at the start of the RK substep
    Inputs: T, q, J, phE, cE: redshifted temperature, heat flux density, redshifted current density and redshifted electric field and conjugate electric field times c (evaluated at location of J_phi) in reduced units
            trC: TransCoeffs object containing thermal conductivity kappa
            thC: TransCoeffs object containing specific heat capacity c_v and neutrino emissivity q_nu
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process object containing information about the current process
    Output: KT: explicit ADI-split part of stiff update
*/
void Compute_KT_IMEX(ScalarField & KT, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH,
                    TransCoeffs & trC, ThermCoeffs & thC, MagCoeffs & mC, TParams & tparams, size_t N_GC, const Domain & dm, const Process & process)
{
    static double JH_pref = tparams.JouleHeating*B_0*B_0/(4.*pi*T_0*s_0); //prefactor to de-dimensionalize Joule heating term
    double JH; //Joule heating term without JH_pref.
    double Fr, Ftheta; //flux difference function in r and theta directions and cross term
    double r_iph, r_imh; //radial coordinate at upper and lower cell face
    double theta_jph, theta_jmh; //theta coordinate at right and left cell face
    double n_r_imh, n_r_iph; //sqrt(g_rr)=n_r at r = r_{i-1/2}, r_{i+1/2} respectively
    double Deltatheta = dm.Deltatheta;
    double gradrT_imhjm1, gradrT_iphjm1, gradrT_imhj, gradrT_iphj;
    double gradrT_ijm1, gradrT_ij, gradrT_ijmh;
    double gradthT_im1jmh, gradthT_im1jph, gradthT_ijmh, gradthT_ijph;
    double gradthT_im1j, gradthT_ij, gradthT_imhj;

    //Computes total q_theta and the part of q_r that includes finite differencing in the non-radial direction
    if(trC.conductivity_anisotropy == true){
        double kappa_delta_imh, kappa_H_imh; //kappa values computed at lower radial cell boundaries (index i-1/2) using harmonic mean
        double kappa_delta_jmh, kappa_perp_jmh, kappa_H_jmh; //kappa values computed at left theta cell boundaries (index j-1/2) using harmonic mean

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            r_imh = dm.r[i] - 0.5*dm.Deltar[i];
            r_iph = dm.r[i] + 0.5*dm.Deltar[i];
            n_r_imh = 0.5*( dm.n_r[i-1] + dm.n_r[i] );
            n_r_iph = 0.5*( dm.n_r[i] + dm.n_r[i+1] );
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){
                theta_jmh = dm.theta[j] - 0.5*Deltatheta;
                theta_jph = dm.theta[j] + 0.5*Deltatheta;

                gradrT_imhjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i,j-1)/n_r_imh;
                gradrT_iphjm1 = gradr(T,0.5*(dm.Deltar[i-1] + dm.Deltar[i]),i+1,j-1)/n_r_iph;
                gradrT_imhj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i,j)/n_r_imh;
                gradrT_iphj = gradr(T,0.5*(dm.Deltar[i] + dm.Deltar[i+1]),i+1,j)/n_r_iph;

                gradrT_ijm1 = MC( 0.5*(gradrT_imhjm1*r_imh+gradrT_iphjm1*r_iph), 2.*gradrT_imhjm1*r_imh, 2.*gradrT_iphjm1*r_iph )/dm.r[i];
                gradrT_ij = MC( 0.5*(gradrT_imhj*r_imh+gradrT_iphj*r_iph), 2.*gradrT_imhj*r_imh, 2.*gradrT_iphj*r_iph )/dm.r[i];

                if( theta_jmh < 1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else if( theta_jmh > pi-1e-10 ){
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1+gradrT_ij), 2.*gradrT_ijm1, 2.*gradrT_ij );
                }
                else{
                    gradrT_ijmh = MC( 0.5*(gradrT_ijm1*sin(dm.theta[j-1])+gradrT_ij*sin(dm.theta[j])), 2.*gradrT_ijm1*sin(dm.theta[j-1]), 2.*gradrT_ij*sin(dm.theta[j]) )/sin(theta_jmh);
                }
                gradthT_im1jmh = gradth(T,Deltatheta*dm.r[i-1],i-1,j);
                gradthT_im1jph = gradth(T,Deltatheta*dm.r[i-1],i-1,j+1);
                gradthT_ijmh = gradth(T,Deltatheta*dm.r[i],i,j);
                gradthT_ijph = gradth(T,Deltatheta*dm.r[i],i,j+1);

                gradthT_im1j = MC( 0.5*(gradthT_im1jmh*sin(theta_jmh)+gradthT_im1jph*sin(theta_jph)), 2.*gradthT_im1jmh*sin(theta_jmh), 2.*gradthT_im1jph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_ij = MC( 0.5*(gradthT_ijmh*sin(theta_jmh)+gradthT_ijph*sin(theta_jph)), 2.*gradthT_ijmh*sin(theta_jmh), 2.*gradthT_ijph*sin(theta_jph) )/sin(dm.theta[j]);
                gradthT_imhj = MC( 0.5*(gradthT_im1j*dm.r[i-1]+gradthT_ij*dm.r[i]), 2.*gradthT_im1j*dm.r[i-1], 2.*gradthT_ij*dm.r[i] )/r_imh;

                kappa_delta_jmh = 2./( 1./trC.kappa_delta[i][j-1] + 1./trC.kappa_delta[i][j] );
                kappa_perp_jmh = 2./( 1./trC.kappa_perp[i][j-1] + 1./trC.kappa_perp[i][j] );
                kappa_H_jmh = 2./( 1./trC.kappa_H[i][j-1] + 1./trC.kappa_H[i][j] );

                kappa_delta_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_delta[i-1][j] + dm.Deltar[i]/trC.kappa_delta[i][j] );
                kappa_H_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_H[i-1][j] + dm.Deltar[i]/trC.kappa_H[i][j] );

                q[0][i][j] = - 0.5*(dm.n_t[i-1]+dm.n_t[i])*( kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[1][i][j] - kappa_H_imh*mC.thBhat[2][i][j] )*
                                            0.5/tan(Deltatheta/2.)*Deltatheta*gradthT_imhj;
                q[1][i][j] = -dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[1][i][j] + kappa_perp_jmh )*gradthT_ijmh
                            - dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[0][i][j] + kappa_H_jmh*mC.rBhat[2][i][j] )*gradrT_ijmh;
            }
        }
    }
    else{
        double kappa_jmh;

        for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
            for(size_t j=N_GC; j<T.shape()[1]-N_GC+1; j++){
                q[0][i][j] = 0.;
                if( j == 0 ){
                    q[1][i][j] = 0.; //sets heat flux in polar direction to zero at left edge of domain. If this edge is not at theta = 0, heat flux will be replaced by exchng2Scalar
                }
                else{
                    kappa_jmh = 2./( 1./trC.kappa[i][j-1] + 1./trC.kappa[i][j] );
                    q[1][i][j] = -dm.n_t[i]*kappa_jmh*gradth(T,dm.Deltatheta*dm.r[i],i,j);
                }
            }
        }
    }

    exchng2Vector(q,N_GC,process,dm.stridetype_Vec_q);

    for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
        r_iph = dm.r[i] + 0.5*dm.Deltar[i];
        r_imh = dm.r[i] - 0.5*dm.Deltar[i];
        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            theta_jph = dm.theta[j] + 0.5*dm.Deltatheta;
            theta_jmh = dm.theta[j] - 0.5*dm.Deltatheta;

            Fr = r_iph*r_iph*q[0][i+1][j] - r_imh*r_imh*q[0][i][j];
            Ftheta = sin(theta_jph)*q[1][i][j+1] - sin(theta_jmh)*q[1][i][j];

            JH = JouleHeating(i,j,J,phE,cE,trC,dm);

            //Divide by overall factor of c_v[i][j]
            KT[i][j] = ( -Fr/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]) - Ftheta/(dm.r[i]*2.*sin(dm.theta[j])*sin(dm.Deltatheta/2.))
                        + JH_pref*JH + 0*q_SH[i][j] + dm.n_t[i]*dm.n_t[i]*thC.q_nuCrust[i][j] )/thC.c_v[i][j];
        }
    }

    exchng2Scalar(KT,N_GC,process,dm.stridetype_Sca);

    return;
}

/*
    Computes intermediate temperature update (ADI, r-direction) using Thomas algorithm.

    Inputs: T: redshifted temperature at previous substep in reduced units
            T_star: redshifted temperature of current IMEX substep, to be updated, in reduced units
            KT: non-stiff/ADI-split part of stiff update, evaluated at T
            trC: TransCoeffs object containing thermal conductivity kappa
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            tparams: TParams object containing information about boundary conditions
            corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            process: Process object containing information about the current process
            Tcore: redshifted core temperature in reduced units
            a_ii: diagonal entry in Butcher tableau for this implicit RK step.
    Output: T_int: intermediate updated redshifted temperature in reduced units. Used as an input for Compute_T_IMEX
*/
void Compute_Tint_IMEX(ScalarField & T_int, ScalarField & T, ScalarField & T_star, ScalarField & KT, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const TParams & tparams, const CoreThermalParams & corethermalparams, size_t N_GC,
                        const Domain & dm, const Process & process, double Tcore, double a_ii)
{
    double Deltat = dm.Deltat;
    size_t Nr = T.shape()[0]-2*N_GC-1; //final cell inside outer radius is i = T.shape()[0]-N_GC-1
    size_t Nth = T.shape()[1]-2*N_GC;
    double r_iph, r_imh; //radial coordinate at upper and lower cell face

    std::vector<double> c_ipr( Nr*Nth, 0. ); //vector of zeros for c_i' coefficients
    std::vector<double> d_ipr( Nr*Nth, 0. ); //vector of zeros for d_i' coefficients

    double ai, bi, ci, di;
    double b_0_term, d_0_term;
    double b_n_term, d_n_term; //coefficient of T^{n-1} and constant in expression for T^{n} from T_OuterBC_r: T^{n} = b_n_term*T^{n-1} + d_n_term
    double qr_prf_p1, qr_prf_m1, Ar_p1, Ar_m1, Ar_c;

    if(trC.conductivity_anisotropy == true){

        double kappa_delta_iph, kappa_delta_imh, kappa_perp_iph, kappa_perp_imh;
        T_InnerBC_r(b_0_term, d_0_term, Tcore, dm, corethermalparams, a_ii); //inner radial boundary condition; reuse same coefficients b_0, d_0 for each theta index j

        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
                r_iph = dm.r[i] + 0.5*dm.Deltar[i];
                r_imh = dm.r[i] - 0.5*dm.Deltar[i];

                kappa_delta_iph = (dm.Deltar[i]+dm.Deltar[i+1])/( dm.Deltar[i]/trC.kappa_delta[i][j] + dm.Deltar[i+1]/trC.kappa_delta[i+1][j] );
                kappa_delta_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_delta[i-1][j] + dm.Deltar[i]/trC.kappa_delta[i][j] );
                kappa_perp_iph = (dm.Deltar[i]+dm.Deltar[i+1])/( dm.Deltar[i]/trC.kappa_perp[i][j] + dm.Deltar[i+1]/trC.kappa_perp[i+1][j] );
                kappa_perp_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa_perp[i-1][j] + dm.Deltar[i]/trC.kappa_perp[i][j] );

                qr_prf_p1 = -( 0.5*(dm.n_t[i]+dm.n_t[i+1])*kappa_delta_iph*mC.thBhat[0][i+1][j]*mC.thBhat[0][i+1][j]
                                + 0.5*(dm.n_t[i]+dm.n_t[i+1])*kappa_perp_iph )/(0.5*(dm.Deltar[i] + dm.Deltar[i+1]))/(0.5*( dm.n_r[i] + dm.n_r[i+1] ));
                qr_prf_m1 = -( 0.5*(dm.n_t[i-1]+dm.n_t[i])*kappa_delta_imh*mC.thBhat[0][i][j]*mC.thBhat[0][i][j]
                                + 0.5*(dm.n_t[i-1]+dm.n_t[i])*kappa_perp_imh )/(0.5*(dm.Deltar[i-1] + dm.Deltar[i]))/(0.5*( dm.n_r[i-1] + dm.n_r[i] ));

                Ar_p1 = -r_iph*r_iph*qr_prf_p1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);
                Ar_m1 = -r_imh*r_imh*qr_prf_m1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);
                Ar_c = r_iph*r_iph*qr_prf_p1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]) + r_imh*r_imh*qr_prf_m1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);

                if (i == N_GC) {
                    //i = N_GC (inner radial boundary) entries
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] + Ar_m1*(2.*b_0_term-1.) );
                    ci = -Deltat*a_ii*Ar_p1;
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] + 2.*Ar_m1*d_0_term );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = ci/bi;
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = di/bi;
                }
                else if (i == T.shape()[0]-N_GC-2) {
                    //i = T.shape()[0]-N_GC-2 (outer radial boundary) entries
                    T_OuterBC_r(j, b_n_term, d_n_term, T, trC, mC, tparams, dm, N_GC);

                    ai = -Deltat*a_ii*Ar_m1;
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] + Ar_p1*b_n_term );
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] + Ar_p1*d_n_term );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = 0.;
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = ( di - ai*d_ipr[i-N_GC+(j-N_GC)*Nr-1] )/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                }
                else{
                    ai = -Deltat*a_ii*Ar_m1;
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] );
                    ci = -Deltat*a_ii*Ar_p1;
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = ci/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = ( di - ai*d_ipr[i-N_GC+(j-N_GC)*Nr-1] )/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                }
            }
        }

    }
    else{

        double kappa_iph, kappa_imh;
        T_InnerBC_r(b_0_term, d_0_term, Tcore, dm, corethermalparams, a_ii);

        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
                r_iph = dm.r[i] + 0.5*dm.Deltar[i];
                r_imh = dm.r[i] - 0.5*dm.Deltar[i];

                kappa_iph = (dm.Deltar[i]+dm.Deltar[i+1])/( dm.Deltar[i]/trC.kappa[i][j] + dm.Deltar[i+1]/trC.kappa[i+1][j] );
                kappa_imh = (dm.Deltar[i-1]+dm.Deltar[i])/( dm.Deltar[i-1]/trC.kappa[i-1][j] + dm.Deltar[i]/trC.kappa[i][j] );

                qr_prf_p1 = -0.5*(dm.n_t[i]+dm.n_t[i+1])*kappa_iph/(0.5*(dm.Deltar[i] + dm.Deltar[i+1]))/(0.5*( dm.n_r[i] + dm.n_r[i+1] ));
                qr_prf_m1 = -0.5*(dm.n_t[i-1]+dm.n_t[i])*kappa_imh/(0.5*(dm.Deltar[i-1] + dm.Deltar[i]))/(0.5*( dm.n_r[i-1] + dm.n_r[i] ));

                Ar_p1 = -r_iph*r_iph*qr_prf_p1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);
                Ar_m1 = -r_imh*r_imh*qr_prf_m1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);
                Ar_c = r_iph*r_iph*qr_prf_p1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]) + r_imh*r_imh*qr_prf_m1/(dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*thC.c_v[i][j]);

                if (i == N_GC) {
                    //i = N_GC (inner radial boundary) entries
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] + Ar_m1*(2.*b_0_term-1.) );
                    ci = -Deltat*a_ii*Ar_p1;
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] + 2.*Ar_m1*d_0_term );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = ci/bi;
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = di/bi;
                }
                else if (i == T.shape()[0]-N_GC-2) {
                    //i = T.shape()[0]-N_GC-2 (outer radial boundary) entries
                    T_OuterBC_r(j, b_n_term, d_n_term, T, trC, mC, tparams, dm, N_GC);

                    ai = -Deltat*a_ii*Ar_m1;
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] + Ar_p1*d_n_term );
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] + Ar_p1*b_n_term );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = 0.;
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = ( di - ai*d_ipr[i-N_GC+(j-N_GC)*Nr-1] )/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                }
                else {
                    ai = -Deltat*a_ii*Ar_m1;
                    bi = 1. - Deltat*a_ii*( Ar_c + dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j] );
                    ci = -Deltat*a_ii*Ar_p1;
                    di = T_star[i][j] + Deltat*a_ii*( KT[i][j] - dm.n_t[i]*thC.dq_nuCrustdT[i][j]/thC.c_v[i][j]*T[i][j] );

                    c_ipr[i-N_GC+(j-N_GC)*Nr] = ci/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                    d_ipr[i-N_GC+(j-N_GC)*Nr] = ( di - ai*d_ipr[i-N_GC+(j-N_GC)*Nr-1] )/( bi - ai*c_ipr[i-N_GC+(j-N_GC)*Nr-1] );
                }
            }
        }
    }

    //Backward sweep
    double T_intprev = 0.;
    for(size_t j=T.shape()[1]-N_GC-1; j>N_GC-1; j--){
        for(size_t i=T.shape()[0]-N_GC-2; i>N_GC-1; i--){
            T_int[i][j] = d_ipr[i-N_GC+(j-N_GC)*Nr] - c_ipr[i-N_GC+(j-N_GC)*Nr]*T_intprev;
            T_intprev = T_int[i][j];
        }
    }

    return;
}

/*
    Computes full temperature update (ADI, theta-direction) using Thomas algorithm.

    Inputs: T_int: ADI-intermediate redshifted temperature in reduced units
            T: redshifted initial temperature (at start of ADI step) in reduced units
            thC: ThermCoeffs object containing specific heat capacity c_v and neutrino emissivity q_nu
            trC: TransCoeffs object containing thermal conductivity kappa
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            process: Process object containing information about the current process
            dm: Domain object containing information about the simulation domain
            a_ii: diagonal entry in Butcher tableau for this implicit RK step.
    Output: T_f: ADI-updated redshifted temperatures in reduced units
*/
void Compute_T_IMEX(ScalarField & T_f, ScalarField & T_int, ScalarField & T, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const Process & process, size_t N_GC, const Domain & dm, double a_ii)
{
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    double Deltat = dm.Deltat;
    double Deltatheta = dm.Deltatheta;
    size_t Nr = T.shape()[0]-2*N_GC-1; //final cell inside outer radius is i = T.shape()[0]-N_GC-1
    size_t Nth = T.shape()[1]-2*N_GC;

    double theta_jph, theta_jmh; //theta coordinate at right and left cell face

    std::vector<double> c_ipr( Nr*Nth, 0. ); //vector of zeros for c_i' coefficients
    std::vector<double> d_ipr( Nr*Nth, 0. ); //vector of zeros for d_i' coefficients

    double T_np1_Outer; //estimate for boundary value of T in N_GCth ghost cell on either side
    double ai, bi, ci, di;
    double qth_prf_p1, qth_prf_m1; //q_th prefactors to (T^{j+1}-T^j) and (T^{j}-T^{j-1}) respectively
    double Ath_p1, Ath_m1, Ath_c; //div(q) prefactors to T^{j+1}, T^{j-1} and T^{j} respectively

    bool converged = false;
    size_t iteration = 0;
    double ConvTol = 1e-6; //fractional tolerance for convergence in Schwarz iteration

    std::vector<double> deltaTLeft_maxVec(process.num_procs);
    std::vector<double> deltaTRight_maxVec(process.num_procs);
    double deltaTLeft_max, deltaTRight_max;

    std::vector<double> T_np1_Outer_Left(T.shape()[0]); //edge (interior) temperature at left interior theta boundary
    std::vector<double> T_np1_Outer_Right(T.shape()[0]); //edge (interior) temperature at right interior theta boundary
    std::vector<double> deltaT_Left(T.shape()[0]); //fractional difference between temperature at left interior edge between iterations
    std::vector<double> deltaT_Right(T.shape()[0]); //fractional difference between temperature at right interior edge between iterations

    //Initialize initial guess for temperature in ghost cells at interior (between partial domain) boundaries. Place these guesses in the ghost cells of T_f
    if( nbrleft > 0 ) {
        T_InteriorBC_theta(N_GC-1, T_int, thC, trC, mC, dm, N_GC, a_ii, T_f);
    }
    if( nbrright > 0) {
        T_InteriorBC_theta(T.shape()[1]-N_GC, T_int, thC, trC, mC, dm, N_GC, a_ii, T_f);
    }

    while( converged == false ){
        //Forward sweep through matrix.
        // std::cout << "Forward sweep!" << std::endl;

        if(trC.conductivity_anisotropy == true){

            double kappa_delta_jph, kappa_delta_jmh, kappa_perp_jph, kappa_perp_jmh;

            for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
                for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++) {
                    theta_jph = dm.theta[j] + 0.5*dm.Deltatheta;
                    theta_jmh = dm.theta[j] - 0.5*dm.Deltatheta;

                    kappa_delta_jph = 2./( 1./trC.kappa_delta[i][j] + 1./trC.kappa_delta[i][j+1] );
                    kappa_delta_jmh = 2./( 1./trC.kappa_delta[i][j] + 1./trC.kappa_delta[i][j-1] );
                    kappa_perp_jph = 2./( 1./trC.kappa_perp[i][j] + 1./trC.kappa_perp[i][j+1] );
                    kappa_perp_jmh = 2./( 1./trC.kappa_perp[i][j] + 1./trC.kappa_perp[i][j-1] );

                    qth_prf_p1 = -dm.n_t[i]*( kappa_delta_jph*mC.rBhat[1][i][j+1]*mC.rBhat[1][i][j+1]
                                                    + kappa_perp_jph )/(Deltatheta*dm.r[i]);
                    qth_prf_m1 = -dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[1][i][j]
                                                    + kappa_perp_jmh )/(Deltatheta*dm.r[i]);

                    Ath_p1 = -sin(theta_jph)*qth_prf_p1/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);
                    Ath_m1 = -sin(theta_jmh)*qth_prf_m1/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);
                    Ath_c = (sin(theta_jph)*qth_prf_p1+sin(theta_jmh)*qth_prf_m1)/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);

                    if (j == N_GC) {
                        // j = N_GC (left theta boundary) entries
                        if( nbrleft < 0 ){
                            //Treatment if at theta=0 boundary
                            bi = 1. - Deltat*a_ii*Ath_c;
                            ci = - Deltat*a_ii*(Ath_p1+Ath_m1);
                            di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        else {
                            //Treatment if at a boundary between partial domains
                            bi = 1. - Deltat*a_ii*Ath_c;
                            ci = - Deltat*a_ii*Ath_p1;
                            di = T_int[i][j] + Deltat*a_ii*Ath_m1*T_f[i][j-1] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1]); //subtract off d_theta^2T term for initial T
                        }
                        c_ipr[j-N_GC+(i-N_GC)*Nth] = ci/bi;
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = di/bi;
                    }
                    else if(j == T.shape()[1]-N_GC-1){
                        // j = T.shape()[1]-N_GC-1 (right theta boundary) entries
                        if( nbrright < 0 ){
                            //Treatment if at theta=pi boundary
                            ai = - Deltat*a_ii*(Ath_m1+Ath_p1);
                            bi = 1. - Deltat*a_ii*Ath_c;
                            di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        else{
                            //Treatment if at a boundary between partial domains
                            ai = - Deltat*a_ii*Ath_m1;
                            bi = 1. - Deltat*a_ii*Ath_c;
                            di = T_int[i][j] + Deltat*a_ii*Ath_p1*T_f[i][j+1] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        c_ipr[j-N_GC+(i-N_GC)*Nth] = 0.;
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = ( di - ai*d_ipr[j-N_GC+(i-N_GC)*Nth-1] )/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] );
                    }
                    else{
                        ai = - Deltat*a_ii*Ath_m1;
                        bi = 1. - Deltat*a_ii*Ath_c;
                        ci = - Deltat*a_ii*Ath_p1;
                        di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T

                        c_ipr[j-N_GC+(i-N_GC)*Nth] = ci/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] );
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = ( di - ai*d_ipr[j-N_GC+(i-N_GC)*Nth-1] )/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] );
                    }
                }
            }
        }
        else{

            double kappa_jph, kappa_jmh;

            for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
                for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++) {
                    theta_jph = dm.theta[j] + 0.5*dm.Deltatheta;
                    theta_jmh = dm.theta[j] - 0.5*dm.Deltatheta;

                    kappa_jph = 2./( 1./trC.kappa[i][j] + 1./trC.kappa[i][j+1] );
                    kappa_jmh = 2./( 1./trC.kappa[i][j-1] + 1./trC.kappa[i][j] );

                    qth_prf_p1 = -dm.n_t[i]*kappa_jph/(Deltatheta*dm.r[i]);
                    qth_prf_m1 = -dm.n_t[i]*kappa_jmh/(Deltatheta*dm.r[i]);

                    Ath_p1 = -sin(theta_jph)*qth_prf_p1/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);
                    Ath_m1 = -sin(theta_jmh)*qth_prf_m1/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);
                    Ath_c = (sin(theta_jph)*qth_prf_p1+sin(theta_jmh)*qth_prf_m1)/(dm.r[i]*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*thC.c_v[i][j]);

                    if(j == N_GC){
                        // j = N_GC (left theta boundary) entries
                        if( nbrleft < 0 ) {
                            //Treatment if at theta=0 boundary
                            bi = 1. - Deltat*a_ii*Ath_c;
                            ci = - Deltat*a_ii*(Ath_p1+Ath_m1);
                            di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        else {
                            //Treatment if at a boundary between partial domains
                            bi = 1. - Deltat*a_ii*Ath_c;
                            ci = - Deltat*a_ii*Ath_p1;
                            di = T_int[i][j] + Deltat*a_ii*Ath_m1*T_f[i][j-1] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        c_ipr[j-N_GC+(i-N_GC)*Nth] = ci/bi;
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = di/bi;
                    }
                    else if(j == T.shape()[1]-N_GC-1){
                        // j = T.shape()[1]-N_GC-1 (right theta boundary) entries
                        if( nbrright < 0 ) {
                            //Treatment if at theta=pi boundary
                            ai = - Deltat*a_ii*(Ath_m1+Ath_p1);
                            bi = 1. - Deltat*a_ii*Ath_c;
                            di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] );
                        }
                        else {
                            //Treatment if at a boundary between partial domains
                            ai = - Deltat*a_ii*Ath_m1;
                            bi = 1. - Deltat*a_ii*Ath_c;
                            di = T_int[i][j] + Deltat*a_ii*Ath_p1*T_f[i][j+1] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] ); //subtract off d_theta^2T term for initial T
                        }
                        c_ipr[j-N_GC+(i-N_GC)*Nth] = 0.;
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = ( di - ai*d_ipr[j-N_GC+(i-N_GC)*Nth-1] )/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] );
                    }
                    else{
                        ai = - Deltat*a_ii*Ath_m1;
                        bi = 1. - Deltat*a_ii*Ath_c;
                        ci = - Deltat*a_ii*Ath_p1;
                        di = T_int[i][j] - Deltat*a_ii*( Ath_m1*T[i][j-1] + Ath_c*T[i][j] + Ath_p1*T[i][j+1] );

                        c_ipr[j-N_GC+(i-N_GC)*Nth] = ci/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] );
                        d_ipr[j-N_GC+(i-N_GC)*Nth] = ( di - ai*d_ipr[j-N_GC+(i-N_GC)*Nth-1] )/( bi - ai*c_ipr[j-N_GC+(i-N_GC)*Nth-1] ); //subtract off d_theta^2T term for initial T
                    }
                }
            }
        }

        // std::cout << "Backward sweep!" << std::endl;
        //Backward sweep
        double T_fprev = 0.;
        for(size_t i=T.shape()[0]-N_GC-2; i>N_GC-1; i--){
            for(size_t j=T.shape()[1]-N_GC-1; j>N_GC-1; j--){
                T_f[i][j] = d_ipr[j-N_GC+(i-N_GC)*Nth] - c_ipr[j-N_GC+(i-N_GC)*Nth]*T_fprev;
                T_fprev = T_f[i][j];
            }

            if(iteration > 0){
                deltaT_Left[i] = fabs((T_f[i][N_GC] - T_np1_Outer_Left[i])/T_f[i][N_GC]);
                deltaT_Right[i] = fabs((T_f[i][T.shape()[1]-N_GC-1] - T_np1_Outer_Right[i])/T_f[i][T.shape()[1]-N_GC-1]);
            }
            T_np1_Outer_Left[i] = T_f[i][N_GC];
            T_np1_Outer_Right[i] = T_f[i][T.shape()[1]-N_GC-1];
        }
        //Determine maximum relative error across all partial domains
        deltaTLeft_max = *std::max_element(deltaT_Left.begin(),deltaT_Left.end());
        deltaTRight_max = *std::max_element(deltaT_Right.begin(),deltaT_Right.end());
        MPI_Allgather(&deltaTLeft_max, 1, MPI_DOUBLE, deltaTLeft_maxVec.data(), 1, MPI_DOUBLE, process.comm1D);
        MPI_Allgather(&deltaTRight_max, 1, MPI_DOUBLE, deltaTRight_maxVec.data(), 1, MPI_DOUBLE, process.comm1D);
        deltaTLeft_max = *std::max_element(deltaTLeft_maxVec.begin(), deltaTLeft_maxVec.end());
        deltaTRight_max = *std::max_element(deltaTRight_maxVec.begin(), deltaTRight_maxVec.end());

        //If using one process, returns true and exits this loop. Otherwise, only exits if temperature at boundary has converged after Schwarz iteration.
        if(nbrleft < 0 && nbrright < 0){
            converged = true;
        }
        else if(iteration > 0 && deltaTLeft_max <= ConvTol && deltaTRight_max <= ConvTol){
            converged = true;
        }
        // std::cout << "iteration = " << iteration << std::endl;
        // std::cout << "deltaT_left = " << *std::max_element(deltaT_Left.begin(),deltaT_Left.end()) << ", deltaT_right = " << *std::max_element(deltaT_Right.begin(),deltaT_Right.end()) << std::endl;

        exchng2Scalar(T_f,N_GC,process,dm.stridetype_Sca); //fill ghost cells after each sweep
        iteration += 1; //increment iteration number by 1
        // std::cout << "iteration = " << iteration << ", converged = " << converged << ", process = " << process.world_rank << std::endl;

    }
    // std::cout << "iteration = " << iteration << std::endl;

    return;
}
