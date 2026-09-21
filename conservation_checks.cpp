#include <vector>
#include <iostream>

#include "common.h"
#include "field_evolution.h"
#include "conservation_checks.h"

inline double MC(double a, double b, double c) noexcept
{
    if ((a > 0.0) && (b > 0.0) && (c > 0.0))
        return std::min({a, b, c});

    if ((a < 0.0) && (b < 0.0) && (c < 0.0))
        return std::max({a, b, c});

    return 0.0;
}

/*
        Computes the divergence of the magnetic field in each cell and returns its total across current process
        Input: B: vector field containing magnetic field
               N_GC: number of ghost cells at each edge
               dm: Domain object containing information about the simulation domain
        Output: mean value of magnetic flux out of each cell. Should be zero to within machine precision.
*/
double divB_Calculator(VectorField & B, size_t N_GC, const Domain & dm)
{
    double Deltatheta = dm.Deltatheta;
    double sum = 0;
    for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){
            sum += 2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -(dm.r[i]-0.5*dm.Deltar[i])*(dm.r[i]-0.5*dm.Deltar[i])*B[0][i][j] + (dm.r[i]+0.5*dm.Deltar[i])*(dm.r[i]+0.5*dm.Deltar[i])*B[0][i+1][j] )
                    + dm.r[i]*dm.Deltar[i]*dm.n_r[i]*( -sin(dm.theta[j]-0.5*Deltatheta)*B[1][i][j] + sin(dm.theta[j]+0.5*Deltatheta)*B[1][i][j+1] );
        }
    }

    return sum;
}

/*
        Collects the mean div(B) from each process and averages them, then prints out the full domain-averaged div(B) value
        Input: B: vector field containing magnetic field
               N_GC: number of ghost cells at each edge
               dm: Domain object containing information about the simulation domain
               comm1D: communicator between processes
               world_rank: rank of current process
               Ntheta_locs: vector containing the extent of the domain in the decomposed direction
*/
void divB_Monitor(VectorField & B, size_t N_GC, const Domain & dm, MPI_Comm comm1D, int world_rank, std::vector<int> & Ntheta_locs)
{
    //Compute divergence of B for each cell, then take average and print out to track this. Should be 0 to within machine precision.
    static int num_procs = std::size(Ntheta_locs); //gives the number of processes
    std::vector<double> divBVec(num_procs);
    double divB = divB_Calculator(B, N_GC, dm);
    MPI_Allgather(&divB, 1, MPI_DOUBLE, divBVec.data(), 1, MPI_DOUBLE, comm1D);
    double divBSum;
    double num_cells = double( num_procs*(B.shape()[1]-2*N_GC)*(B.shape()[2]-2*N_GC) );
    if(world_rank == 0){
        divBSum = std::reduce(divBVec.begin(), divBVec.end())/double(num_procs);
        std::cout << "Average div(B) = " << divBSum/num_cells << std::endl;
    }
}

/*
    Calculates total magnetic field energy, Joule heating rate and outward Poynting flux integrated over decomposed domain
    Inputs: B, phE, cE, J: magnetic field, electric field and conjugate electric field times c (evaluated at location of J_phi), current density in reduced units
            tC: TransCoeffs instance containing either eta_O (isotropic conductivity) or eta_O_par and eta_O_perp (anisotropic conductivity)
            N_GC: number of ghost cells
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            dm: Domain object containing information about the simulation domain
    Output: U_B, JouleH, PoyntingF: magnetic field energy, Joule heating rate, Poynting flux integrated over decomposed domain
*/
void EnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & q_SH, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & JouleH, double & PoyntingF, double & SH)
{
    double Deltatheta = dm.Deltatheta;
    double U_BSum = 0.; //running sum of B^2/(8*pi) integrated over each cell's volume
    double JH_density, JouleSum = 0.; //J^2/sigma and running sum of J^2/sigma integrated over each cell's volume
    double PoyntingSum = 0.; //Poynting flux out of cell
    double SHSum = 0.;
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double EthBph, EphBth; //product of E_th*B_ph and E_ph*B_th used in computing radial Poynting flux
    double gr_m, gr_p; // n_t*n_r
    double Vij_r, Vij_th, Vij_ph;

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);//*dm.n_r[i]*dm.n_t[i];
            Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
            Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
            U_BSum += 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                        + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                    - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                        + Vij_ph*B[2][i][j]*B[2][i][j] );

            double JthEthH_ph_ij, JthEthH_ph_ip1j, JthEthH_ph_ijp1, JthEthH_ph_ip1jp1; //J_theta*E_theta_H evaluated at locations of J_phi
            double JrErH_ph_ij, JrErH_ph_ip1j, JrErH_ph_ijp1, JrErH_ph_ip1jp1; //J_r*E_r_H evaluated at locations of J_phi
            double JrEr_l, JrEr_r; //difference between J_r*E_r_H at location of J_r and the sum of half its value at each of the neighbouring J_phi locations
            double JthEth_b, JthEth_t; //difference between J_theta*E_theta_H at location of J_theta and the sum of half its value at each of the neighbouring J_phi locations

            if( dm.theta[j]-0.5*dm.Deltatheta > 0. ){
                JthEthH_ph_ij = 0.5*(sin(dm.theta[j-1])*J[1][i][j-1]+sin(dm.theta[j])*J[1][i][j])/sin(dm.theta[j]-0.5*dm.Deltatheta)*mC.phE_H[1][i][j];
                JthEthH_ph_ip1j = 0.5*(sin(dm.theta[j-1])*J[1][i+1][j-1]+sin(dm.theta[j])*J[1][i+1][j])/sin(dm.theta[j]-0.5*dm.Deltatheta)*mC.phE_H[1][i+1][j];
            }
            else{
                JthEthH_ph_ij = 0.;
                JthEthH_ph_ip1j = 0.;
            }

            if( dm.theta[j]+0.5*dm.Deltatheta + 1e-9 < pi ){
                JthEthH_ph_ijp1 = 0.5*(sin(dm.theta[j])*J[1][i][j]+sin(dm.theta[j+1])*J[1][i][j+1])/sin(dm.theta[j]+0.5*dm.Deltatheta)*mC.phE_H[1][i][j+1];
                JthEthH_ph_ip1jp1 = 0.5*(sin(dm.theta[j])*J[1][i+1][j]+sin(dm.theta[j+1])*J[1][i+1][j+1])/sin(dm.theta[j]+0.5*dm.Deltatheta)*mC.phE_H[1][i+1][j+1];
            }
            else{
                JthEthH_ph_ijp1 = 0.;
                JthEthH_ph_ip1jp1 = 0.;
            }

            JrErH_ph_ij = 0.5*(dm.r[i-1]*J[0][i-1][j]+dm.r[i]*J[0][i][j])/(dm.r[i]-0.5*dm.Deltar[i])*mC.phE_H[0][i][j];
            JrErH_ph_ijp1 = 0.5*(dm.r[i-1]*J[0][i-1][j+1]+dm.r[i]*J[0][i][j+1])/(dm.r[i]-0.5*dm.Deltar[i])*mC.phE_H[0][i][j+1];
            JrErH_ph_ip1j = 0.5*(dm.r[i]*J[0][i][j]+dm.r[i+1]*J[0][i+1][j])/(dm.r[i]+0.5*dm.Deltar[i])*mC.phE_H[0][i+1][j];
            JrErH_ph_ip1jp1 = 0.5*(dm.r[i]*J[0][i][j+1]+dm.r[i+1]*J[0][i+1][j+1])/(dm.r[i]+0.5*dm.Deltar[i])*mC.phE_H[0][i+1][j+1];

            JrEr_l = J[0][i][j]*mC.E_H_Pol[0][i][j] - 0.5*( JrErH_ph_ij + JrErH_ph_ip1j );
            JrEr_r = J[0][i][j+1]*mC.E_H_Pol[0][i][j+1] - 0.5*( JrErH_ph_ijp1 + JrErH_ph_ip1jp1 );
            JthEth_b = J[1][i][j]*mC.E_H_Pol[1][i][j] - 0.5*( JthEthH_ph_ij + JthEthH_ph_ijp1 );
            JthEth_t = J[1][i+1][j]*mC.E_H_Pol[1][i+1][j] - 0.5*( JthEthH_ph_ip1j + JthEthH_ph_ip1jp1 );

            JH_density = JouleHeating(i,j,J,phE,cE,tC,dm);// + ( JrEr_l + JrEr_r + JthEth_b + JthEth_t )/2.;

            JouleSum += -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*JH_density;
            SHSum += -q_SH[i][j]; //shock heating.

            if( i == N_GC ){
                //Poynting flux inward at inner radial surface
                // EthBph = 0.25*( phE[1][i-1][j] + phE[1][i][j] )*( B[2][i-1][j] + B[2][i][j] );
                EthBph = 0.25*( phE[1][i][j] + phE[1][i][j+1] )*( B[2][i-1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i][j]*( B[1][i-1][j] + B[1][i][j] ) + phE[2][i][j+1]*( B[1][i-1][j+1] + B[1][i][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( EthBph - EphBth ) );
            }
            if( i == B.shape()[1]-N_GC-2 ) {
                //Poynting flux outward at outer radial surface. Cell i = B.shape()[1]-N_GC-2 contains the upper boundary of the final interior cell, so it mostly depends on cells i = B.shape()[1]-N_GC-1
                // EthBph = 0.25*( phE[1][i][j]+ phE[1][i+1][j] )*( B[2][i+1][j] + B[2][i][j] );
                EthBph = 0.25*( phE[1][i+1][j]+ phE[1][i+1][j+1] )*( B[2][i+1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i+1][j]*( B[1][i][j] + B[1][i+1][j] ) + phE[2][i+1][j+1]*( B[1][i][j+1] + B[1][i+1][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( EthBph - EphBth );
            }
        }
    }

    //Set output variables equal to values computed by integrating over simulation domain
    U_B = U_BSum; //magnetic energy density integrated over volume; multiply by B_0^2*L_0^2 to get quantity in erg/cm
    JouleH = JouleSum; //Joule heating rate integrated over volume: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    PoyntingF = PoyntingSum; //Poynting flux through inner and outer radial surfaces : multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    SH = SHSum; //Shock heating rate integrated over volume: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time

    return;
}

/*
    Calculates total magnetic field energy, quasi-Joule heating rate (curl(B).E != J.E if B != H ), Joule heating rate and outward Poynting flux integrated over decomposed domain
    Inputs: B, phE, cE, J: magnetic field, electric field and conjugate electric field times c (evaluated at location of J_phi), current density
            tC: TransCoeffs instance containing either eta_O (isotropic conductivity) or eta_O_par and eta_O_perp (anisotropic conductivity)
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
    Output: U_B, QuasiJouleH, JouleH, PoyntingF: magnetic field energy, quasi-Joule heating rate, Joule heating rate, Poynting flux integrated over decomposed domain, all in reduced units
*/
void QuasiEnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & q_SH, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & QuasiJouleH, double & PoyntingF,
    double & JouleH, double & SH)
{
    double Deltatheta = dm.Deltatheta;
    double U_BSum = 0.; //B^2/(8*pi)*cell volume
    double QuasiJH_density, QuasiJouleSum = 0.; //J_B*E, running sum of J_B*E integrated over each cell's volume
    double JH_density, JouleSum = 0.; //J^2/sigma, running sum of J^2/sigma integrated over each cell's volume
    double PoyntingSum = 0.; //Poynting flux out of cell
    double SHSum = 0.;
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double EthBph, EphBth; //product of E_th*B_ph and E_ph*B_th used in computing radial Poynting flux
    double gr_m, gr_p; // n_t*n_r
    double Vij_r, Vij_th, Vij_ph;

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);//*dm.n_r[i]*dm.n_t[i];
            Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
            Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
            U_BSum += 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                        + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                    - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                        + Vij_ph*B[2][i][j]*B[2][i][j] );

            QuasiJH_density = QuasiJouleHeating(i,j,mC.J_B,phE,cE,tC,mC,dm);
            QuasiJouleSum += -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*QuasiJH_density;

            JH_density = JouleHeating(i,j,J,phE,cE,tC,dm);
            JouleSum += -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*JH_density;
            SHSum += -q_SH[i][j]; //shock heating.

            if( i == N_GC ){
                //Poynting flux inward at inner radial surface
                EthBph = 0.25*( phE[1][i-1][j] + phE[1][i][j] )*( B[2][i-1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i][j]*( B[1][i-1][j] + B[1][i][j] ) + phE[2][i][j+1]*( B[1][i-1][j+1] + B[1][i][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( EthBph - EphBth ) );
            }
            else if( i == B.shape()[1]-N_GC-2 ){
                //Poynting flux outward at outer radial surface. Cell i = B.shape()[1]-N_GC-2 contains the upper boundary of the final interior cell, so it mostly depends on cells i = B.shape()[1]-N_GC-1
                EthBph = 0.25*( phE[1][i][j]+ phE[1][i+1][j] )*( B[2][i+1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i+1][j]*( B[1][i][j] + B[1][i+1][j] ) + phE[2][i+1][j+1]*( B[1][i][j+1] + B[1][i+1][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( ( EthBph - EphBth ) );            }
        }
    }

    //Set output variables equal to values computed by integrating over simulation domain
    U_B = U_BSum; //magnetic energy density integrated over x-y area; multiply by B_0^2*L_0^2 to get quantity in erg/cm
    QuasiJouleH = QuasiJouleSum; //quasi-Joule heating rate integrated over x-y area: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    JouleH = JouleSum; //Joule heating rate integrated over x-y area: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    PoyntingF = PoyntingSum; //Poynting flux through x and y : multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    SH = SHSum; //Shock heating rate integrated over volume: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time

    return;
}

/*
    Initializes local energy in each cell. For testing code.
    Inputs: B: magnetic field  in reduced units
            UB_Init: Initial magnetic field energy in each cell
            JHPF_Cumulative: cumulative time integral of Joule heating and Poynting flux out of each cell
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            t: current time in reduced units
            world_rank: current process
*/
void EnergyConservationLocal_Initialize(VectorField & B, ScalarField & UB_Init, ScalarField & JHPF_Cumulative, size_t N_GC, const Domain & dm, double t, int world_rank)
{
    double Deltatheta = dm.Deltatheta;
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double Vij_r, Vij_th, Vij_ph;
    double gr_m, gr_p; // n_t*n_r

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);//*dm.n_r[i]*dm.n_t[i];
            Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
            Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
            UB_Init[i][j] = 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                        + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                    - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                        + Vij_ph*B[2][i][j]*B[2][i][j] );

            JHPF_Cumulative[i][j] = 0.;
        }
    }

    return;
}

/*
    Calculates local energy conservation and prints out to data file. For testing code.
    Inputs: B, phE, cE, J: magnetic field, electric field and conjugate electric field times c (evaluated at location of J_phi), current density in reduced units
            DeltaE: cumulative change in energy conservation in each cell. Ideally will be zero.
            UB_Init: Initial magnetic field energy in each cell
            JHPF_Cumulative: cumulative time integral of Joule heating and Poynting flux out of each cell
            tC: TransCoeffs instance containing either eta_O (isotropic conductivity) or eta_O_par and eta_O_perp (anisotropic conductivity)
            mC: MagCoeffs object containing properties of magnetic field and magnetization
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
            t: current time in reduced units
            world_rank: current process
*/
void EnergyConservationLocal(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & DeltaE, ScalarField & UB_Init, ScalarField & JHPF_Cumulative, TransCoeffs & tC, MagCoeffs & mC,
                            size_t N_GC, const Domain & dm, double t, const Process & process, int world_rank, size_t iter, size_t ECons_cadence, const BParams & bparams)
{
    double Deltatheta = dm.Deltatheta;
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double thJacobian1, thJacobian2; //Jacobian factors from radial part of (angular) volume integral. Should be ~ 1
    double PF, JH_density, JH, UB_new;
    double Vij_r, Vij_th, Vij_ph, Vij;
    double gr_m, gr_p; // n_t*n_r
    double Deltarimh, Deltariph, rimh, riph;
    double ar, arr, a0;
    double br, bth, brth, bthth, b0;
    double sigma_Bth_m, sigma_Bth_p, sigma_Br_m, sigma_Br_p;
    double Br_av, Bth_av;

    if( (iter == 1) || (iter%(ECons_cadence*100) == 0) ){

        // std::fstream Dataout;
        // std::fstream Dataout_Losses;
        // std::fstream Dataout_fractional;
        // std::string tstring = std::to_string(t*t_0/yr);
        // std::string RANK_NAME = std::to_string(world_rank);
        // Dataout.open("LocalECon/EConLocalCheck"+tstring+"_"+RANK_NAME+".dat",std::fstream::out);
        // Dataout_Losses.open("LocalECon/EConLossesLocalCheck"+tstring+"_"+RANK_NAME+".dat",std::fstream::out);
        // Dataout_fractional.open("LocalECon/EConFractionalErrorLocalCheck"+tstring+"_"+RANK_NAME+".dat",std::fstream::out);

        for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
            rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
            gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
            gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
            for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);//*dm.n_r[i]*dm.n_t[i];
                Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
                Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
                UB_new = 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                        + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                    - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                        + Vij_ph*B[2][i][j]*B[2][i][j] );

                JH_density = JouleHeating(i,j,J,phE,cE,tC,dm);
                JH = -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*JH_density;

                PF = -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( 0.25*(phE[1][i-1][j]+phE[1][i][j])*(B[2][i-1][j]+B[2][i][j]) - 0.25*( phE[2][i][j]*(B[1][i-1][j]+B[1][i][j]) + phE[2][i][j+1]*(B[1][i-1][j+1]+B[1][i][j+1]) ) ) )
                    -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( ( 0.25*(phE[1][i][j]+phE[1][i+1][j])*(B[2][i][j]+B[2][i+1][j]) - 0.25*( phE[2][i+1][j]*(B[1][i][j]+B[1][i+1][j]) + phE[2][i+1][j+1]*(B[1][i][j+1]+B[1][i+1][j+1]) ) ) )
                    -0.5*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*dm.n_t[i]*sin(dm.theta[j]-Deltatheta/2.)*( -( 0.25*( phE[2][i][j]*(B[0][i][j]+B[0][i][j-1]) + phE[2][i+1][j]*(B[0][i+1][j]+B[0][i+1][j-1]) ) - 0.25*(phE[0][i][j]+phE[0][i][j-1])*(B[2][i][j]+B[2][i][j-1]) ) )
                    -0.5*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*dm.n_t[i]*sin(dm.theta[j]+Deltatheta/2.)*( ( 0.25*( phE[2][i][j+1]*(B[0][i][j+1]+B[0][i][j]) + phE[2][i+1][j+1]*(B[0][i+1][j+1]+B[0][i+1][j]) ) - 0.25*(phE[0][i][j+1]+phE[0][i][j])*(B[2][i][j+1]+B[2][i][j]) ) );

                JHPF_Cumulative[i][j] = JHPF_Cumulative[i][j] + ( JH + PF )*dm.Deltat;
                DeltaE[i][j] = (UB_new - UB_Init[i][j]) - JHPF_Cumulative[i][j];

                // Dataout << std::setprecision(10) << DeltaE[i][j] << "   ";
                // Dataout_Losses << std::setprecision(10) << -JHPF_Cumulative[i][j] << "   ";
                // // Dataout_fractional << std::setprecision(10) << 1.-(UB_new-JHPF_Cumulative[i][j])/UB_Init[i][j] << "   ";
                // Dataout_fractional << std::setprecision(10) << 1.-(JHPF_Cumulative[i][j])/(UB_new-UB_Init[i][j]) << "   ";
            }
            // Dataout << std::endl;
            // Dataout << " " << std::endl;
            // Dataout_Losses << std::endl;
            // Dataout_Losses << " " << std::endl;
            // Dataout_fractional << std::endl;
            // Dataout_fractional << " " << std::endl;
        }
        // Dataout.close();
        // Dataout_Losses.close();
        // Dataout_fractional.close();
    }
    else{

        for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

                Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
                Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
                Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
                UB_new = 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                    + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                    + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                                + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                            - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                                + Vij_ph*B[2][i][j]*B[2][i][j] );

                JH_density = JouleHeating(i,j,J,phE,cE,tC,dm);
                JH = -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*JH_density;

                PF = -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( 0.25*(phE[1][i-1][j]+phE[1][i][j])*(B[2][i-1][j]+B[2][i][j]) - 0.25*( phE[2][i][j]*(B[1][i-1][j]+B[1][i][j]) + phE[2][i][j+1]*(B[1][i-1][j+1]+B[1][i][j+1]) ) ) )
                    -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( ( 0.25*(phE[1][i][j]+phE[1][i+1][j])*(B[2][i][j]+B[2][i+1][j]) - 0.25*( phE[2][i+1][j]*(B[1][i][j]+B[1][i+1][j]) + phE[2][i+1][j+1]*(B[1][i][j+1]+B[1][i+1][j+1]) ) ) )
                    -0.5*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*dm.n_t[i]*sin(dm.theta[j]-Deltatheta/2.)*( -( 0.25*( phE[2][i][j]*(B[0][i][j]+B[0][i][j-1]) + phE[2][i+1][j]*(B[0][i+1][j]+B[0][i+1][j-1]) ) - 0.25*(phE[0][i][j]+phE[0][i][j-1])*(B[2][i][j]+B[2][i][j-1]) ) )
                    -0.5*dm.r[i]*dm.Deltar[i]*dm.n_r[i]*dm.n_t[i]*sin(dm.theta[j]+Deltatheta/2.)*( ( 0.25*( phE[2][i][j+1]*(B[0][i][j+1]+B[0][i][j]) + phE[2][i+1][j+1]*(B[0][i+1][j+1]+B[0][i+1][j]) ) - 0.25*(phE[0][i][j+1]+phE[0][i][j])*(B[2][i][j+1]+B[2][i][j]) ) );

                JHPF_Cumulative[i][j] = JHPF_Cumulative[i][j] + ( JH + PF )*dm.Deltat;
            }
        }

    }

    return;
}

/*
    Computes average temperature over a decomposed domain. Simple average = sum( cell-centered T in each cell )/( number of cells )
    Inputs: T: temperature in reduced units
            N_GC: number of ghost cells
    Output: T_av: averaged temperature over decomposed domain in reduced units.
            T_min, T_max: minimum and maximum temperature across the partial domain
*/
void T_avCalc(ScalarField & T, size_t N_GC, double & T_av, double & T_min, double & T_max){

    double TSum = 0.;
    size_t cellcount = 0;
    T_min = T[N_GC][N_GC];
    T_max = T[N_GC][N_GC];

    for(size_t i=N_GC; i<T.shape()[0]-N_GC-1; i++){
        for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
            cellcount += 1;
            TSum += T[i][j];
            if(T[i][j] < T_min) T_min = T[i][j];
            if(T[i][j] > T_max) T_max = T[i][j];
        }
    }

    T_av = TSum/double(cellcount);

    return;
}

/*
    Computes crust volume-integrated neutrino emissivity in reduced units
    Inputs: q_nuCrust: neutrino emissivity in reduced units
            N_GC: number of ghost cells
            dm: Domain object containing information about the simulation domain
    Output: q_nuCrustIntegrated: volume-integrated neutrino emissivity of the crust in reduced units

*/
void Qnu_CrustIntegratedCalc(ScalarField & q_nuCrust, size_t N_GC, const Domain & dm, double & Q_nuCrustIntegrated){

    Q_nuCrustIntegrated  = 0.;
    double rJacobian;

    for(size_t i=N_GC; i<q_nuCrust.shape()[0]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        for(size_t j=N_GC; j<q_nuCrust.shape()[1]-N_GC; j++){
            Q_nuCrustIntegrated += q_nuCrust[i][j]*2.*pi*2.*sin(dm.theta[j])*sin(dm.Deltatheta/2.)*dm.n_r[i]*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian;
        }
    }

    return;
}
