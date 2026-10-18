#include <vector>
#include <iostream>

#include "common.h"
#include "field_evolution.h"
#include "conservation_checks.h"

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
void EnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & JouleH, double & PoyntingF)
{
    double Deltatheta = dm.Deltatheta;
    double U_BSum = 0.; //running sum of B^2/(8*pi) integrated over each cell's volume
    double JH_density, JouleSum = 0.; //J^2/sigma and running sum of J^2/sigma integrated over each cell's volume
    double PoyntingSum = 0.; //Poynting flux out of cell
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double EthBph, EphBth; //product of E_th*B_ph and E_ph*B_th used in computing radial Poynting flux
    double gr_m, gr_p; // n_t*n_r
    double Vij_r, Vij_th, Vij_ph;

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);
            Vij_th = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*dm.n_r[i]*dm.n_t[i];
            Vij_ph = 2.*pi*dm.r[i]*dm.r[i]*dm.Deltar[i]*rJacobian*2.*sin(0.5*Deltatheta)*sin(dm.theta[j])*dm.n_r[i]*dm.n_t[i];
            U_BSum += 1./(8.*pi)*( Vij_r*dm.Deltar[i]*( 0.5*dm.r[i]*dm.r[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )
                                + dm.r[i]*dm.Deltar[i]*( B[0][i+1][j]*B[0][i+1][j]*gr_p - B[0][i][j]*B[0][i][j]*gr_m )/6.
                                + dm.Deltar[i]*dm.Deltar[i]*( B[0][i][j]*B[0][i][j]*gr_m + B[0][i+1][j]*B[0][i+1][j]*gr_p )/24. )
                        + Vij_th*( sin(dm.theta[j])*sin(0.5*Deltatheta)*( B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1] )
                                    - ( B[1][i][j+1]*B[1][i][j+1] - B[1][i][j]*B[1][i][j] )*cos(dm.theta[j])*( cos(0.5*Deltatheta) - 2.*sin(0.5*Deltatheta)/Deltatheta ) )
                        + Vij_ph*B[2][i][j]*B[2][i][j] );

            JH_density = JouleHeating(i,j,J,phE,cE,tC,dm);

            JouleSum += -0.5*dm.Deltar[i]*dm.r[i]*dm.r[i]*dm.n_r[i]*rJacobian*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*JH_density;

            if( i == N_GC ){
                //Poynting flux inward at inner radial surface
                EthBph = 0.25*( phE[1][i][j] + phE[1][i][j+1] + cE[1][i][j] - cE[1][i][j+1] )*( B[2][i-1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i][j] + phE[2][i][j+1] )*( B[1][i-1][j] + B[1][i][j] + B[1][i-1][j+1] + B[1][i][j+1] );
                PoyntingSum += -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( EthBph - EphBth ) );
            }
            if( i == B.shape()[1]-N_GC-2 ) {
                //Poynting flux outward at outer radial surface. Cell i = B.shape()[1]-N_GC-2 contains the upper boundary of the final interior cell, so it mostly depends on cells i = B.shape()[1]-N_GC-1
                EthBph = 0.25*( phE[1][i+1][j] + phE[1][i+1][j+1] + cE[1][i+1][j] - cE[1][i+1][j+1] )*( B[2][i+1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i+1][j]*( B[1][i][j] + B[1][i+1][j] ) + phE[2][i+1][j+1]*( B[1][i][j+1] + B[1][i+1][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( EthBph - EphBth );
            }
        }
    }

    //Set output variables equal to values computed by integrating over simulation domain
    U_B = U_BSum; //magnetic energy density integrated over volume; multiply by B_0^2*L_0^2 to get quantity in erg/cm
    JouleH = JouleSum; //Joule heating rate integrated over volume: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    PoyntingF = PoyntingSum; //Poynting flux through inner and outer radial surfaces : multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time

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
void QuasiEnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & QuasiJouleH, double & PoyntingF,
    double & JouleH)
{
    double Deltatheta = dm.Deltatheta;
    double U_BSum = 0.; //B^2/(8*pi)*cell volume
    double QuasiJH_density, QuasiJouleSum = 0.; //J_B*E, running sum of J_B*E integrated over each cell's volume
    double JH_density, JouleSum = 0.; //J^2/sigma, running sum of J^2/sigma integrated over each cell's volume
    double PoyntingSum = 0.; //Poynting flux out of cell
    double rJacobian; //Jacobian factor from radial part of volume integral. Should be ~ 1
    double EthBph, EphBth; //product of E_th*B_ph and E_ph*B_th used in computing radial Poynting flux
    double gr_m, gr_p; // n_t*n_r
    double Vij_r, Vij_th, Vij_ph;

    for(size_t i=N_GC; i<B.shape()[1]-N_GC-1; i++){
        rJacobian = 1. + dm.Deltar[i]*dm.Deltar[i]/(dm.r[i]*dm.r[i])/12.;
        gr_m = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i-1]*dm.n_r[i-1]);
        gr_p = 0.5*(dm.n_t[i]*dm.n_r[i] + dm.n_t[i+1]*dm.n_r[i+1]);
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            Vij_r = 2.*pi*2.*sin(0.5*Deltatheta)*sin(dm.theta[j]);
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

            if( i == N_GC ){
                //Poynting flux inward at inner radial surface
                EthBph = 0.25*( phE[1][i-1][j] + phE[1][i][j] + cE[1][i][j] - cE[1][i][j+1] )*( B[2][i-1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i][j]*( B[1][i-1][j] + B[1][i][j] ) + phE[2][i][j+1]*( B[1][i-1][j+1] + B[1][i][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]-0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i-1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( -( EthBph - EphBth ) );
            }
            else if( i == B.shape()[1]-N_GC-2 ){
                //Poynting flux outward at outer radial surface. Cell i = B.shape()[1]-N_GC-2 contains the upper boundary of the final interior cell, so it mostly depends on cells i = B.shape()[1]-N_GC-1
                EthBph = 0.25*( phE[1][i][j]+ phE[1][i+1][j] + cE[1][i+1][j] - cE[1][i+1][j+1] )*( B[2][i+1][j] + B[2][i][j] );
                EphBth = 0.25*( phE[2][i+1][j]*( B[1][i][j] + B[1][i+1][j] ) + phE[2][i+1][j+1]*( B[1][i][j+1] + B[1][i+1][j+1] ) );
                PoyntingSum += -0.5*pow(dm.r[i]+0.5*dm.Deltar[i],2.)*0.5*(dm.n_t[i]+dm.n_t[i+1])*2.*sin(dm.theta[j])*sin(Deltatheta/2.)*( ( EthBph - EphBth ) );            }
        }
    }

    //Set output variables equal to values computed by integrating over simulation domain
    U_B = U_BSum; //magnetic energy density integrated over x-y area; multiply by B_0^2*L_0^2 to get quantity in erg/cm
    QuasiJouleH = QuasiJouleSum; //quasi-Joule heating rate integrated over x-y area: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    JouleH = JouleSum; //Joule heating rate integrated over x-y area: multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time
    PoyntingF = PoyntingSum; //Poynting flux through x and y : multiply by B_0^2*L_0^2 to get quantity in erg/cm/reduced time

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
