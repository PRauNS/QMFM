/*
        2D Finite volume electron MHD solver in axisymmetric spherical coordinates
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include <span>
#include <limits>
#include <list>
#include <filesystem>
namespace fs = std::filesystem;

#include "boost/multi_array.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <fftw3.h>
#include <shtns.h>

#include "common.h"
#include "initial_conditions.h"
#include "boundary_conditions.h"
#include "conservation_checks.h"
#include "microphysics.h"
#include "microphysics_Bdep.h"
#include "field_evolution.h"

const H5std_string FILE_EXT( ".h5" );

void RK_Step(VectorField & B, ScalarField & T, VectorField & phE, VectorField & cE, VectorField & J, VectorField & q, VectorField & B_np1, ScalarField & T_np1, ScalarField & q_SH,
    TransCoeffs & transCoeffs, ThermCoeffs & thermCoeffs, MagCoeffs & magCoeffs, const Domain & dm, const Process & process, const BParams & bparams, TParams & tparams, CoreThermalParams & corethermalparams, double t);

int main(int argc, char **argv) {
    int world_rank, num_procs;

    MPI_Init(&argc, &argv); /* Initialize the infrastructure necessary for communication */
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); /* Identify this process */
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); /* Find out how many total processes are active */

    std::string RANK_NAME = std::to_string(world_rank);
    std::string PROCS = std::to_string(num_procs);

    /*
        Load information from "SimSetup.in"
    */
    SimParams simparams; //object of SimParams class containing simulation parameters (see common.h for full list)
    BParams bparams; //object of BParams class containing initial magnetic field configuration and magnetic field boundary conditions (see common.h for full list)
    TParams tparams; //object of Tparams class containing initial temperature profile and thermal boundary conditions (see common.h for full list)
    ConductParams cparams; //object of Conductparams class containing conductivity calculation parameters (see common.h for full list)
    cparams.EorT = "Both";
    load_params(simparams, bparams, tparams, cparams, world_rank);

    const H5std_string FILE_NAME( simparams.OutputFile );

    /*
        Sets number of ghost cells needed at end of each simulation domain. Using
        monotonised central difference limiter (MC) method, so need two at each edge.
    */
    const size_t N_GC = 2;

    //Define problem domain in reduced units
    size_t Nr = simparams.Nr; //number of finite volume cells in r-direction, excluding ghost cells
    size_t Ntheta = simparams.Ntheta; //number of finite volume cells in theta-direction, excluding ghost cells
    double r_min = simparams.r_min, r_max = simparams.r_max;
    double theta_min = simparams.theta_min, theta_max = simparams.theta_max;

    if( (Nr%2 !=0 || Ntheta%2 !=0) && world_rank == 0 ){
        printf("Please choose number of finite volume cells in each direction to be even numbers\n");
    }
    if( (Ntheta % (size_t)num_procs != 0) && world_rank == 0 ){
        printf("Number of finite volume cells in theta direction cannot be evenly divided among processes\n");
    }

    double Lr = r_max-r_min; //extent of simulation domain in r-direction
    double Ltheta = theta_max-theta_min; //extent of simulation domain in theta-direction
    double Deltatheta = Ltheta/double(Ntheta); //angular size of each cell in theta-direction in reduced units

    double t = 0.; //current simulation time. Initialize to 0.
    const double t_max = simparams.t_max; //maximum simulation time

    /*
        Define finite volume cell centres and spacing between them.
    */
    std::vector<double> rFull; //Central points of cells over full domain
    std::vector<double> Deltar; //size of each cell in r-direction in reduced units
    Deltar.insert(Deltar.end(), { 0, 0 }); //ghost cell Deltar values. Will fill these with correct values later.
    if(simparams.varying_mesh == true){
        double r_scale = double(Nr)/2; //determines how far apart the cell centres are spaced in r. Larger values are more evenly spaced, while smaller values give more cells near the upper boundary.
        for(size_t i = 0; i<Nr-1; i++){
            if(i == 0){
                rFull.push_back( 0.5*r_min + 0.5*( r_min + Lr*( 1.-exp(-(double(i+1))/r_scale) )/( 1.-exp(-(double(Nr-1))/r_scale) ) ) ); //divides domain into Nr equally-spaced cells with r as their central points
                Deltar.push_back( 2.*( rFull.back() - r_min ) );
            }
            else{
                rFull.push_back( 0.5*( r_min + Lr*( 1.-exp(-(double(i+1))/r_scale) )/( 1.-exp(-(double(Nr-1))/r_scale) ) ) + 0.5*( rFull.back() + Deltar.back()/2. ) );
                Deltar.push_back( 2.*( rFull.back() - rFull[rFull.size()-2] ) - Deltar.back() );
            }
        }
        Deltar.push_back( Deltar.back() );
        rFull.push_back( rFull.back() + Deltar.back() );
    }
    else{
        for(size_t i = 0; i<Nr; i++){
            rFull.push_back( r_min + Lr*(double(i)+0.5)/double(Nr-1) ); //divides domain into Nr equally-spaced cells with r as their central points. Center of the final cell is outside the simulation domain!
            Deltar.push_back( Lr/double(Nr-1) ); //size of each cell in r-direction in reduced units. Subtract off one from denominator to account for non-periodic cells
        }
    }
    //Fill ghost cells in Deltar appropriately i.e., mirrored across the simulation boundaries so Deltar[0]=Deltar[3], Deltar[1]=Deltar[2] for N_GC=2
    for(size_t i=0; i<N_GC; i++){
        Deltar[N_GC-1-i] = Deltar[N_GC+i];
        Deltar.push_back( Deltar[Deltar.size()-3-2*i] );
    }
    std::vector<double> thetaFull; //Central points of cells over full domain
    thetaFull = linspace(theta_min+Deltatheta/2., theta_max-Deltatheta/2., Ntheta); //divides domain into Ntheta equally-spaced cells with theta as their central points

    /*
        Cartesian decomposition of the domain using MPI
    */
    int n_dims = 1; //number of spatial dimensions in the decomposition of the problem (not necessarily the spatial dimension of the problem)
    int dims[1] = {num_procs};
    int isperiodic[1] = {false};
    int reorder = true; //whether to allow MPI to reorder the partial domains
    MPI_Comm comm1D;
    MPI_Cart_create( MPI_COMM_WORLD, n_dims, dims, isperiodic, reorder, &comm1D ); //creates communicator comm1D

    int MyID;
    MPI_Comm_rank( comm1D, &MyID );
    int nbrleft, nbrright;
    MPI_Cart_shift( comm1D, 0, 1, &nbrleft, &nbrright );

    size_t MyS = 0, MyE = 0; //start and end indices of each partial domain, excluding ghost cells. Determined using MPE_Decomp1D
    MPE_Decomp1D(Ntheta, (size_t)num_procs, (size_t)MyID, MyS, MyE); //decomposes domain into num_procs partial domains in the theta-direction with Ntheta cells

    std::vector<int> starts(num_procs); //Vector of the starting indices (in decomposed domain direction) of each partial domain
    std::vector<int> Ntheta_locs(num_procs); //Vectors of the extent (in decomposed domain direction) of each partial domain
    //Send starts and Ntheta_locs to each process.
    MPI_Allgather(&MyS, 1, MPI_INT, starts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    int Ntheta_loc = int(MyE - MyS);
    MPI_Allgather(&Ntheta_loc, 1, MPI_INT, Ntheta_locs.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<double> r( rFull.begin(), rFull.end() ); //Local domain cell mid-points in r-direction. Unsplit.
    std::vector<double> theta( thetaFull.begin() + MyS, thetaFull.begin() + MyE ); //Local domain cell mid-points in theta-direction. Split based on the start and end indices computed by MPE_Decomp1D

    //Add ghost cell entries to r and theta vectors. For both r and theta, these are given by mirroring across the boundary. The mirror is
    for(size_t i=0; i<N_GC; i++){
        r.insert(r.begin(),r[2*i]);
        theta.insert(theta.begin(),theta[2*i]);
        r.push_back(r.back()+0.5*(Deltar[Nr+N_GC-1+i]+Deltar[Nr+N_GC+i]));
        theta.push_back(theta[theta.size()-1-2*i]);
    }
    r[Nr+N_GC-1] = r[Nr+N_GC-2]; //Since i=Nr+N_GC-1 cell is outside of simulation domain, set its r coordinate equal to that of i=Nr+N_GC-2 cell

    std::vector<double> rDomain( rFull.begin(), rFull.end() ); //Local domain cell mid-points in r-direction. Unsplit. Excludes ghost cells
    std::vector<double> thetaDomain( thetaFull.begin() + MyS, thetaFull.begin() + MyE ); //Local domain cell mid-points in theta-direction. Split based on the start and end indices computed by MPE_Decomp1D. Excludes ghost cells

    /*
        Initialize magnetic field and temperature. Define fields for updated magnetic field B_np1, electric field E, current density J
    */

    const double rho_cutoff = simparams.rho_cutoff; //cutoff (energy) density in g/cm^3 (lowest density to include in simulation domain).
    EOSInterpolation EOS_Interps; //object of EOSInterpolation class to contain EOS interpolating functions as function of "radial" coordinate r
    double g_s14, r_s, n_t_s, n_b_nd; //gravitational acceleration at the surface of the star in cm/s^2, outer radius of the star in reduced units, lapse function n_t = sqrt(-g_{tt}) = exp(nu/2) at outer radius of star,
    //neutron drip density n_b_nd in fm^-3
    load_EOS(simparams.CrustEOS, rho_cutoff, g_s14, r_s, n_t_s, n_b_nd, EOS_Interps); //loads EOS/stellar structure and generates interpolators for crust microphysics in EOS_Interps
    tparams.g_s14 = g_s14; //sets gravitational acceleration at surface of star based on EOS data
    tparams.r_s = r_s; //sets outer (surface) radius of star based on EOS data
    tparams.n_t_s = n_t_s; //sets lapse function at outer radius of star based on EOS data
    cparams.n_b_nd = n_b_nd; //sets neutron drip density in fm^-3 based on EOS data

    //Determine metric coefficients within crust from background stellar structure model
    std::vector<double> n_t, n_r; //n_t = sqrt(-g_{tt}), n_r = sqrt(g_{rr})
    GR_Factor_Initialize(simparams.GR, EOS_Interps, r, Nr, N_GC, n_t, n_r); //initialize n_t and n_r values. If params.GR == false, these are set to 1.

    double n_router;
    if(simparams.GR == true){
        n_router = sqrt( gsl_spline_eval(EOS_Interps.grr_spline, r_max, EOS_Interps.grr_acc) );
    }
    else n_router = 1;

    //Define struct process, containing data about the current individual process, respectively.
    struct Process process = {simparams.Timestep_method, world_rank, num_procs, nbrleft, nbrright, MyS, MyE, comm1D};
    exchng2Array(theta,N_GC,process); //Exchange theta values across partial domains

    int count_Vec = static_cast<int>(3*(Nr+2*N_GC)); //number of vector field components (3-component) times extent of vector field in r-dimension (including ghost cells).
    int count_Vec_q = static_cast<int>(2*(Nr+2*N_GC)); //number of vector field components (2-component) times extent of vector field in r-dimension (including ghost cells).
    int count_Sca = static_cast<int>(Nr+2*N_GC); //Extent of vector field in r-dimension (including ghost cells).
    int blocklength = static_cast<int>(N_GC);
    int stride = static_cast<int>(MyE-MyS+2*N_GC); //Extent of vector field in theta-dimension. includes ghost cells.
    MPI_Datatype stridetype_Vec;
    MPI_Type_vector( count_Vec, blocklength, stride, MPI_DOUBLE, &stridetype_Vec);
    MPI_Type_commit( &stridetype_Vec );
    MPI_Datatype stridetype_Vec_q;
    MPI_Type_vector( count_Vec_q, blocklength, stride, MPI_DOUBLE, &stridetype_Vec_q);
    MPI_Type_commit( &stridetype_Vec_q );
    MPI_Datatype stridetype_Sca;
    MPI_Type_vector( count_Sca, blocklength, stride, MPI_DOUBLE, &stridetype_Sca);
    MPI_Type_commit( &stridetype_Sca );

    //Define structs domain, containing data about the overall simulation domain. Set Deltat = 0 here since have not used CFL condition to compute it yet. Deltatheta_dfactor = 1. - 0.5*Deltatheta/tan(0.5*Deltatheta) is a common factor
    // used repeatedly in determining the theta-direction cell centroids
    struct Domain domain = {Nr, Ntheta, N_GC, Lr, Ltheta, Deltatheta, 1. - 0.5*Deltatheta/tan(0.5*Deltatheta), Deltar, r, theta, n_t, n_r, n_router,
                            r_min, r_max, EOS_Interps.RStar, 0., Ntheta_locs, starts, 0., simparams.InitiallyEquatoriallySymmetric, stridetype_Vec, stridetype_Vec_q, stridetype_Sca};

    //Compute initial magnetic field and create VectorField objects to hold electric field and updated magnetic field
    VectorField B(boost::extents[3][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of B across partial domain
    InitializeB(r, theta, bparams, domain, N_GC, Deltar, Deltatheta, B); //Generate initial values of B components by cell-face averaging over initial functional form
    B_BoundaryConditions(B, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, domain); //Impose boundary conditions on magnetic field
    exchng2Vector(B, N_GC, process, stridetype_Vec);
    VectorField B_np1(boost::extents[3][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of B in reduced units at next time step
    VectorField phE(boost::extents[3][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of redshifted E times c in reduced units (E*c*exp(nu/2)/(B_0*L_0/t_0)). Evaluated at cell corners (location of J_phi)
    VectorField cE(boost::extents[3][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of redshifted conjugate E times c in reduced units (E*c*exp(nu/2)/(B_0*L_0/t_0)). Evaluated at cell corners (location of J_phi)
    VectorField J(boost::extents[3][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-edge average values of redshifted J in reduced units (J*exp(nu/2)*L_0/B_0*4*pi/c)

    /*
        Define physical parameters and simulation domain
    */

    //Declare objects containing ScalarFields/RadialScalarFields for transport, thermodynamic and magnetization coefficients
    //These objects are referenced throughout the code.
    TransCoeffs transCoeffs(Nr+2*N_GC, MyE-MyS+2*N_GC, tparams.conductivity_anisotropy, simparams.C_hyp);
    ThermCoeffs thermCoeffs(Nr+2*N_GC, MyE-MyS+2*N_GC, simparams.IMEX);
    MagCoeffs magCoeffs(Nr+2*N_GC, MyE-MyS+2*N_GC, bparams.quantization, bparams.PolarAxisQuantization);

    ScalarField n_e(boost::extents[Nr+2*N_GC][MyE-MyS+2*N_GC]); //electron number density in fm^{-3}
    RadialScalarField A(boost::extents[Nr+2*N_GC]); //mass number
    RadialScalarField Z(boost::extents[Nr+2*N_GC]); //atomic number
    RadialScalarField n_i(boost::extents[Nr+2*N_GC]); //number density of ions in fm^{-3}
    RadialScalarField n_b(boost::extents[Nr+2*N_GC]); //total baryon number density in fm^{-3}
    RadialScalarField n_nf(boost::extents[Nr+2*N_GC]); //free neutron number density in fm^{-3}
    RadialScalarField mu_nf(boost::extents[Nr+2*N_GC]); //free neutron chemical potential in fm^{-3}

    for(size_t i=0; i<Nr+2*N_GC; i++){
        if( i < N_GC ){
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, EOS_Interps.R_cc, EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, EOS_Interps.R_cc, EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, EOS_Interps.R_cc, EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, EOS_Interps.R_cc, EOS_Interps.n_i_acc);
            n_nf[i] = gsl_spline_eval(EOS_Interps.n_nf_spline, EOS_Interps.R_cc, EOS_Interps.n_nf_acc);
            mu_nf[i] = gsl_spline_eval(EOS_Interps.mu_nf_spline, EOS_Interps.R_cc, EOS_Interps.mu_nf_acc);
        }
        else if( i > Nr+N_GC-2 ){
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, EOS_Interps.R_rhocutoff, EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, EOS_Interps.R_rhocutoff, EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_i_acc);
            n_nf[i] = gsl_spline_eval(EOS_Interps.n_nf_spline, EOS_Interps.R_rhocutoff, EOS_Interps.n_nf_acc);
            mu_nf[i] = gsl_spline_eval(EOS_Interps.mu_nf_spline, EOS_Interps.R_rhocutoff, EOS_Interps.mu_nf_acc);
        }
        else{
            A[i] = gsl_spline_eval(EOS_Interps.A_spline, r[i], EOS_Interps.A_acc);
            Z[i] = gsl_spline_eval(EOS_Interps.Z_spline, r[i], EOS_Interps.Z_acc);
            n_b[i] = gsl_spline_eval(EOS_Interps.n_b_spline, r[i], EOS_Interps.n_b_acc);
            n_i[i] = gsl_spline_eval(EOS_Interps.n_i_spline, r[i], EOS_Interps.n_i_acc);
            n_nf[i] = gsl_spline_eval(EOS_Interps.n_nf_spline, r[i], EOS_Interps.n_nf_acc);
            mu_nf[i] = gsl_spline_eval(EOS_Interps.mu_nf_spline, r[i], EOS_Interps.mu_nf_acc);
        }
        for(size_t j=0; j<MyE-MyS+2*N_GC; j++){
            if( i < N_GC ){
                //Mirrored values across inner boundary. Boundary is at mid-point between cells i=N_GC-1 and N_GC
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, r[2*N_GC-1-i], EOS_Interps.n_e_acc);
            }
            else if( i > Nr+N_GC-2 ){
                //Mirrored values across outer boundary. Boundary is at mid-point between cells i=Nr+N_GC-1 and Nr+N_GC-2
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, r[2*Nr-i+1], EOS_Interps.n_e_acc);
            }
            else{
                n_e[i][j] = gsl_spline_eval(EOS_Interps.n_e_spline, r[i], EOS_Interps.n_e_acc);
            }
        }
    }

    //    double sigma_const = 1e23, ne_const = 2e-6;
    //Set cell-averaged value of Hall diffusivity. Also includes ghost cells, which should be filled identically to nearest cell.
    for(size_t i=0; i<Nr+2*N_GC; i++){
        for(size_t j=0; j<MyE-MyS+2*N_GC; j++){
            //            transCoeffs.eta_O[i][j] = ( 6.*pow( (r[i]-r_min)/(r_max-r_min),3.5 ) + 3.*pow( (r_max-r[i])/(r_max-r_min),4. ) )*(1e5*1e5/(1e6*yr))*t_0/(L_0*L_0); //Ohmic diffusivity in reduced units from Vigano et al 2021
            // transCoeffs.eta_O[i][j] = 0.;//c*c*t_0/(4.*pi*sigma_const*L_0*L_0); //Ohmic diffusivity in reduced units
            //            transCoeffs.eta_H[i][j] = c*B_0*t_0/(4.*pi*unit_e*ne_const*1e39*L_0*L_0); //Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
            transCoeffs.eta_H[i][j] = c*B_0*t_0/(4.*pi*unit_e*n_e[i][j]*1e39*L_0*L_0); //Hall diffusivity in reduced units. 1e39 factor converts n_e from fm^-3 to cm^-3
            transCoeffs.psi_H[i][j] = transCoeffs.eta_H[i][j]/( domain.r[i]*domain.r[i]*domain.n_t[i] ); //Hall stream function for calculation of advection velocities in toroidal B evolution
            magCoeffs.chi[i][j] = 0.; //initializes differential magnetic susceptibility values to zero- needed for non-quantizing case.
            magCoeffs.vtheta_phi_lb_Lag[i][j] = 0.;
            magCoeffs.vtheta_phi_lt_Lag[i][j] = 0.;
            magCoeffs.vtheta_phi_rb_Lag[i][j] = 0.;
            magCoeffs.vtheta_phi_rt_Lag[i][j] = 0.;
            if(magCoeffs.quantization == true){
                magCoeffs.C_m[i][j] = 0.; //initializes magnetocaloric coefficient values to zero- needed for fixed temperature case.
            }
        }
    }
    exchng2Scalar(transCoeffs.eta_H,N_GC,process,stridetype_Sca);

    transCoeffs.max_eta_local = 0.; //maximum value of magnetic diffusivity (both Hall and Ohmic parts) in reduced units, including correction for including correction factor for nonzero magnetic susceptibility (which is 1 if Landau quantization is disabled)
    if(transCoeffs.conductivity_anisotropy == true) Bmag_and_BhatCalc(B, magCoeffs, transCoeffs.eta_H, transCoeffs.eta_O_perp, transCoeffs.max_eta_local, domain, N_GC, process); //calculates magnitude of magnetic field and magnetic field unit vector at cell centers
    else Bmag_and_BhatCalc(B, magCoeffs, transCoeffs.eta_H, transCoeffs.eta_O, transCoeffs.max_eta_local, domain, N_GC, process);

    ScalarField T(boost::extents[Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of redshifted temperature across partial domain
    ScalarField T_np1(boost::extents[Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of temperature in reduced units at next time step
    VectorField q(boost::extents[2][Nr+2*N_GC][MyE-MyS+2*N_GC]); //Cell-face average values of heat flux density times lapse function squared across partial domain in reduced units q*exp(nu)*t_0/(s_0*T_0*L_0) . Computed once transport coefficients are computed

    InitializeT(r, theta, tparams, domain, N_GC, Deltar, Deltatheta, T); //Generate initial values of T by cell-face averaging over initial functional form
    tparams.Tcore = tparams.T_init; //sets Tcore equal to T_init initially.

    //Calculate electron chemical potential given an electron number density profile. Assume mu_e is fixed afterward.
    Compute_mu_e(n_e, magCoeffs.mu_e);

    //Compute magnetization coefficients and/or specific heat capacity
    OmegaB_deriv_Interpolation OmB_Interpolators; //instantiate OmegaB_deriv_Interpolation object OmB_Interpolators.
    Compute_OmegaB_Interpolators(OmB_Interpolators); //generate interpolating functions used when computing magnetization-related coefficients. Always need this because use Landau-quantized electron specific heat capacity.
    if( magCoeffs.quantization == true ){
        if( tparams.fixed_T == false ){
            Omega_xyVaryingT(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, magCoeffs.M, magCoeffs.chi, magCoeffs.M_mu, thermCoeffs.c_v, magCoeffs.C_m);
            c_vfuncCalcB(T, n_e, A, Z, n_i, n_b, n_nf, mu_nf, tparams, N_GC, domain.n_t, thermCoeffs.c_v);

            //Exchange transport coefficients on periodic dimension.
            exchng2Scalar(thermCoeffs.c_v,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.M,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.chi,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.M_mu,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.C_m,N_GC,process,stridetype_Sca);

        }
        else{
            Omega_xyFixedT(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, magCoeffs.M, magCoeffs.chi, magCoeffs.M_mu);
            Omega_xyVaryingTC_veOnly(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, thermCoeffs.c_v);
            c_vfuncCalc(T, n_e, A, Z, n_i, n_b, n_nf, mu_nf, tparams, N_GC, domain.n_t, thermCoeffs.c_v); //Don't need to compute c_v if T is held fixed

            //Exchange transport coefficients on periodic dimension.
            exchng2Scalar(magCoeffs.M,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.chi,N_GC,process,stridetype_Sca);
            exchng2Scalar(magCoeffs.M_mu,N_GC,process,stridetype_Sca);
        }
    }
    else{
        Omega_xyVaryingTC_veOnly(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, thermCoeffs.c_v);
        c_vfuncCalc(T, n_e, A, Z, n_i, n_b, n_nf, mu_nf, tparams, N_GC, domain.n_t, thermCoeffs.c_v);

        //Exchange transport coefficients on periodic dimension.
        exchng2Scalar(thermCoeffs.c_v,N_GC,process,stridetype_Sca);
    }

    double max_eta_T_local = 0.; //maximum value of thermal diffusivity in reduced units
    //Compute transport coefficients, depending on whether anisotropic conductivity is turned on or off
    if( tparams.conductivity_anisotropy == true ){
        ConductivityCalcB(T, magCoeffs.Bmag, n_e, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, thermCoeffs.c_v, transCoeffs.eta_O_delta, transCoeffs.eta_O_perp, transCoeffs.kappa_delta, transCoeffs.kappa_perp, transCoeffs.kappa_H, max_eta_T_local, process, domain);

        //Exchange transport coefficients on periodic dimension.
        exchng2Scalar(transCoeffs.eta_O_delta,N_GC,process,stridetype_Sca);
        exchng2Scalar(transCoeffs.eta_O_perp,N_GC,process,stridetype_Sca);
        exchng2Scalar(transCoeffs.kappa_delta,N_GC,process,stridetype_Sca);
        exchng2Scalar(transCoeffs.kappa_perp,N_GC,process,stridetype_Sca);
        exchng2Scalar(transCoeffs.kappa_H,N_GC,process,stridetype_Sca);
    }
    else{
        ConductivityCalc(T, n_e, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, thermCoeffs.c_v, transCoeffs.eta_O, transCoeffs.kappa, max_eta_T_local);

        //Exchange transport coefficients on periodic dimension.
        exchng2Scalar(transCoeffs.eta_O,N_GC,process,stridetype_Sca);
        exchng2Scalar(transCoeffs.kappa,N_GC,process,stridetype_Sca);
    }

    //Exchange transport coefficients on periodic dimension.
    q_nuCrustCalc(T, n_e, magCoeffs.Bmag, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, simparams.IMEX, thermCoeffs.q_nuCrust, thermCoeffs.dq_nuCrustdT); //Compute neutrino emissivity (this is the same for all simulations)
    exchng2Scalar(thermCoeffs.q_nuCrust, N_GC, process, stridetype_Sca);
    if(simparams.IMEX == true){
        exchng2Scalar(thermCoeffs.dq_nuCrustdT, N_GC, process, stridetype_Sca);
    }

    double B_pol; //magnetic field magnitude at poles in reduced units. Determined from process including theta=0 cell and broadcast to the others.
    if(process.nbrleft < 0) B_pol = B[0][B.shape()[1]-N_GC-1][N_GC]; //determines B_pol from process with nbrleft = -2 (i.e., process containing north pole theta=0)
    MPI_Bcast(&B_pol, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of B_pol to all processes
    tparams.B_pol = B_pol; //assign updated value of B_pol to tparams object, to be used in implementing anisotropic (B-dependent) thermal boundary condition

    //Compute core thermal properties interpolators to be used in temperature inner boundary condition
    size_t N_tCore = 400; //number of core temperatures to use in constructing core thermal property interpolating functions
    double T_core_min = 7e6, T_core_max = 3e10; //maximum and minimum values of uniform redshifted core temperature in K to use in constructing core thermal property interpolating functions
    std::vector<double> T_coreVec; //vector of uniform redshifted core temperatures in reduced units to use in constructing core thermal property interpolating functions
    double C_vCore, Q_nu;
    std::vector<double> C_vCoreVec, Q_nuCoreVec; //vector of core-integrated specific heat capacity and neutrino emissivity to use in constructing core thermal property interpolating functions
    double kappa_eCC, kappa_nCC;

    for(size_t j=0; j<N_tCore; j++){
        T_coreVec.push_back( pow( 10.,log10(T_core_min)+double(j)/double(N_tCore-1)*(log10(T_core_max)-log10(T_core_min)) )/T_0 );
        ThermPhysCore(simparams.CoreEOS, simparams.GR, tparams, T_coreVec.back(), C_vCore, Q_nu, kappa_eCC, kappa_nCC);
        C_vCoreVec.push_back(C_vCore);
        Q_nuCoreVec.push_back(Q_nu);
    }
    CoreThermalParams corethermalparams; //define CoreThermalParams object "corethermalparams"
    Compute_CoreThermalParams_Interpolators(T_coreVec,C_vCoreVec,Q_nuCoreVec,corethermalparams); //create interpolation functions for C_vCore and Q_nu and store them in corethermalparams
    transCoeffs.kappa_eCore = kappa_eCC;
    transCoeffs.kappa_nCore = kappa_nCC;
    corethermalparams.Q_nuCoreVal = gsl_spline_eval( corethermalparams.Q_nuCore, tparams.Tcore, corethermalparams.Q_nuCore_acc );

    //Impose boundary conditions on T. Sets core temperature equal to initial (uniform) temperature, sets outer temperature using heat-blanketing envelope model. Needs thermal conductivity, so
    //must be done after these are initially computed
    T_BoundaryConditions(T, q, tparams.Tcore, transCoeffs, magCoeffs, tparams, domain, N_GC, process);
    exchng2Scalar(T, N_GC, process, stridetype_Sca);
    Compute_q(T, q, transCoeffs, magCoeffs, N_GC, domain, process); //compute initial heat flux
    Compute_J(B, J, T, magCoeffs, N_GC, domain, process); //compute initial current density
    if( nbrright < 0 ){
        for(size_t i=N_GC-1; i<B.shape()[1]-N_GC+1; i++){
            J[2][i][B.shape()[2]-N_GC] = 0.; //set Jphi on theta=pi boundary equal to zero.
        }
    }
    Compute_E(B, B, phE, cE, J, transCoeffs, magCoeffs, N_GC, domain, process); //compute initial electric field
    E_BoundaryConditions(phE, cE, bparams, N_GC, process, domain); //impose boundary conditions on initial electric field

    //Determine time step using cell centre spacing and Hall diffusivity
    double min_Deltar = *std::min_element( Deltar.begin(), Deltar.end() );
    double DeltaL = 1./std::sqrt( 1./(min_Deltar*min_Deltar) + 1./(r_min*r_min*Deltatheta*Deltatheta) ); //spatial step in reduced units

    std::vector<double> max_eta_Vec(num_procs); //Vector of the maximum values of eta_H*Bmag+eta_O in each process
    MPI_Allgather(&transCoeffs.max_eta_local, 1, MPI_DOUBLE, max_eta_Vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of eta_T on each process to all processes.
    double max_eta = *std::max_element(max_eta_Vec.data(), max_eta_Vec.data() + max_eta_Vec.size()); //maxmimum eta_H*Bmag+eta_O across whole simulation domain
    double Deltat_B = DeltaL*DeltaL/max_eta; //CFL limited time step for magnetic field evolution (Hall term)

    std::vector<double> max_eta_T_Vec(num_procs); //Vector of the maximum values of eta_T in each process
    MPI_Allgather(&max_eta_T_local, 1, MPI_DOUBLE, max_eta_T_Vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of eta_T on each process to all processes.
    double max_eta_T = *std::max_element(max_eta_T_Vec.data(), max_eta_T_Vec.data() + max_eta_T_Vec.size()); //maxmimum eta_T across whole simulation domain
    double Deltat_T = DeltaL*DeltaL/max_eta_T; //CFL limited time step for thermal evolution (heat conduction term)

    //double Tcore = tparams.T_init; //sets core redshifted temperature equal to t_init
    double k_Cnu = 1e-3;
    double Deltat_nu = k_Cnu*C_vCore/( -gsl_spline_eval( corethermalparams.Q_nuCore, tparams.Tcore, corethermalparams.Q_nuCore_acc ) );//core neutrino luminosity-limited time step

    double Deltat; //timestep determined using CFL limit
    if(tparams.fixed_T == false) Deltat = std::min( std::min(simparams.k_CB*Deltat_B, simparams.k_CT*Deltat_T), Deltat_nu ); //smallest of two CFL limited time steps, multiplied by Courant number < 1
    else Deltat = simparams.k_CB*Deltat_B;
    domain.Deltat = Deltat;

    /*
        Create h5 file
    */

    //Create directory for H5 files holding simulation data
    //Also copy SimSetup.in to this directory, so can check what simulation parameters (e.g., initial conditions) were used for any particular run
    if( world_rank == 0 ){
        std::ostringstream ndir;
        ndir << FILE_NAME + "/" + FILE_NAME + "/";
        std::filesystem::create_directories(ndir.str().c_str());

        fs::path sourceFile = "SimSetup.in";
        std::ostringstream ndir_up;
        ndir_up << FILE_NAME + "/";
        fs::path targetParent = ndir_up.str();
        auto target = targetParent / sourceFile.filename();

        try{
            fs::copy_file(sourceFile, target, fs::copy_options::overwrite_existing);
        }
        catch (std::exception& e){
            std::cout << e.what();
        }
    }
    MPI_Barrier( MPI_COMM_WORLD );

    H5::H5File file( FILE_NAME+"/"+FILE_NAME+"/"+FILE_NAME+"_"+RANK_NAME+FILE_EXT, H5F_ACC_TRUNC );

    //Write B to output H5 file and prepare to write data from additional time steps
    hsize_t B_dims[4] = {1,3,Nr,MyE-MyS};
    hsize_t B_maxdims[4] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t B_size[4] = {1,3,Nr,MyE-MyS}; //Total size of extended array. Increases with every time step. Initially should equal dimsf.
    hsize_t B_offset[4] = {0,0,0,0}; //offset to print new values of u at each timestep. Progressively gets larger in direction that output is extended.
    hsize_t B_dimsext[4] = {1,3,Nr,MyE-MyS}; //size of extended output each time it is extended. Should be the same as the original output i.e. equal to dimsf
    const size_t RANK = 4; //rank (dimension) of arrays to write to H5 file (equals 2 plus number of spatial dimensions)
    H5::DataSpace *dataspace_B = new DataSpace( RANK, B_dims, B_maxdims );

    H5::DSetCreatPropList B_prop;
    hsize_t B_chunk_dims[4] = {1,3,Nr,MyE-MyS};
    B_prop.setChunk(RANK, B_chunk_dims);
    H5::DataSet *B_set = new DataSet(file.createDataSet("B", H5::PredType::NATIVE_DOUBLE, *dataspace_B, B_prop));

    typedef VectorField::index_range range;
    VectorField::array_view<3>::type B_view = B[ boost::indices[range()][range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
    VectorField Phys_B(boost::extents[3][Nr][MyE-MyS]); //Cell-face average values of B across partial domain. Excludes ghost cells
    hsize_t B_memdims[3] = {3, Nr, MyE-MyS}; //the dimensions of each new Phys_B entry- used when writing subsequent time slices to output
    Phys_B = B_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
    B_set->write( Phys_B.data(), PredType::NATIVE_DOUBLE ); //does not write ghost cells.

    //Write temperature to output H5 file and prepare to write data from additional time steps. If tparams.fixed_T == true, then only include initial temperature profile.
    hsize_t T_dims[3] = {1,Nr,MyE-MyS};
    hsize_t T_maxdims[3] = {H5S_UNLIMITED,H5S_UNLIMITED,H5S_UNLIMITED};
    hsize_t T_size[3] = {1,Nr,MyE-MyS}; //Total size of extended array. Increases with every time step. Initially should equal dimsf.
    hsize_t T_offset[3] = {0,0,0}; //offset to print new values of u at each timestep. Progressively gets larger in direction that output is extended.
    hsize_t T_dimsext[3] = {1,Nr,MyE-MyS}; //size of extended output each time it is extended. Should be the same as the original output i.e. equal to dimsf
    const size_t RANK_T = 3; //rank (dimension) of scalar arrays to write to H5 file (equals 1 plus number of spatial dimensions)
    H5::DataSpace *dataspace_T = new DataSpace( RANK_T, T_dims, T_maxdims );

    H5::DSetCreatPropList T_prop;
    hsize_t T_chunk_dims[3] = {1,Nr,MyE-MyS};
    T_prop.setChunk(RANK_T, T_chunk_dims);
    H5::DataSet *T_set = new DataSet(file.createDataSet("T", H5::PredType::NATIVE_DOUBLE, *dataspace_T, T_prop));

    typedef VectorField::index_range range;
    VectorField::array_view<2>::type T_view = T[ boost::indices[range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
    ScalarField Phys_T(boost::extents[Nr][MyE-MyS]); //Cell-face average values of B across partial domain. Excludes ghost cells
    hsize_t T_memdims[2] = {Nr, MyE-MyS}; //the dimensions of each new Phys_T entry- used when writing subsequent time slices to output
    Phys_T = T_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
    T_set->write( Phys_T.data(), PredType::NATIVE_DOUBLE ); //does not write ghost cells.

    //Write J to output H5 file and prepare to write data from additional time steps
    H5::DataSpace *dataspace_J = new DataSpace( RANK, B_dims, B_maxdims );

    H5::DSetCreatPropList J_prop;
    hsize_t J_chunk_dims[4] = {1,3,Nr,MyE-MyS};
    J_prop.setChunk(RANK, J_chunk_dims);
    H5::DataSet *J_set = new DataSet(file.createDataSet("J", H5::PredType::NATIVE_DOUBLE, *dataspace_J, J_prop));

    //typedef VectorField::index_range range;
    VectorField::array_view<3>::type J_view = J[ boost::indices[range()][range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
    VectorField Phys_J(boost::extents[3][Nr][MyE-MyS]); //Cell-face average values of J across partial domain. Excludes ghost cells
    Phys_J = J_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
    J_set->write( Phys_J.data(), PredType::NATIVE_DOUBLE ); //does not write ghost cells.

    // Spatial coordinates (cell centers) for output H5 file
    hsize_t r_dims[2] = {1,Nr};
    hsize_t r_maxdims[2] = {1,Nr};
    const size_t RANK_coords = 2; //rank (dimension) of coordinate arrays to write to H5 file
    H5::DataSpace *dataspace_r = new DataSpace( RANK_coords, r_dims, r_maxdims );

    H5::DSetCreatPropList r_prop;
    hsize_t r_chunk_dims[2] = {1,Nr};
    r_prop.setChunk(RANK_coords, r_chunk_dims);

    H5::DataSet *r_set = new DataSet(file.createDataSet("r", H5::PredType::NATIVE_DOUBLE, *dataspace_r, r_prop));
    r_set->write( rDomain.data(), PredType::NATIVE_DOUBLE );

    hsize_t theta_dims[2] = {1,MyE-MyS};
    hsize_t theta_maxdims[2] = {1,MyE-MyS};
    H5::DataSpace *dataspace_theta = new DataSpace( RANK_coords, theta_dims, theta_maxdims );
    H5::DSetCreatPropList theta_prop;
    hsize_t theta_chunk_dims[2] = {1,MyE-MyS};
    theta_prop.setChunk(RANK_coords, theta_chunk_dims);

    H5::DataSet *theta_set = new DataSet(file.createDataSet("theta", H5::PredType::NATIVE_DOUBLE, *dataspace_theta, theta_prop));
    theta_set->write( thetaDomain.data(), PredType::NATIVE_DOUBLE );

    // Time data for output H5 file
    hsize_t t_dims[2] = {1,1};
    hsize_t t_maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t t_size[2] = {1,1};
    hsize_t t_offset[2] = {0,0};
    hsize_t t_dimsext[2] = {1,1};
    H5::DataSpace *dataspace_t = new DataSpace( RANK_coords, t_dims, t_maxdims );

    H5::DSetCreatPropList t_prop;
    hsize_t t_chunk_dims[2] = {1,1};
    t_prop.setChunk(RANK_coords, t_chunk_dims);

    H5::DataSet *t_set = new DataSet(file.createDataSet("t", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    hsize_t t_memdims[1] = {1}; //the dimensions of each new t entry- used when writing subsequent time slices to output
    t_set->write( &t, PredType::NATIVE_DOUBLE );

    //Energy conservation data for output H5 file. Only include this in root process H5 file
    //Uses same dataspace and sizes as the time data, so don't need to define separate hsize_t and DataSpace for this
    H5::DataSet *U_B_set = new DataSet(file.createDataSet("U_B", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *JH_set = new DataSet(file.createDataSet("JH", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *QJH_set = new DataSet(file.createDataSet("QJH", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *PF_set = new DataSet(file.createDataSet("PF", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *SH_set = new DataSet(file.createDataSet("SH", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *DeltaEInt_set = new DataSet(file.createDataSet("DeltaEInt", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop));
    H5::DataSet *T_s_eff_set = new DataSet(file.createDataSet("T_s_eff", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop)); //average local surface temperature in reduced units
    H5::DataSet *T_core_set = new DataSet(file.createDataSet("T_core", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop)); //uniform redshifted core temperature in reduced units
    H5::DataSet *Qnu_crust_set = new DataSet(file.createDataSet("Qnu_crust", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop)); //volume-integrated crust neutrino emissivity in reduced units
    H5::DataSet *Qnu_core_set = new DataSet(file.createDataSet("Qnu_core", H5::PredType::NATIVE_DOUBLE, *dataspace_t, t_prop)); //volume-integrated core neutrino emissivity in reduced units

    // Add attributes to h5 file
    H5::DataSpace att_space(H5S_SCALAR);
    H5::Attribute att1 = file.createAttribute( "B_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att1.write( H5::PredType::NATIVE_DOUBLE, &B_0 );
    H5::Attribute att2 = file.createAttribute( "t_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att2.write( H5::PredType::NATIVE_DOUBLE, &t_0 );
    H5::Attribute att3 = file.createAttribute( "L_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att3.write( H5::PredType::NATIVE_DOUBLE, &L_0 );
    H5::Attribute att4 = file.createAttribute( "T_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att4.write( H5::PredType::NATIVE_DOUBLE, &T_0 );
    H5::Attribute att5 = file.createAttribute( "s_0", H5::PredType::NATIVE_DOUBLE, att_space );
    att5.write( H5::PredType::NATIVE_DOUBLE, &s_0 );
    H5::Attribute att6 = file.createAttribute( "r_min", H5::PredType::NATIVE_DOUBLE, att_space );
    att6.write( H5::PredType::NATIVE_DOUBLE, &r_min );
    H5::Attribute att7 = file.createAttribute( "r_max", H5::PredType::NATIVE_DOUBLE, att_space );
    att7.write( H5::PredType::NATIVE_DOUBLE, &r_max );
    H5::Attribute att8 = file.createAttribute( "r_s", H5::PredType::NATIVE_DOUBLE, att_space );
    att8.write( H5::PredType::NATIVE_DOUBLE, &tparams.r_s );
    H5::Attribute att9 = file.createAttribute( "n_t_s", H5::PredType::NATIVE_DOUBLE, att_space );
    att9.write( H5::PredType::NATIVE_DOUBLE, &tparams.n_t_s );

    /*
       Run simulation until t >= t_max
    */

    size_t iter = 0; //iteration counter
    std::vector<double> t_next_Vec; //vector of times in reduced units at which to save data file
    save_timesCalc(t_next_Vec,t_max,simparams.saves_number); //times at which data is saved to h5 file
    double t_next = t_next_Vec[0]; //time in reduced units at which to save data file
    size_t N_t = 1; //index of current "next time" at which to save data to h5 file
    const size_t ECons_cadence = simparams.ECons_cadence; //cadence for printing energy conservation information to terminal
    MPI_Barrier( MPI_COMM_WORLD );

    double start_t = 0., end_t;
    if(world_rank == 0){
        start_t = MPI_Wtime();
    }

    //Energy conservation checking
    // U_BSum0: Initial volume-integrated energy in magnetic field in reduced units
    // U_BSum: Current volume-integrated energy in magnetic field in reduced units < U_BSum0
    // JouleSum, QuasiJouleSum, PoyntingSum: time-and volume-integrated energy losses due to Joule heating, quasi-Joule heating, and Poynting flux out of the simulation domain in reduced units. Values <0 are losses.
    // DeltaEInt: difference between (U_BSum - U_BSum0) and JouleSum + PoyntingSum in reduced units
    // T_av, T_s_eff: average redshifted temperature in simulation domain in reduced units (crude average) and average (effective) local surface temperature in reduced units
    // Qnu_CrustSum: crust volume-integrated neutrino emissivity in reduced units.
    double U_BSum0 = 0., U_BSum, JouleSum, QuasiJouleSum, PoyntingSum, SHSum, DeltaEInt = 0., JHInt = 0., PoyntingFInt, SHInt, T_av, T_s_eff, Qnu_CrustSum;
    std::vector<double> U_BSumVec(num_procs);
    std::vector<double> JouleSumVec(num_procs);
    std::vector<double> QuasiJouleSumVec(num_procs);
    std::vector<double> PoyntingSumVec(num_procs);
    std::vector<double> SHSumVec(num_procs);
    std::vector<double> T_avVec(num_procs);
    std::vector<double> T_minVec(num_procs);
    std::vector<double> T_maxVec(num_procs);
    std::vector<double> Ts4_avVec(num_procs);
    std::vector<double> Qnu_CrustVec(num_procs);

    std::vector<double> JouleSumHist;
    std::vector<double> QuasiJouleSumHist;
    std::vector<double> PoyntingSumHist;
    std::vector<double> SHSumHist;
    std::vector<double> DeltaESumHist;
    std::vector<double> Deltat_vec; //vector of time step values

    ScalarField q_SH(boost::extents[Nr+2*N_GC][MyE-MyS+2*N_GC]); //shock heating density (energy per unit volume per unit time)

    if( magCoeffs.quantization == true ){
        exchng2Vector(magCoeffs.J_B, N_GC, process, stridetype_Vec); //Exchange J_B cells in a periodic manner in the theta-direction into ghost cells
        QuasiEnergyConservation(B, phE, cE, J, q_SH, transCoeffs, magCoeffs, N_GC, domain, U_BSum, QuasiJouleSum, PoyntingSum, JouleSum, SHSum);
    }
    else EnergyConservation(B, phE, cE, J, q_SH, transCoeffs, magCoeffs, N_GC, domain, U_BSum, JouleSum, PoyntingSum, SHSum);
    double T_min, T_max;
    T_avCalc(T, N_GC, T_av, T_min, T_max);
    Qnu_CrustIntegratedCalc(thermCoeffs.q_nuCrust, N_GC, domain, Qnu_CrustSum);

    MPI_Allgather(&U_BSum, 1, MPI_DOUBLE, U_BSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&JouleSum, 1, MPI_DOUBLE, JouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
    if(magCoeffs.quantization == true) MPI_Allgather(&QuasiJouleSum, 1, MPI_DOUBLE, QuasiJouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&PoyntingSum, 1, MPI_DOUBLE, PoyntingSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&SHSum, 1, MPI_DOUBLE, SHSumVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&T_av, 1, MPI_DOUBLE, T_avVec.data(), 1, MPI_DOUBLE, comm1D);
    // MPI_Allgather(&T_min, 1, MPI_DOUBLE, T_minVec.data(), 1, MPI_DOUBLE, comm1D);
    // MPI_Allgather(&T_max, 1, MPI_DOUBLE, T_maxVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&tparams.Ts4_av, 1, MPI_DOUBLE, Ts4_avVec.data(), 1, MPI_DOUBLE, comm1D);
    MPI_Allgather(&Qnu_CrustSum, 1, MPI_DOUBLE, Qnu_CrustVec.data(), 1, MPI_DOUBLE, comm1D);

    if(world_rank == 0){
        U_BSum = std::reduce(U_BSumVec.begin(), U_BSumVec.end());
        JouleSum = std::reduce(JouleSumVec.begin(), JouleSumVec.end());
        if( magCoeffs.quantization == true ) QuasiJouleSum = std::reduce(QuasiJouleSumVec.begin(), QuasiJouleSumVec.end());
        PoyntingSum = std::reduce(PoyntingSumVec.begin(), PoyntingSumVec.end());
        SHSum = std::reduce(SHSumVec.begin(), SHSumVec.end());
        Qnu_CrustSum = std::reduce(Qnu_CrustVec.begin(), Qnu_CrustVec.end());
        T_av = std::reduce(T_avVec.begin(), T_avVec.end())/double(num_procs);
        T_s_eff = pow( std::reduce(Ts4_avVec.begin(), Ts4_avVec.end())/(4.*pi),0.25 );

        U_BSum0 = U_BSum;
        JouleSumHist.push_back(JouleSum);
        SHSumHist.push_back(SHSum);
        if( magCoeffs.quantization == true ){
            QuasiJouleSumHist.push_back(QuasiJouleSum);
            DeltaESumHist.push_back(QuasiJouleSum+PoyntingSum+SHSum);
        }
        else DeltaESumHist.push_back(JouleSum+PoyntingSum+SHSum);

        //Add energy conservation, core temperature and surface temperature data to H5 file
        U_B_set->write( &U_BSum, PredType::NATIVE_DOUBLE );
        JH_set->write( &JouleSum, PredType::NATIVE_DOUBLE );
        if( magCoeffs.quantization == true ) QJH_set->write( &QuasiJouleSum, PredType::NATIVE_DOUBLE );
        else QJH_set->write( &JouleSum, PredType::NATIVE_DOUBLE );
        PF_set->write( &PoyntingSum, PredType::NATIVE_DOUBLE );
        SH_set->write( &SHSum, PredType::NATIVE_DOUBLE );
        DeltaEInt_set->write( &DeltaEInt, PredType::NATIVE_DOUBLE );
        T_s_eff_set->write( &T_s_eff, PredType::NATIVE_DOUBLE );
        T_core_set->write( &tparams.Tcore, PredType::NATIVE_DOUBLE );
        Qnu_crust_set->write( &Qnu_CrustSum, PredType::NATIVE_DOUBLE );
        Qnu_core_set->write( &corethermalparams.Q_nuCoreVal, PredType::NATIVE_DOUBLE );
    }

    /*
        Main time-evolution simulation loop
    */
    while(t <= t_max){

        iter++;

        domain.t = t;

        RK_Step(B, T, phE, cE, J, q, B_np1, T_np1, q_SH, transCoeffs, thermCoeffs, magCoeffs, domain, process, bparams, tparams, corethermalparams, t);
        exchng2Vector(phE,N_GC,process,stridetype_Vec);
        exchng2Vector(cE,N_GC,process,stridetype_Vec);

        //Update time, magnetic field and temperature
        t = t + Deltat;
        B = B_np1;
        T = T_np1;

        // Check energy conservation
        if( magCoeffs.quantization == true ){
            QuasiEnergyConservation(B, phE, cE, J, q_SH, transCoeffs, magCoeffs, N_GC, domain, U_BSum, QuasiJouleSum, PoyntingSum, JouleSum, SHSum);
        }
        else EnergyConservation(B, phE, cE, J, q_SH, transCoeffs, magCoeffs, N_GC, domain, U_BSum, JouleSum, PoyntingSum, SHSum);
        //EnergyPrintOut(J, phE, cE, transCoeffs, magCoeffs, N_GC, domain);

        T_avCalc(T, N_GC, T_av, T_min, T_max);
        Qnu_CrustIntegratedCalc(thermCoeffs.q_nuCrust, N_GC, domain, Qnu_CrustSum);
        corethermalparams.Q_nuCoreVal = gsl_spline_eval( corethermalparams.Q_nuCore, tparams.Tcore, corethermalparams.Q_nuCore_acc );

        // EnergyConservationLocal(B, rE, thE, phE, J, DeltaE, UB_Init, JHPF_Cumulative, transCoeffs, magCoeffs, N_GC, domain, t, process, world_rank, iter, ECons_cadence, bparams);

        //Update microphysics. If evolving temperature (tparams.fixed_T == false), then update magnetization coefficients, thermodynamic coefficients and transport coefficients.
        // If not evolving temperature (tparams.fixed_T == true), then updated magnetization coefficients only.
        if( tparams.fixed_T == false ){

            if( bparams.quantization == true ){
                Omega_xyVaryingT(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, magCoeffs.M, magCoeffs.chi, magCoeffs.M_mu, thermCoeffs.c_v, magCoeffs.C_m);
                c_vfuncCalcB(T, n_e, A, Z, n_i, n_b, n_nf, mu_nf, tparams, N_GC, domain.n_t, thermCoeffs.c_v);
                exchng2Scalar(thermCoeffs.c_v,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.M,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.chi,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.M_mu,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.C_m,N_GC,process,stridetype_Sca);
            }
            else{
                Omega_xyVaryingTC_veOnly(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, thermCoeffs.c_v);
                c_vfuncCalc(T, n_e, A, Z, n_i, n_b, n_nf, mu_nf, tparams, N_GC, domain.n_t, thermCoeffs.c_v);
                exchng2Scalar(thermCoeffs.c_v,N_GC,process,stridetype_Sca);
            }

            max_eta_T_local = 0.;
            if( tparams.conductivity_anisotropy == true ){
                ConductivityCalcB(T, magCoeffs.Bmag, n_e, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, thermCoeffs.c_v, transCoeffs.eta_O_delta, transCoeffs.eta_O_perp, transCoeffs.kappa_delta, transCoeffs.kappa_perp, transCoeffs.kappa_H, max_eta_T_local, process, domain);
                q_nuCrustCalc(T, n_e, magCoeffs.Bmag, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, simparams.IMEX, thermCoeffs.q_nuCrust, thermCoeffs.dq_nuCrustdT);
                exchng2Scalar(transCoeffs.eta_O_delta,N_GC,process,stridetype_Sca);
                exchng2Scalar(transCoeffs.eta_O_perp,N_GC,process,stridetype_Sca);
                exchng2Scalar(transCoeffs.kappa_delta,N_GC,process,stridetype_Sca);
                exchng2Scalar(transCoeffs.kappa_perp,N_GC,process,stridetype_Sca);
                exchng2Scalar(transCoeffs.kappa_H,N_GC,process,stridetype_Sca);
                exchng2Scalar(thermCoeffs.q_nuCrust,N_GC,process,stridetype_Sca);
                if(simparams.IMEX == true){
                    exchng2Scalar(thermCoeffs.dq_nuCrustdT,N_GC,process,stridetype_Sca);
                }
            }
            else{
                ConductivityCalc(T, n_e, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, thermCoeffs.c_v, transCoeffs.eta_O, transCoeffs.kappa, max_eta_T_local);
                q_nuCrustCalc(T, n_e, magCoeffs.Bmag, A, Z, n_i, n_b, cparams, tparams, N_GC, domain.n_t, simparams.IMEX, thermCoeffs.q_nuCrust, thermCoeffs.dq_nuCrustdT);
                exchng2Scalar(transCoeffs.eta_O,N_GC,process,stridetype_Sca);
                exchng2Scalar(transCoeffs.kappa,N_GC,process,stridetype_Sca);
                exchng2Scalar(thermCoeffs.q_nuCrust,N_GC,process,stridetype_Sca);
                if(simparams.IMEX == true){
                    exchng2Scalar(thermCoeffs.dq_nuCrustdT,N_GC,process,stridetype_Sca);
                }
            }
        }
        else{
            if( magCoeffs.quantization == true ){
                Omega_xyFixedT(magCoeffs.Bmag, magCoeffs.mu_e, T, N_GC, domain.n_t, OmB_Interpolators, magCoeffs.M, magCoeffs.chi, magCoeffs.M_mu);
                exchng2Scalar(magCoeffs.M,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.chi,N_GC,process,stridetype_Sca);
                exchng2Scalar(magCoeffs.M_mu,N_GC,process,stridetype_Sca);
            }
        }

        MPI_Allgather(&U_BSum, 1, MPI_DOUBLE, U_BSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&JouleSum, 1, MPI_DOUBLE, JouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
        if( magCoeffs.quantization == true ) MPI_Allgather(&QuasiJouleSum, 1, MPI_DOUBLE, QuasiJouleSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&PoyntingSum, 1, MPI_DOUBLE, PoyntingSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&SHSum, 1, MPI_DOUBLE, SHSumVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&T_av, 1, MPI_DOUBLE, T_avVec.data(), 1, MPI_DOUBLE, comm1D);
        // MPI_Allgather(&T_min, 1, MPI_DOUBLE, T_minVec.data(), 1, MPI_DOUBLE, comm1D);
        // MPI_Allgather(&T_max, 1, MPI_DOUBLE, T_maxVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&tparams.Ts4_av, 1, MPI_DOUBLE, Ts4_avVec.data(), 1, MPI_DOUBLE, comm1D);
        MPI_Allgather(&Qnu_CrustSum, 1, MPI_DOUBLE, Qnu_CrustVec.data(), 1, MPI_DOUBLE, comm1D);

        if(world_rank == 0){
            U_BSum = std::reduce(U_BSumVec.begin(), U_BSumVec.end());
            JouleSum = std::reduce(JouleSumVec.begin(), JouleSumVec.end());
            if( magCoeffs.quantization == true ) QuasiJouleSum = std::reduce(QuasiJouleSumVec.begin(), QuasiJouleSumVec.end());
            PoyntingSum = std::reduce(PoyntingSumVec.begin(), PoyntingSumVec.end());
            SHSum = std::reduce(SHSumVec.begin(), SHSumVec.end());
            Qnu_CrustSum = std::reduce(Qnu_CrustVec.begin(), Qnu_CrustVec.end());
            T_av = std::reduce(T_avVec.begin(), T_avVec.end())/double(num_procs);
            T_s_eff = pow( std::reduce(Ts4_avVec.begin(), Ts4_avVec.end())/(4.*pi),0.25 ); //divides sum of T_s^4 values integrated over partial domains by total solid angle of star (i.e., a sphere 4*pi), then takes quartic root to get T_s_eff

            // T_max = *std::max_element(T_maxVec.begin(), T_maxVec.end());
            // T_min = *std::min_element(T_minVec.begin(), T_minVec.end());

            JouleSumHist.push_back(JouleSum);
            PoyntingSumHist.push_back(PoyntingSum);
            SHSumHist.push_back(SHSum);
            if( magCoeffs.quantization == true ){
                QuasiJouleSumHist.push_back(QuasiJouleSum);
                DeltaESumHist.push_back(QuasiJouleSum+PoyntingSum+SHSum);
            }
            else DeltaESumHist.push_back(JouleSum+PoyntingSum+SHSum);
            Deltat_vec.push_back(Deltat);

            //print out temperature and energy conservation status to terminal
            DeltaEInt = TrapezoidIntegrator(DeltaESumHist,Deltat_vec);
            PoyntingFInt = TrapezoidIntegrator(PoyntingSumHist,Deltat_vec);
            SHInt = TrapezoidIntegrator(SHSumHist,Deltat_vec);
            if( (iter == 1) || (iter%ECons_cadence == 0) ){
                std::cout << "t = " << std::setprecision(5) << t*t_0/yr << " yr, T_av = " << T_av*T_0 << " K, ";
                std::cout << "T_s = " << T_s_eff*T_0 << " K, ";
                // std::cout << "T_min = " << T_min*T_0 << " K, ";
                std::cout << "ΔU_B = " << (U_BSum-U_BSum0)*pow(B_0,2.)*pow(L_0,3.) << " erg, ";
                std::cout << "ΔE = " << DeltaEInt*pow(B_0,2.)*pow(L_0,3.) << " erg, ";
                if( magCoeffs.quantization == true ){
                    JHInt = TrapezoidIntegrator(JouleSumHist,Deltat_vec);
                    std::cout << "JH = " << (JHInt+SHInt)*pow(B_0,2.)*pow(L_0,3.) << " erg, ";
                }
                std::cout << "PF Int = " << PoyntingFInt*pow(B_0,2.)*pow(L_0,3.) << " erg, % error = " << abs(1.-(U_BSum-DeltaEInt)/U_BSum0)*100. << std::endl;
            }
        }
        if(simparams.divBCheck == true){
            if( (iter == 1) || (iter%ECons_cadence == 0) ) divB_Monitor(B, N_GC, domain, comm1D, world_rank, Ntheta_locs); //check that divB is still zero
        }
        //Adjust time step based on updated diffusivities
        if(tparams.fixed_T == false){
            MPI_Allgather(&transCoeffs.max_eta_local, 1, MPI_DOUBLE, max_eta_Vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of eta_T on each process to all processes.
            max_eta = *std::max_element(max_eta_Vec.data(), max_eta_Vec.data() + max_eta_Vec.size()); //maxmimum eta_H*Bmag across whole simulation domain
            Deltat_B = DeltaL*DeltaL/max_eta; //CFL limited time step for magnetic field evolution (Hall term)

            MPI_Allgather(&max_eta_T_local, 1, MPI_DOUBLE, max_eta_T_Vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of eta_T on each process to all processes.
            max_eta_T = *std::max_element(max_eta_T_Vec.data(), max_eta_T_Vec.data() + max_eta_T_Vec.size()); //maximum eta_T across whole simulation domain
            Deltat_T = DeltaL*DeltaL/max_eta_T; //CFL limited time step for thermal evolution (heat conduction term)

            MPI_Allgather(&T_max, 1, MPI_DOUBLE, T_maxVec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of T on each process to all processes.
            T_max = *std::max_element(T_maxVec.data(), T_maxVec.data() + T_maxVec.size());

            Deltat_nu = k_Cnu*C_vCore/( -gsl_spline_eval( corethermalparams.Q_nuCore, tparams.Tcore, corethermalparams.Q_nuCore_acc ) );//core neutrino luminosity-limited time step

            Deltat = std::min( std::min(simparams.k_CB*Deltat_B, simparams.k_CT*Deltat_T), Deltat_nu ); //smallest of three CFL limited time steps, multiplied by Courant number < 1
        }
        else{
            MPI_Allgather(&transCoeffs.max_eta_local, 1, MPI_DOUBLE, max_eta_Vec.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD); //Send maximum values of eta_T on each process to all processes.
            max_eta = *std::max_element(max_eta_Vec.data(), max_eta_Vec.data() + max_eta_Vec.size()); //maxmimum eta_H*Bmag across whole simulation domain
            Deltat_B = DeltaL*DeltaL/max_eta; //CFL limited time step for magnetic field evolution (Hall term)

            Deltat = simparams.k_CB*Deltat_B;
        }
        domain.Deltat = Deltat;

        if( t > t_next ){
            t_next = t_next_Vec[N_t];
            N_t += 1;

            // Add updated B values to H5 file
            B_size[0] += B_dimsext[0];
            B_size[1] = B_dims[1];
            B_size[2] = B_dims[2];
            B_size[3] = B_dims[3];
            B_set->extend(B_size);

            H5::DataSpace *B_filespace = new H5::DataSpace(B_set->getSpace());
            B_offset[0] += 1;
            B_filespace->selectHyperslab(H5S_SELECT_SET, B_dimsext, B_offset);
            H5::DataSpace *B_memspace = new H5::DataSpace( RANK-1, B_memdims, NULL );

            B_view = B[ boost::indices[range()][range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
            Phys_B = B_view; //copy view of B which removes ghost cells into smaller array which can be written into h5 file
            B_set->write( Phys_B.data(), PredType::NATIVE_DOUBLE, *B_memspace, *B_filespace );
            delete B_filespace;
            delete B_memspace;

            // Add updated J values to H5 file
            J_set->extend(B_size);

            H5::DataSpace *J_filespace = new H5::DataSpace(J_set->getSpace());
            J_filespace->selectHyperslab(H5S_SELECT_SET, B_dimsext, B_offset);
            H5::DataSpace *J_memspace = new H5::DataSpace( RANK-1, B_memdims, NULL );

            J_view = J[ boost::indices[range()][range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
            Phys_J = J_view; //copy view of J which removes ghost cells into smaller array which can be written into h5 file
            J_set->write( Phys_J.data(), PredType::NATIVE_DOUBLE, *J_memspace, *J_filespace );
            delete J_filespace;
            delete J_memspace;

            // Add updated T values to H5 file
            if(tparams.fixed_T == false){
                T_size[0] += T_dimsext[0];
                T_size[1] = T_dims[1];
                T_size[2] = T_dims[2];
                T_set->extend(T_size);

                H5::DataSpace *T_filespace = new H5::DataSpace(T_set->getSpace());
                T_offset[0] += 1;
                T_filespace->selectHyperslab(H5S_SELECT_SET, T_dimsext, T_offset);
                H5::DataSpace *T_memspace = new H5::DataSpace( RANK_T-1, T_memdims, NULL );

                T_view = T[ boost::indices[range(N_GC,Nr+N_GC)][range(N_GC,MyE-MyS+N_GC)] ];
                Phys_T = T_view; //copy view of T which removes ghost cells into smaller array which can be written into h5 file
                T_set->write( Phys_T.data(), PredType::NATIVE_DOUBLE, *T_memspace, *T_filespace );
                delete T_filespace;
                delete T_memspace;
            }

            // Add updated t value to H5 file
            t_size[0] += t_dimsext[0];
            t_size[1] = t_dims[0];
            t_set->extend(t_size);

            H5::DataSpace *t_filespace = new H5::DataSpace(t_set->getSpace());
            t_offset[0] += 1;
            t_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
            H5::DataSpace *t_memspace = new H5::DataSpace( 1, t_memdims, NULL );
            t_set->write( &t, PredType::NATIVE_DOUBLE, *t_memspace, *t_filespace );
            delete t_filespace;
            delete t_memspace;

            //Add updated energy conservation values to H5 file for root process
            if(world_rank == 0){

                U_B_set->extend(t_size);
                JH_set->extend(t_size);
                QJH_set->extend(t_size);
                PF_set->extend(t_size);
                SH_set->extend(t_size);
                DeltaEInt_set->extend(t_size);
                T_s_eff_set->extend(t_size);
                T_core_set->extend(t_size);
                Qnu_crust_set->extend(t_size);
                Qnu_core_set->extend(t_size);

                H5::DataSpace *U_B_filespace = new H5::DataSpace(U_B_set->getSpace());
                H5::DataSpace *JH_filespace = new H5::DataSpace(JH_set->getSpace());
                H5::DataSpace *QJH_filespace = new H5::DataSpace(QJH_set->getSpace());
                H5::DataSpace *PF_filespace = new H5::DataSpace(PF_set->getSpace());
                H5::DataSpace *SH_filespace = new H5::DataSpace(SH_set->getSpace());
                H5::DataSpace *DeltaEInt_filespace = new H5::DataSpace(DeltaEInt_set->getSpace());
                H5::DataSpace *T_s_eff_filespace = new H5::DataSpace(T_s_eff_set->getSpace());
                H5::DataSpace *T_core_filespace = new H5::DataSpace(T_core_set->getSpace());
                H5::DataSpace *Qnu_crust_filespace = new H5::DataSpace(Qnu_crust_set->getSpace());
                H5::DataSpace *Qnu_core_filespace = new H5::DataSpace(Qnu_core_set->getSpace());
                U_B_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                JH_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                QJH_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                PF_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                SH_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                DeltaEInt_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                T_s_eff_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                T_core_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                Qnu_crust_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                Qnu_core_filespace->selectHyperslab(H5S_SELECT_SET, t_dimsext, t_offset);
                H5::DataSpace *U_B_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *JH_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *QJH_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *PF_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *SH_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *DeltaEInt_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *T_s_eff_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *T_core_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *Qnu_crust_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                H5::DataSpace *Qnu_core_memspace = new H5::DataSpace( 1, t_dimsext, NULL );
                U_B_set->write( &U_BSum, PredType::NATIVE_DOUBLE, *U_B_memspace, *U_B_filespace );
                JH_set->write( &JouleSum, PredType::NATIVE_DOUBLE, *JH_memspace, *JH_filespace );
                if( magCoeffs.quantization == true ) QJH_set->write( &QuasiJouleSum, PredType::NATIVE_DOUBLE, *QJH_memspace, *QJH_filespace );
                else QJH_set->write( &JouleSum, PredType::NATIVE_DOUBLE, *QJH_memspace, *QJH_filespace );
                PF_set->write( &PoyntingSum, PredType::NATIVE_DOUBLE, *PF_memspace, *PF_filespace );
                SH_set->write( &SHSum, PredType::NATIVE_DOUBLE, *SH_memspace, *SH_filespace );
                DeltaEInt_set->write( &DeltaEInt, PredType::NATIVE_DOUBLE, *DeltaEInt_memspace, *DeltaEInt_filespace );
                T_s_eff_set->write( &T_s_eff, PredType::NATIVE_DOUBLE, *T_s_eff_memspace, *T_s_eff_filespace );
                T_core_set->write( &tparams.Tcore, PredType::NATIVE_DOUBLE, *T_core_memspace, *T_core_filespace );
                Qnu_crust_set->write( &Qnu_CrustSum, PredType::NATIVE_DOUBLE, *Qnu_crust_memspace, *Qnu_crust_filespace );
                Qnu_core_set->write( &corethermalparams.Q_nuCoreVal, PredType::NATIVE_DOUBLE, *Qnu_core_memspace, *Qnu_core_filespace );
                delete U_B_filespace;
                delete U_B_memspace;
                delete JH_filespace;
                delete JH_memspace;
                delete QJH_filespace;
                delete QJH_memspace;
                delete PF_filespace;
                delete PF_memspace;
                delete SH_filespace;
                delete SH_memspace;
                delete DeltaEInt_filespace;
                delete DeltaEInt_memspace;
                delete T_s_eff_filespace;
                delete T_s_eff_memspace;
                delete T_core_filespace;
                delete T_core_memspace;
                delete Qnu_crust_filespace;
                delete Qnu_crust_memspace;
                delete Qnu_core_filespace;
                delete Qnu_core_memspace;
            }
        }
    }

    B_prop.close();
    J_prop.close();
    r_prop.close();
    theta_prop.close();
    t_prop.close();
    delete dataspace_B;
    delete dataspace_J;
    delete dataspace_T;
    delete dataspace_t;
    delete dataspace_r;
    delete dataspace_theta;
    delete B_set;
    delete J_set;
    if(tparams.fixed_T == false) delete T_set;
    delete r_set;
    delete theta_set;
    delete t_set;
    delete U_B_set;
    delete JH_set;
    delete QJH_set;
    delete PF_set;
    delete SH_set;
    delete DeltaEInt_set;
    delete T_s_eff_set;
    delete T_core_set;
    delete Qnu_crust_set;
    delete Qnu_core_set;
    file.close();

    MPI_Type_free( &stridetype_Vec );
    MPI_Type_free( &stridetype_Vec_q );
    MPI_Type_free( &stridetype_Sca );

    if(world_rank == 0){
        end_t = MPI_Wtime();
        std::cout << "Solver done. Total run time is " << end_t-start_t << " s" << std::endl;
        std::cout << "Wrote simulation data to " << FILE_NAME << ".h5" << std::endl;

        //Append to the copy of SimSetup.in within the folder containing the h5 file for this run the time it took the run to complete and the number of processes used
        std::ostringstream filenameOutputFinal;
        filenameOutputFinal << FILE_NAME << "/SimSetup.in";
        std::fstream DatafileOutputFinal;
        DatafileOutputFinal.open(filenameOutputFinal.str().c_str(),std::fstream::app);
        DatafileOutputFinal << std::endl; //blank line
        DatafileOutputFinal << "Total run time was " << end_t-start_t << " s using " << num_procs << " processes" << std::endl;
        DatafileOutputFinal.close();
    }

    /* Tear down the communication infrastructure */
    MPI_Finalize();

    return 0;
}

/*
        Runge-Kutta timestep. Uses a two-step advance as described in Vigano et al Com. Phys. Commun. 183 (2012), 2042.
        Input: B, T, phE, cE, J, q: magnetic field, temperature, redshifted electric field and conjugate electric field times c (evaluated at location of J_phi), redshifted current density and heat flux times lapse function squared
               transCoeffs: Hall and Ohmic diffusivities and thermal conductivities in reduced units
               thermCoeffs: specific heat capacity and crust neutrino emissivity in reduced units
               magCoeffs: magnetization, differential magnetic susceptibility, mixed susceptibility and magnetocaloric coefficient
               dm, process: Domain and Process objects containing information about the simulation domain and the current process
               bparams: BParams object containing information about the initial magnetic field and its boundary conditions
               tparams: TParams object containing information about the initial temperature profile and its boundary conditions
               corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
               t: current time in reduced units
        Output: B_np1: updated magnetic field
                T_np1: updated temperature
*/
void RK_Step(VectorField & B, ScalarField & T, VectorField & phE, VectorField & cE, VectorField & J, VectorField & q, VectorField & B_np1, ScalarField & T_np1, ScalarField & q_SH,
    TransCoeffs & trC, ThermCoeffs & thC, MagCoeffs & mC, const Domain & dm, const Process & process, const BParams & bparams, TParams & tparams, CoreThermalParams & corethermalparams, double t)
{
    size_t Ntheta = dm.Ntheta;
    size_t N_GC = dm.N_GC;
    double Deltat = dm.Deltat;
    std::vector<int> Ntheta_locs = dm.Ntheta_locs;
    std::vector<int> starts = dm.starts;

    double B_pol; //magnetic field magnitude at poles in reduced units. Determined from process including theta=0 cell and broadcast to the others.
    if(process.nbrleft < 0) B_pol = B[0][B.shape()[1]-N_GC-1][N_GC]; //determines B_pol from process with nbrleft = -2 (i.e., process containing north pole theta=0)
    MPI_Bcast(&B_pol, 1, MPI_DOUBLE, 0, process.comm1D); //communicate value of B_pol to all processes
    tparams.B_pol = B_pol; //assign updated value of B_pol to tparams object, to be used in implementing anisotropic (B-dependent) thermal boundary condition

    //Implemented RK methods are SSP2 (Heun's method, strong stability-preserving RK2) and SSP3 (strong stability-preserving RK3)
    //IMEX methods are IMEX-SSP2
    std::string method = process.Timestep_method;

    static double MC_pref = B_0*B_0/(T_0*s_0); //prefactor to de-dimensionalize magnetocaloric term in heat equation
    static double gam = 1. - 1./sqrt(2.); //diagonal entry in Butcher tableau for IMEX_SSP2 timestepping

    if(method == "SSP2"){

        static ScalarField Qr(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), r-component
        static ScalarField Qtheta(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), theta-component
        static ScalarField Qphi(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law, phi-component (toroidal field)
        static ScalarField QT(boost::extents[T.shape()[0]][T.shape()[1]]); //RHS of generalized heat equation
        static VectorField B_1(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after first step
        static ScalarField T_1(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated T after first step
        double Tcore_1; //Updated Tcore after first time step

        //Step 1 of SSP2

        Compute_q(T, q, trC, mC, N_GC, dm, process);
        Compute_J(B, J, T, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B, J, T, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[2][i][j] = B[2][i][j] + Deltat*Qphi[i][j];
            }
        }
        Bphi_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_1, J, T, mC, N_GC, dm, process);
        Compute_E(B, B_1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);
        TEvolve(QT, T, q, J, phE, cE, q_SH, trC, thC, tparams, N_GC, dm);

        for(size_t i=1; i<B.shape()[1]-1; i++){
            for(size_t j=1; j<B.shape()[2]-1; j++){
                B_1[0][i][j] = B[0][i][j] + Deltat*Qr[i][j];
                B_1[1][i][j] = B[1][i][j] + Deltat*Qtheta[i][j];
                if(mC.quantization == true && tparams.fixed_T == false){
                    T_1[i][j] = T[i][j] + Deltat*( QT[i][j] - T[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                                            + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                                            + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] ) );
                }
                else{
                    T_1[i][j] = T[i][j] + Deltat*QT[i][j];
                }
                if(std::isnan(T_1[i][j])==true) std::cout << "At t = " << t << ", i = " << i << ", j = " << j << " T_1 is nan" << std::endl;
            }
        }
        Tcore_1 = tparams.Tcore + Deltat*TCoreEvolve(T, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams);
//        B_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm, stridetype_Vec);
        B_BoundaryConditionsIntermediate(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        //Updates magnitude of magnetic field at cell centers
        T_BoundaryConditions(T_1, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);

        exchng2Vector(phE, N_GC, process, dm.stridetype_Vec);
        exchng2Vector(cE, N_GC, process, dm.stridetype_Vec);
        exchng2Vector(J, N_GC, process, dm.stridetype_Vec);

        //Step 2 of SSP2
        Compute_J(B_1, J, T_1, mC, N_GC, dm, process);
        Compute_q(T_1, q, trC, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B_1, J, T_1, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[2][i][j] = 0.5*( B[2][i][j] + B_1[2][i][j] + Deltat*Qphi[i][j] );
            }
        }
        Bphi_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_np1, J, T_1, mC, N_GC, dm, process);
        Compute_E(B_1, B_np1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);
        TEvolve(QT, T_1, q, J, phE, cE, q_SH, trC, thC, tparams, N_GC, dm);

        for(size_t i=1; i<B.shape()[1]-1; i++){
            for(size_t j=1; j<B.shape()[2]-1; j++){
                B_np1[0][i][j] = 0.5*( B[0][i][j] + B_1[0][i][j] + Deltat*Qr[i][j] );
                B_np1[1][i][j] = 0.5*( B[1][i][j] + B_1[1][i][j] + Deltat*Qtheta[i][j] );
                if(mC.quantization == true && tparams.fixed_T == false){
                    T_np1[i][j] = 0.5*( T[i][j] + T_1[i][j] + Deltat*( QT[i][j] - T_1[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                                                       + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                                                       + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] ) ) );
                }
                else if(tparams.fixed_T == false){
                    T_np1[i][j] = 0.5*( T[i][j] + T_1[i][j] + Deltat*QT[i][j] );
                }
                else T_np1[i][j] = T[i][j];
            }
        }
        tparams.Tcore = 0.5*( tparams.Tcore + Tcore_1 + Deltat*TCoreEvolve(T_1, Tcore_1, trC, mC, N_GC, dm, process, corethermalparams) );

        B_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        T_BoundaryConditions(T_np1, q, tparams.Tcore, trC, mC, tparams, dm, N_GC, process);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process); //Updates magnitude of magnetic field at cell centers
        else Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process); //Updates magnitude of magnetic field at cell centers
        exchng2Vector(B_np1, N_GC, process, dm.stridetype_Vec);
        exchng2Scalar(T_np1, N_GC, process, dm.stridetype_Sca);

    }

    else if(method == "SSP3"){

        static ScalarField Qr(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), r-component
        static ScalarField Qtheta(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), theta-component
        static ScalarField Qphi(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law, phi-component (toroidal field)
        static ScalarField QT(boost::extents[T.shape()[0]][T.shape()[1]]); //RHS of generalized heat equation
        static VectorField B_1(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after first step
        static ScalarField T_1(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated T after first step
        static VectorField B_2(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after second step
        static ScalarField T_2(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated T after second step
        static double Tcore_1; //Updated Tcore after first time step
        static double Tcore_2; //Updated Tcore after second time step

        Compute_q(T, q, trC, mC, N_GC, dm, process);
        Compute_J(B, J, T, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B, J, T, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[2][i][j] = B[2][i][j] + Deltat*Qphi[i][j];
            }
        }
        Bphi_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_1, J, T, mC, N_GC, dm, process);
        Compute_E(B, B_1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);
        TEvolve(QT, T, q, J, phE, cE, q_SH, trC, thC, tparams, N_GC, dm);

        for(size_t i=1; i<B.shape()[1]-1; i++){
            for(size_t j=1; j<B.shape()[2]-1; j++){
                B_1[0][i][j] = B[0][i][j] + Deltat*Qr[i][j];
                B_1[1][i][j] = B[1][i][j] + Deltat*Qtheta[i][j];
                if(mC.quantization == true && tparams.fixed_T == false){
                    T_1[i][j] = T[i][j] + Deltat*( QT[i][j] - T[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                                            + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                                            + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] ) );
                }
                else{
                    T_1[i][j] = T[i][j] + Deltat*QT[i][j];
                }
            }
        }
        Tcore_1 = tparams.Tcore + Deltat*TCoreEvolve(T, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams);
//        B_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        B_BoundaryConditionsIntermediate(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        T_BoundaryConditions(T_1, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);

        Compute_J(B_1, J, T_1, mC, N_GC, dm, process);
        Compute_q(T_1, q, trC, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B_1, J, T_1, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_2[2][i][j] = 0.25*( 3.*B[2][i][j] + B_1[2][i][j] + Deltat*Qphi[i][j]);
            }
        }
        Bphi_BoundaryConditions(B_2, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_2, J, T_1, mC, N_GC, dm, process);
        Compute_E(B_1, B_2, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);
        TEvolve(QT, T_1, q, J, phE, cE, q_SH, trC, thC, tparams, N_GC, dm);

        for(size_t i=1; i<B.shape()[1]-1; i++){
            for(size_t j=1; j<B.shape()[2]-1; j++){
                B_2[0][i][j] = 0.25*( 3.*B[0][i][j] + B_1[0][i][j] + Deltat*Qr[i][j]);
                B_2[1][i][j] = 0.25*( 3.*B[1][i][j] + B_1[1][i][j] + Deltat*Qtheta[i][j]);
                if(mC.quantization == true && tparams.fixed_T == false){
                    T_2[i][j] = 0.25*( 3.*T[i][j] + T_1[i][j] + Deltat*( QT[i][j] - T_1[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                                                                    + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                                                                    + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] ) ) );
                }
                else{
                    T_2[i][j] = 0.25*( 3.*T[i][j] + T_1[i][j] + Deltat*QT[i][j] );
                }
            }
        }
        Tcore_2 = 0.25*( 3.*tparams.Tcore + Tcore_1 + Deltat*TCoreEvolve(T_1, Tcore_1, trC, mC, N_GC, dm, process, corethermalparams));
//        B_BoundaryConditions(B_2, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm, stridetype_Vec);
        B_BoundaryConditionsIntermediate(B_2, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_2, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_2, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        T_BoundaryConditions(T_2, q, Tcore_2, trC, mC, tparams, dm, N_GC, process);

        Compute_J(B_2, J, T_2, mC, N_GC, dm, process);
        Compute_q(T_2, q, trC, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B_2, J, T_2, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[2][i][j] = 1./3.*( B[2][i][j] + 2.*B_2[2][i][j] + 2.*Deltat*Qphi[i][j] );
            }
        }
        Bphi_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_np1, J, T_2, mC, N_GC, dm, process);
        Compute_E(B_2, B_np1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);
        TEvolve(QT, T_2, q, J, phE, cE, q_SH, trC, thC, tparams, N_GC, dm);

        for(size_t i=1; i<B.shape()[1]-1; i++){
            for(size_t j=1; j<B.shape()[2]-1; j++){
                B_np1[0][i][j] = 1./3.*( B[0][i][j] + 2.*B_2[0][i][j] + 2.*Deltat*Qr[i][j] );
                B_np1[1][i][j] = 1./3.*( B[1][i][j] + 2.*B_2[1][i][j] + 2.*Deltat*Qtheta[i][j] );
                if(mC.quantization == true && tparams.fixed_T == false){
                    T_np1[i][j] = 1./3.*( T[i][j] + 2.*T_2[i][j] + 2.*Deltat*( QT[i][j] - T_2[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                                                               + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                                                               + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] ) ) );
                }
                else if(tparams.fixed_T == false){
                    T_np1[i][j] = 1./3.*( T[i][j] + 2.*T_2[i][j] + 2.*Deltat*QT[i][j] );
                }
                else T_np1[i][j] = T[i][j];
            }
        }
        tparams.Tcore = 1./3.*( tparams.Tcore + 2.*Tcore_2 + 2.*Deltat*TCoreEvolve(T_2, Tcore_2, trC, mC, N_GC, dm, process, corethermalparams) );

        B_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        T_BoundaryConditions(T_np1, q, tparams.Tcore, trC, mC, tparams, dm, N_GC, process);
        exchng2Vector(B_np1, N_GC, process, dm.stridetype_Vec);
        exchng2Scalar(T_np1, N_GC, process, dm.stridetype_Sca);

    }

    else if(method == "IMEX-SSP2"){

        static ScalarField Qr(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), r-component
        static ScalarField Qtheta(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law (EMF), theta-component
        static ScalarField Qphi(boost::extents[B.shape()[1]][B.shape()[2]]); //RHS of Faraday's Law, phi-component (toroidal field)
        static VectorField B_1(boost::extents[3][B.shape()[1]][B.shape()[2]]); //Updated B after first step
        static ScalarField T_star(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated T after first IMEX step without implicit correction
        static ScalarField T_1(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated and corrected T after first IMEX step
        static ScalarField T_2(boost::extents[T.shape()[0]][T.shape()[1]]); //Updated and corrected T after second IMEX step
        static ScalarField QT1(boost::extents[T.shape()[0]][T.shape()[1]]); //RHS of heat equation evaluated with updated and corrected T_1
        static ScalarField QT2(boost::extents[T.shape()[0]][T.shape()[1]]); //RHS of heat equation evaluated with updated and corrected T_2
        static ScalarField FT1(boost::extents[T.shape()[0]][T.shape()[1]]); //RHS of heat equation-magnetocaloric term
        static ScalarField KT(boost::extents[T.shape()[0]][T.shape()[1]]); //K_T (explicit term of stiff part in thermal evolution)
        static ScalarField T_int(boost::extents[T.shape()[0]][T.shape()[1]]); //T_int (intermediate temperature update in ADI splitting)

        double Tcore_star, Tcore_1, Tcore_2, QTCore1, QTCore2; //IMEX-intermediate Tcore, updated Tcore after two IMEX steps, and RHS's for Tcore update

        //Coefficients in Butcher tableau for IMEX-SSP2
        static double atilde21 = 1.;
        static double a11 = gam, a21 = 1.-2.*gam, a22 = gam;
        static double wtilde1 = 0.5, wtilde2 = 0.5, w1 = 0.5, w2 = 0.5;

        //Step 1 of IMEX-SSP2
        Compute_J(B, J, T, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B, J, T, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_1[2][i][j] = ( B[2][i][j] + Deltat*Qphi[i][j] );
            }
        }
        Bphi_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_1, J, T, mC, N_GC, dm, process);
        Compute_E(B, B_1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);

        // Compute_KT_IMEX(KT, T, q, J, phE, cE, q_SH, trC, thC, mC, tparams, N_GC, dm, process);
        // Compute_Tint_IMEX(T_int, T, T, KT, thC, trC, mC, tparams, corethermalparams, N_GC, dm, process, tparams.Tcore, a11); //ADI-intermediate temperature (r-direction)
        // T_BoundaryConditions(T_int, q, tparams.Tcore, trC, mC, tparams, dm, N_GC, process);
        // Compute_T_IMEX(T_1, T_int, T, thC, trC, mC, process, N_GC, dm, a11); //ADI-final temperature (theta-direction)
        // Tcore_1 = TCoreEvolveImplicit(T, tparams.Tcore, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams, a11);

        Compute_KT_IMEX(KT, T, q, J, phE, cE, q_SH, trC, thC, mC, tparams, N_GC, dm, process);
        Tcore_1 = TCoreEvolveImplicit(T, tparams.Tcore, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams, a11);
        Compute_Tint_IMEX(T_int, T, T, KT, thC, trC, mC, tparams, corethermalparams, N_GC, dm, process, Tcore_1, a11); //ADI-intermediate temperature (r-direction)
        T_BoundaryConditions(T_int, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);
        Compute_T_IMEX(T_1, T_int, T, thC, trC, mC, process, N_GC, dm, a11); //ADI-final temperature (theta-direction)

//        T_BoundaryConditions(T_1, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);

        for(size_t i=N_GC-1; i<B.shape()[1]-(N_GC-1); i++){
            for(size_t j=N_GC-1; j<B.shape()[2]-(N_GC-1); j++) {
                B_1[0][i][j] = B[0][i][j] + Deltat*Qr[i][j];
                B_1[1][i][j] = B[1][i][j] + Deltat*Qtheta[i][j];
                if( mC.quantization == true ){
                    FT1[i][j] = -T_1[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] );
                }
                else FT1[i][j]= 0.;
                QT1[i][j] = wtilde1*FT1[i][j] + w1*( T_1[i][j] - T[i][j] )/(a11*Deltat);
            }
        }

        // B_BoundaryConditions(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        B_BoundaryConditionsIntermediate(B_1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        T_BoundaryConditions(T_1, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);

        //Step 2 of IMEX-SSP2
        Compute_J(B_1, J, T_1, mC, N_GC, dm, process);
        B_torEvolve(Qphi, q_SH, B_1, J, T_1, trC, mC, N_GC, dm, process);
        for(size_t i=0; i<B.shape()[1]; i++){
            for(size_t j=0; j<B.shape()[2]; j++){
                B_np1[2][i][j] = 0.5*( B[2][i][j] + B_1[2][i][j] + Deltat*Qphi[i][j] );
            }
        }
        Bphi_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        Compute_J_Poloidal(B_np1, J, T_1, mC, N_GC, dm, process);
        Compute_E(B_1, B_np1, phE, cE, J, trC, mC, N_GC, dm, process);
        E_BoundaryConditions(phE, cE, bparams, N_GC, process, dm);
        Compute_EMF(Qr, Qtheta, phE, N_GC, dm, process);

        Compute_KT_IMEX(KT, T_1, q, J, phE, cE, q_SH, trC, thC, mC, tparams, N_GC, dm, process);
        Compute_Tstar_IMEX(T_star, T, T_1, q, J, phE, cE, q_SH, thC, trC, mC, tparams, N_GC, dm, process, a21); //IMEX-intermediate temperature
        Tcore_star = TCore_starEvolve(T, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams, a21); //IMEX-intermediate core temperature
        Compute_Tint_IMEX(T_int, T_1, T_star, KT, thC, trC, mC, tparams, corethermalparams, N_GC, dm, process, Tcore_1, a22); //ADI-intermediate temperature (r-direction)
        T_BoundaryConditions(T_int, q, Tcore_1, trC, mC, tparams, dm, N_GC, process);
        Compute_T_IMEX(T_2, T_int, T_1, thC, trC, mC, process, N_GC, dm, a22); //ADI-final temperature (theta-direction)
        Tcore_2 = TCoreEvolveImplicit(T_1, Tcore_1, Tcore_star, trC, mC, N_GC, dm, process, corethermalparams, a22);
        // Compute_KT_IMEX(KT, T_1, q, J, phE, cE, q_SH, trC, thC, mC, tparams, N_GC, dm, process);
        // Compute_Tstar_IMEX(T_star, T, T_1, q, J, phE, cE, q_SH, thC, trC, mC, tparams, N_GC, dm, process, a21); //IMEX-intermediate temperature
        // Tcore_star = TCore_starEvolve(T, tparams.Tcore, trC, mC, N_GC, dm, process, corethermalparams, a21); //IMEX-intermediate core temperature
        // Tcore_2 = TCoreEvolveImplicit(T_1, Tcore_1, Tcore_star, trC, mC, N_GC, dm, process, corethermalparams, a22);
        // Compute_Tint_IMEX(T_int, T_1, T_star, KT, thC, trC, mC, tparams, corethermalparams, N_GC, dm, process, Tcore_2, a22); //ADI-intermediate temperature (r-direction)
        // T_BoundaryConditions(T_int, q, Tcore_2, trC, mC, tparams, dm, N_GC, process);
        // Compute_T_IMEX(T_2, T_int, T_1, thC, trC, mC, process, N_GC, dm, a22); //ADI-final temperature (theta-direction)

//        T_BoundaryConditions(T_2, q, Tcore_2, trC, mC, tparams, dm, N_GC, process);

        for(size_t i=N_GC-1; i<B.shape()[1]-(N_GC-1); i++){
            for(size_t j=N_GC-1; j<B.shape()[2]-(N_GC-1); j++){
                B_np1[0][i][j] = 0.5*( B[0][i][j] + B_1[0][i][j] + Deltat*Qr[i][j] );
                B_np1[1][i][j] = 0.5*( B[1][i][j] + B_1[1][i][j] + Deltat*Qtheta[i][j] );

                if( tparams.fixed_T == false ){
                    QT2[i][j] = w2*( ( ( T_2[i][j] - T[i][j] ) - a21/a11*( T_1[i][j] - T[i][j] ) )/(a22*Deltat) - atilde21/a22*FT1[i][j] );
                    if( mC.quantization == true ){
                        QT2[i][j] += -wtilde2*T_2[i][j]*MC_pref*mC.C_m[i][j]/thC.c_v[i][j]*( 0.5*(mC.thBhat[0][i][j]*Qr[i][j]+mC.thBhat[0][i+1][j]*Qr[i+1][j])
                                                                                    + 0.5*(mC.rBhat[1][i][j]*Qtheta[i][j]+mC.rBhat[1][i][j+1]*Qtheta[i][j+1])
                                                                                    + 0.25*(mC.rBhat[2][i][j]+mC.rBhat[2][i][j+1]+mC.thBhat[2][i][j]+mC.thBhat[2][i+1][j])*Qphi[i][j] );
                    }
                    T_np1[i][j] = T[i][j] + Deltat*( QT1[i][j] + QT2[i][j] );
                }
                else T_np1[i][j] = T[i][j];
            }
        }
        if( tparams.fixed_T == false ){
            QTCore1 = w1*( Tcore_1 - tparams.Tcore )/(a11*Deltat);
            QTCore2 = w2*( ( Tcore_2 - tparams.Tcore ) - a21/a11*( Tcore_1 - tparams.Tcore ) )/(a22*Deltat);
            tparams.Tcore = tparams.Tcore + Deltat*( QTCore1 + QTCore2 );
        }

        if(dm.InitiallyEquatoriallySymmetric == true){
            B_Symmetrize(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
            // MPI_Barrier( MPI_COMM_WORLD );
        }
        B_BoundaryConditions(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(dm.InitiallyEquatoriallySymmetric == true) B_Symmetrize(B_np1, bparams, Ntheta, N_GC, process, Ntheta_locs, starts, dm);
        if(trC.conductivity_anisotropy == true) Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O_perp, trC.max_eta_local, dm, N_GC, process);
        else Bmag_and_BhatCalc(B_np1, mC, trC.eta_H, trC.eta_O, trC.max_eta_local, dm, N_GC, process);
        T_BoundaryConditions(T_np1, q, tparams.Tcore, trC, mC, tparams, dm, N_GC, process);
        exchng2Vector(B_np1, N_GC, process, dm.stridetype_Vec);
        exchng2Scalar(T_np1, N_GC, process, dm.stridetype_Sca);
    }

    else{
        std::cout << "Invalid Runge-Kutta timestepping scheme order" << std::endl;
    }

    return;

}
