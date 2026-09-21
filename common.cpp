#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <mpi.h>
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "common.h"
#include "microphysics.h"

constexpr double pi = 3.141592653589793238463;
constexpr double c = 29979245800; //speed of light in cm/s
constexpr double hbarc = 197.3269804; //hbar times c in MeV*fm
constexpr double unit_e = 4.80320425e-10; //elementary charge in statcoulomb
constexpr double yr = 3600*24*365; //1 year in seconds
constexpr double G = 6.67430e-8; //gravitational constant in dyn*cm^2/g^2
constexpr double k_B = 8.61733034e-11; //Boltzmann constant in MeV/K
constexpr double k_Bcgs = 1.380649e-16; //Boltzmann constant in cgs units (erg/K)
constexpr double M_e = 0.51099895000; //electron mass in MeV
constexpr double M_m = 105.658375; //muon mass in MeV
constexpr double M_N = 938.92; //mean nucleon mass in MeV (938.91875434 MeV)
constexpr double M_n = 939.56542052; //(bare) neutron mass in MeV
constexpr double M_p = 938.27208943; //(bare) proton mass in MeV
constexpr double M_solar = 1.98847e33; //solar mass in g
constexpr double alpha_e = 0.0072973525693; //electromagnetic fine structure constant (dimensionless)
constexpr double B_crit = 4.41400564e13; //quantum critical magnetic field (Schwinger field) in G
constexpr double eB_crit = M_e*M_e; //critical magnetic field times elementary charge in MeV^2
constexpr double sigma_SB = 5.670374419e-5; //Stefan-Boltzmann constant in erg/cm^2/s/K^4

//1 sqrt(MeV/fm^3) = 1 sqrt(  ) = 4.002719869e16 G
constexpr double MeVtoErg = 1/6.24150907e5; //conversion factor from MeV to erg
constexpr double statCGtoMeV2 = 1.2316185036858e-5; //6.241509074e-8*hbarc; //conversion factor from statCoulomb*Gauss to MeV^2: 1 statCoulomb*Gauss = 1 erg/cm times hbar*c in erg*cm converted to MeV^2
constexpr double MeV2toGauss = 1.444027592e13; //conversion factor from MeV^2 to Gauss: 1 MeV^2 = 4.002719869e16/(197.3269804)^(3/2) G = 1.444027592e13 G (hbarc=c=1 all energy units). When using this, express e as powers of alpha_e

constexpr double B_0 = 1e15; //characteristic magnetic field (G)
constexpr double n_e0 = 1e-4; //characteristic electron density (fm^{-3})
constexpr double L_0 = 1e5; //characteristic length scale (cm)
constexpr double t_0 = 4.*pi*unit_e*n_e0*1e39*L_0*L_0/(B_0*c); //characteristic timescale (s). Taken as the Hall time-scale with length scale L_0 and field B_0
constexpr double T_0 = 1e8; //characteristic temperature (K)
constexpr double s_0 = 1e18; //characteristic entropy density (erg/K/cm^3)
constexpr double E_0 = B_0*L_0/(c*t_0); //characteristic electric field (statV/cm)

/*
    Load simulation parameters from SimProfile.in
    Inputs: world_rank: rank of current process
    Output: params: object of class SimParams containing general parameters of simulation
            Bparams: object of BParams class containing initial field configuration and magnetic field boundary conditions (see common.h for full list)
            Tparams: object of class TParams containing initial temperature configuration and thermal boundary conditions (see common.h for full list)
*/
void load_params(SimParams & params, BParams & Bparams, TParams & Tparams, ConductParams & Cparams, int world_rank)
{
    std::ifstream file("SimSetup.in");
    std::string line, entry;

    std::vector<std::string> labels;
    std::vector<std::string> values;

    while(std::getline(file, line)){
        if( line.rfind("#",0) != 0 && line.size() != 0 ){ //skips comments starting with "#" and blank lines
            std::vector<std::string> current_entries;
            std::string::iterator end_pos = std::remove(line.begin(), line.end(), ' ');
            line.erase(end_pos, line.end()); //removes white space from "line"
            std::stringstream current_line( line );
            while(std::getline(current_line, entry, ':')){
                current_entries.push_back( entry );
            }
            labels.push_back( current_entries[0] ); //label to left of colon
            values.push_back( current_entries[1] ); //value to right of colon
        }
    }
    file.close();

    //Include more substrings here if using different equations of state
    //Only BSk24 and SLy4 currently implemented
    std::string subBSk24 = "BSk24";
    std::string subSLy4 = "SLy4";

    //Set values of parameters in params object
    for(size_t j=0; j<labels.size(); j++){

        if(labels[j] == "OutputFile"){
            params.OutputFile = values[j];
        }
        if(labels[j] == "CrustEOS"){
            params.CrustEOS = values[j];

            // Fill EOS parameter of Cparams
            unsigned int resBSk24 = values[j].find(subBSk24);
            unsigned int resSLy4 = values[j].find(subSLy4);
            if(resBSk24 != std::string::npos)
                Cparams.EOS = subBSk24;
            else if(resSLy4 != std::string::npos)
                Cparams.EOS = subSLy4;
            else{
                if(world_rank == 0) std::cout << "USING UNRECOGNIZED CRUST EOS!";
            }
        }
        if(labels[j] == "CoreEOS")
            params.CoreEOS  = values[j];
        if(labels[j] == "Timestep_method"){
            params.Timestep_method = values[j];
            if(params.Timestep_method == "IMEX-SSP2" || params.Timestep_method == "IMEX-SSP3")
                params.IMEX = true;
            else params.IMEX = false;
        }
        if(labels[j] == "varying_mesh"){
            if(values[j] == "true")
                params.varying_mesh = true;
            else params.varying_mesh = false;
        }
        if(labels[j] == "Nr"){
            params.Nr = stoul(values[j]);
        }
        if(labels[j] == "Ntheta"){
            params.Ntheta = stoul(values[j]);
        }
        if(labels[j] == "r_min"){
            params.r_min = stod(values[j])/L_0; //Converts from cm to reduced units
        }
        if(labels[j] == "r_max"){
            params.r_max = stod(values[j])/L_0; //Converts from cm to reduced units
        }
        if(labels[j] == "theta_min"){
            params.theta_min = stod(values[j])*pi/180.; //Converts from degrees to radians
        }
        if(labels[j] == "theta_max"){
            params.theta_max = stod(values[j])*pi/180.; //Converts from degrees to radians
        }
        if(labels[j] == "t_max"){
            params.t_max = stod(values[j])*yr/t_0; //Converts from yr to reduced units
        }
        if(labels[j] == "k_CB"){
            params.k_CB = stod(values[j]);
        }
        if(labels[j] == "k_CT"){
            params.k_CT = stod(values[j]);
        }
        if(labels[j] == "rho_cutoff"){
            params.rho_cutoff = stod(values[j]);
        }
        if(labels[j] == "saves_number"){
            params.saves_number = stoul(values[j]);
        }
        if(labels[j] == "ECons_cadence"){
            params.ECons_cadence = stoul(values[j]);
        }
        if(labels[j] == "divBCheck"){
            if(values[j] == "true")
                params.divBCheck = true;
            else params.divBCheck = false;
        }

        if(labels[j] == "B_pol_init"){
            Bparams.B_pol_init = stod(values[j])/B_0; //Converts from G to reduced units
            if(world_rank == 0) std::cout << "Initial poloidal magnetic field: dipole with " << Bparams.B_pol_init*B_0 << " G at poles" << std::endl;
        }
        if(labels[j] == "B_tor_init"){
            Bparams.B_tor_init = stod(values[j])/B_0; //Converts from G to reduced units
            if(world_rank == 0) std::cout << "Initial toroidal magnetic field with " << Bparams.B_tor_init*B_0 << " G maximum strength" << std::endl;
        }
        if(labels[j] == "InitiallyEquatoriallySymmetric"){
            if(values[j] == "true")
                params.InitiallyEquatoriallySymmetric = true;
            else params.InitiallyEquatoriallySymmetric = false;
        }
        if(labels[j] == "T_init"){
            Tparams.T_init = stod(values[j])/T_0; //Converts from K to reduced units
            if(world_rank == 0) std::cout << "Initial temperature of " << Tparams.T_init*T_0 << " K" << std::endl;
        }
        if(labels[j] == "uniform_T"){
            if(values[j] == "true")
                Tparams.uniform_T = true;
            else Tparams.uniform_T = false;
        }
        //Whether to allow the temperature to vary (false) or not (true)
        if(labels[j] == "fixed_T"){
            if(values[j] == "true"){
                Tparams.fixed_T = true;
                if(world_rank == 0) std::cout << "FIXED temperature simulation" << std::endl;
            }
            else{
                Tparams.fixed_T = false;
                if(world_rank == 0) std::cout << "EVOLVING temperature simulation" << std::endl;
            }
        }
        //Whether to turn on (true) or not (false) the Joule heating term in the heat equation
        if(labels[j] == "JouleHeating"){
            if(values[j] == "true"){
                Tparams.JouleHeating = 1.;
                if(world_rank == 0) std::cout << "JOULE HEATING turned ON" << std::endl;
            }
            else{
                Tparams.JouleHeating = 0.;
                if(world_rank == 0) std::cout << "JOULE HEATING turned OFF" << std::endl;
            }
        }
        //Whether to include Landau quantization effects (true) or not (false)
        if(labels[j] == "quantization"){
            if(values[j] == "true"){
                Bparams.quantization = true;
                Cparams.quantization = true;
                if(world_rank == 0) std::cout << "Landau quantization effects turned ON" << std::endl;
            }
            else{
                Bparams.quantization = false;
                Cparams.quantization = false;
                if(world_rank == 0) std::cout << "Landau quantization effects turned OFF" << std::endl;
            }
        }
        //Whether to include anisotropic conductivity (true) or not (false)
        if(labels[j] == "conductivity_anisotropy"){
            if(values[j] == "true"){
                Tparams.conductivity_anisotropy = true;
                if(world_rank == 0) std::cout << "ANISOTROPIC transport coefficients" << std::endl;
            }
            else{
                Tparams.conductivity_anisotropy = false;
                if(world_rank == 0) std::cout << "ISOTROPIC transport coefficients" << std::endl;
            }
        }

        //Whether to include general relativity (true) or not (false)
        if(labels[j] == "GR"){
            if(values[j] == "true"){
                params.GR = true;
                if(world_rank == 0) std::cout << "General relativity turned ON" << std::endl;
            }
            else{
                params.GR = false;
                if(world_rank == 0) std::cout << "General relativity turned OFF" << std::endl;
            }
        }

        //Whether to include superfluidity (true) or not (false)
        if(labels[j] == "SF"){
            if(values[j] == "true"){
                Tparams.SF = true;
                if(world_rank == 0) std::cout << "Superfluidity turned ON" << std::endl;
            }
            else{
                Tparams.SF = false;
                if(world_rank == 0) std::cout << "Superfluidity turned OFF" << std::endl;
            }
        }

        //Whether to include in-medium corrections to modified Urca neutrino emissivities (true) or not (false)
        if(labels[j] == "InMediumMUrca"){
            if(values[j] == "true"){
                Tparams.InMediumMUrca = true;
                if(world_rank == 0) std::cout << "Using in-medium MUrca rates" << std::endl;
            }
            else{
                Tparams.InMediumMUrca = false;
                if(world_rank == 0) std::cout << "NO in-medium MUrca corrections" << std::endl;
            }
        }

        //Whether to include Landau quantization effects in current at polar axis (true) or not (false)
        if(labels[j] == "PolarAxisQuantization"){
            if(values[j] == "true"){
                Bparams.PolarAxisQuantization = true;
            }
            else{
                Bparams.PolarAxisQuantization = false;
            }
        }

        if(labels[j] == "B_perp_lower"){
            Bparams.B_perp_lower = values[j];
        }
        if(labels[j] == "B_parallel_lower"){
            Bparams.B_parallel_lower = values[j];
        }
        if(labels[j] == "B_pol_upper"){
            Bparams.B_pol_upper = values[j];
        }
        if(labels[j] == "B_tor_upper"){
            Bparams.B_tor_upper = values[j];
        }
        if(labels[j] == "E_perp_lower"){
            Bparams.E_perp_lower = values[j];
        }
        if(labels[j] == "E_parallel_lower"){
            Bparams.E_parallel_lower = values[j];
        }
        if(labels[j] == "E_perp_upper"){
            Bparams.E_perp_upper = values[j];
        }
        if(labels[j] == "E_parallel_upper"){
            Bparams.E_parallel_upper = values[j];
        }
        if(labels[j] == "T_lower"){
            Tparams.T_lower = values[j];
        }
        if(labels[j] == "T_upper"){
            Tparams.T_upper = values[j];
            if(values[j] == "magnetic_P15" && world_rank == 0) std::cout << "Magnetic field-DEPENDENT heat blanketing envelope from Potekhin+2015" << std::endl;
            else if(values[j] == "magnetic_P01" && world_rank == 0) std::cout << "Magnetic field-DEPENDENT heat blanketing envelope from Potekhin+2001" << std::endl;
            else if(values[j] == "non-magnetic_G83" && world_rank == 0) std::cout << "Magnetic field-INDEPENDENT heat blanketing envelope from Gudmundsson+1983" << std::endl;
            else if(values[j] == "non-magnetic_P01" && world_rank == 0) std::cout << "Magnetic field-INDEPENDENT heat blanketing envelope from Potekhin+2001" << std::endl;
            else if(world_rank == 0) std::cout << "!!! WARNING: INVALID HEAT BLANKETING ENVELOPE SELECTED !!!" << std::endl;
        }
        if(labels[j] == "Z_impurity"){
            Bparams.Z_impurity = values[j];
            Cparams.Z_impurity = values[j];
            if(values[j] == "zero" && world_rank == 0) std::cout << "Using impurity parameter Q=0 " << std::endl;
            else if(values[j] == "Vigano2013Q100" && world_rank == 0) std::cout << "Using impurity parameter model Q100 from Vigano 2013 PhD thesis" << std::endl;
            else if(values[j] == "Carreau2020BSk24" && world_rank == 0) std::cout << "Using impurity parameter model for BSk24 EoS from Carreau+2020" << std::endl;
        }
        if(labels[j] == "C_hyp"){
            params.C_hyp = stod(values[j]);
            if(world_rank == 0) std::cout << "Using toroidal B hyperdiffusivity prefactor of " << values[j] << std::endl;
        }

        //Determines superfluid gap model. Options are: SF_n_crust = "SFB" or "Ho2012"; SF_p_core = "CCDK" or "Ho2012"; SF_n_core = "TToa" or "SYHHP"
        //Gaps are parametrized based on Andersson, Comer and Glampedakis, Nucl. Phys. A 763, 212 (2005), Kaminker, Haensel, and Yakovlev, A&A 373, L17 (2001), and Ho et al. PRC 91, 015806 (2015) model
        if(labels[j] == "SF_n_crust"){
            Tparams.SF_n_crust = values[j];
        }
        if(labels[j] == "SF_p_core"){
            Tparams.SF_p_core = values[j];
        }
        if(labels[j] == "SF_n_core"){
            Tparams.SF_n_core = values[j];
        }

    }

    //Print error if certain entries in params are empty
    if( world_rank == 0 ){
        if( params.CrustEOS.empty() )
            std::cout << "!!! WARNING: Missing crust EOS data table name !!!" << std::endl;
        if( params.CoreEOS.empty() && Tparams.fixed_T == false )
            std::cout << "!!! WARNING: Missing core EOS data table name; needed for temperature boundary condition !!!" << std::endl;
        if( params.Timestep_method.empty() )
            std::cout << "!!! WARNING: Missing timestepping method !!!" << std::endl;
        if( params.Nr == 0 )
            std::cout << "!!! WARNING: Missing r-direction resolution !!!" << std::endl;
        if( params.Ntheta == 0)
            std::cout << "!!! WARNING: Missing theta-direction resolution !!!" << std::endl;
        if( params.r_min < 0. || params.r_max < 0. || params.theta_min > 0. || params.theta_max < 0. )
            std::cout << "!!! WARNING: Missing r_min, r_max, theta_min or theta_max !!!" << std::endl;
        if( params.t_max < 1e-20 )
            std::cout << "!!! WARNING: Missing maximum run time !!!" << std::endl;
        if( params.k_CB < 1e-20 )
            std::cout << "!!! WARNING: Missing Courant number for magnetic evolution !!!" << std::endl;
        if( params.k_CT < 1e-20 )
            std::cout << "!!! WARNING: Missing Courant number for thermal evolution !!!" << std::endl;
        if( params.rho_cutoff < 1e-20 )
            std::cout << "!!! WARNING: Missing cutoff density !!!" << std::endl;
        if( params.saves_number == 0 )
            std::cout << "!!! WARNING: Missing maximum number of snapshots to save !!!" << std::endl;
        if( params.ECons_cadence == 0 )
            std::cout << "!!! WARNING: Missing energy conservation print-out cadence !!!" << std::endl;
        if( params.OutputFile.empty() )
            std::cout << "!!! WARNING: Missing output file name !!!" << std::endl;

        if( Bparams.B_perp_lower.empty() ){
            Bparams.B_perp_lower = "perfect_conducting";
            std::cout << "!!! WARNING: Missing lower boundary condition on perpendicular magnetic field; assuming perfect conducting !!!" << std::endl;
        }
        if( Bparams.B_parallel_lower.empty() ){
            Bparams.B_parallel_lower = "perfect_conducting";
            std::cout << "!!! WARNING: Missing lower boundary condition on parallel magnetic field; assuming perfect conducting !!!" << std::endl;
        }
        if( Bparams.E_perp_lower.empty() ){
            Bparams.E_perp_lower = "perfect_conducting";
            std::cout << "!!! WARNING: Missing lower boundary condition on perpendicular electric field; assuming perfect conducting !!!" << std::endl;
        }
        if( Bparams.E_parallel_lower.empty() ){
            Bparams.E_parallel_lower = "perfect_conducting";
            std::cout << "!!! WARNING: Missing upper boundary condition on parallel electric field; assuming perfect conducting !!!" << std::endl;
        }
        if( Bparams.B_pol_upper.empty() ){
            Bparams.B_pol_upper = "vacuum";
            std::cout << "!!! WARNING: Missing upper boundary condition on poloidal magnetic field; assuming vacuum !!!" << std::endl;
        }
        if( Bparams.B_tor_upper.empty() ){
            Bparams.B_tor_upper = "vacuum";
            std::cout << "!!! WARNING: Missing upper boundary condition on toroidal magnetic field; assuming vacuum !!!" << std::endl;
        }
        if( Bparams.E_perp_upper.empty() ){
            Bparams.E_perp_upper = "vacuum";
            std::cout << "!!! WARNING: Missing upper boundary condition on perpendicular electric field; assuming vacuum !!!" << std::endl;
        }
        if( Bparams.E_parallel_upper.empty() ){
            Bparams.E_parallel_upper = "vacuum";
            std::cout << "!!! WARNING: Missing upper boundary condition on parallel electric field; assuming vacuum !!!" << std::endl;
        }
        if( Tparams.T_lower.empty() ){
            Tparams.T_lower = "fake_core";
            std::cout << "!!! WARNING: Missing lower boundary condition on temperature; assuming fake core !!!" << std::endl;
        }
        if( Tparams.T_upper.empty() ){
            Tparams.T_upper = "isotropic";
            std::cout << "!!! WARNING: Missing upper boundary condition on temperature; assuming isotropic !!!" << std::endl;
        }

    }

    //Set values of SF gap parameters using SFgaps function
    if(Tparams.SF == true){
        SFgaps(Tparams);
        if(world_rank == 0){
            std::cout << "Using pairing gap models " << Tparams.SF_n_crust << " (1S0 n crust), " << Tparams.SF_p_core << " (1S0 p, core), " << Tparams.SF_n_core << " (3P2 n, core)" << std::endl;
        }
    }

    return;

}

/*

    Partitions initial domain into partial domains to be evolved by each process, in 1 dimension (e.g., strips for a 2D domain decomposition)
    Input: N: size of domain in direction to be partitions
           num_procs: number of processes
           MyID: rank of this process
    Output: s, e: start and end indices of partition of array to be evolved by this process

*/
void MPE_Decomp1D( size_t N, int num_procs, int MyID, size_t &s, size_t &e )
{

    size_t nlocal, deficit;

    nlocal = N / num_procs;
    s = MyID * nlocal;
    deficit = N % num_procs;
    s = s + std::min(MyID,(int)deficit);

    if (MyID < (int)deficit){
      nlocal = nlocal;
    }
    e = s + nlocal;
    if( (e > N) || (MyID == num_procs-1)){
     e = N;
    }

    return;
}

/*
    Fills ghost cells by exchanging data between cells- for VectorField objects
    Input: A: local array A
           N_GC: number of ghost cells to exchange
           process: Process object containing information about the simulation domain and the current process
           stridetype: MPI_Datatype defining the stride of the data to be exchanged between parallelized domains
*/
void exchng2Vector(VectorField & A, size_t N_GC, const Process & process, MPI_Datatype stridetype)
{
    MPI_Comm comm1D = process.comm1D; //MPI communicator for decomposed domain
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t N_comps = A.shape()[0]; //number of vector field components
    size_t Nr = A.shape()[1]; //Extent of vector field in r-dimension. Includes ghost cells.
    size_t Ntheta_loc = A.shape()[2]; //Extent of vector field in theta-dimension. includes ghost cells.

    // int count = static_cast<int>(N_comps*Nr);
    // int blocklength = static_cast<int>(N_GC);
    // int stride = static_cast<int>(Ntheta_loc);
    // MPI_Datatype stridetype_int;
    // MPI_Type_vector( N_comps*Nr, N_GC, Ntheta_loc, MPI_DOUBLE, &stridetype_int);
    // MPI_Type_vector( count, blocklength, stride, MPI_DOUBLE, &stridetype_int);
    // MPI_Type_commit( &stridetype_int );
    //
    // MPI_Sendrecv( &(A[0][0][Ntheta_loc-2*N_GC]), 1, stridetype_int, nbrright, 0, &(A[0][0][0]), 1, stridetype_int, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    // MPI_Sendrecv( &(A[0][0][N_GC]), 1, stridetype_int, nbrleft, 1, &(A[0][0][Ntheta_loc-N_GC]), 1, stridetype_int, nbrright, 1, comm1D, MPI_STATUS_IGNORE);
    //
    // MPI_Type_free( &stridetype_int );

    MPI_Sendrecv( &(A[0][0][Ntheta_loc-2*N_GC]), 1, stridetype, nbrright, 0, &(A[0][0][0]), 1, stridetype, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    MPI_Sendrecv( &(A[0][0][N_GC]), 1, stridetype, nbrleft, 1, &(A[0][0][Ntheta_loc-N_GC]), 1, stridetype, nbrright, 1, comm1D, MPI_STATUS_IGNORE);

    return;
}

/*
    Fills ghost cells by exchanging data between cells- for ScalarField objects
    Input: A: local array A
           N_GC: number of ghost cells to exchange
           process: Process object containing information about the simulation domain and the current process
           stridetype: MPI_Datatype defining the stride of the data to be exchanged between parallelized domains
*/
void exchng2Scalar(ScalarField & A, size_t N_GC, const Process & process, MPI_Datatype stridetype)
{
    MPI_Comm comm1D = process.comm1D; //MPI communicator for decomposed domain
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t Nr = A.shape()[0]; //Extent of vector field in r-dimension. Includes ghost cells.
    size_t Ntheta_loc = A.shape()[1]; //Extent of vector field in theta-dimension. includes ghost cells.

    // int count = static_cast<int>(Nr);
    // int blocklength = static_cast<int>(N_GC);
    // int stride = static_cast<int>(Ntheta_loc);
    // MPI_Datatype stridetype_int;
    // MPI_Type_vector( count, blocklength, stride, MPI_DOUBLE, &stridetype_int);
    // MPI_Type_commit( &stridetype_int );
    //
    // MPI_Sendrecv( &(A[0][Ntheta_loc-2*N_GC]), 1, stridetype_int, nbrright, 0, &(A[0][0]), 1, stridetype_int, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    // MPI_Sendrecv( &(A[0][N_GC]), 1, stridetype_int, nbrleft, 1, &(A[0][Ntheta_loc-N_GC]), 1, stridetype_int, nbrright, 1, comm1D, MPI_STATUS_IGNORE);
    //
    // MPI_Type_free( &stridetype_int );

    MPI_Sendrecv( &(A[0][Ntheta_loc-2*N_GC]), 1, stridetype, nbrright, 0, &(A[0][0]), 1, stridetype, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    MPI_Sendrecv( &(A[0][N_GC]), 1, stridetype, nbrleft, 1, &(A[0][Ntheta_loc-N_GC]), 1, stridetype, nbrright, 1, comm1D, MPI_STATUS_IGNORE);

    return;
}

/*
    Fills ghost cells by exchanging data between cells- for std::vector<double> objects
    Input: A: local array A
           N_GC: number of ghost cells to exchange
           comm1D: MPI communicator for processes organized into Cartesian grid
           nbrleft, nbrright: which processes are to the left and right of this one
*/
void exchng2Array(std::vector<double> & A, size_t N_GC, const Process & process)
{
    MPI_Comm comm1D = process.comm1D; //MPI communicator for decomposed domain
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t Ntheta_loc = A.size(); //Extent of vector field in theta-dimension. includes ghost cells.

    MPI_Sendrecv( &A[Ntheta_loc-2*N_GC], 2, MPI_DOUBLE, nbrright, 0, &A[0], 2, MPI_DOUBLE, nbrleft, 0, comm1D, MPI_STATUS_IGNORE);
    MPI_Sendrecv( &A[N_GC], 2, MPI_DOUBLE, nbrleft, 1, &A[Ntheta_loc-N_GC], 2, MPI_DOUBLE, nbrright, 1, comm1D, MPI_STATUS_IGNORE);

    return;
}

/*
    Computes cadence for saving simulation snapshots to H5 file
    Inputs: t_max: maximum simulation time in reduced units
            saves_number: unsigned integer containing maximum number of snapshots saved to H5 file
    Output: t_next_Vec: vector of times after which to save simulation snapshots to H5 file
*/
void save_timesCalc(std::vector<double> & t_next_Vec, double t_max, size_t saves_number){

    double Delta_t = t_max/double(saves_number);

    //Save extra time steps at early times to better understand early time thermal evolution behaviour
    size_t ETS = 20; //number of additional early time saves minus one
    double t_extra_minimum = 1e-2*yr/t_0;
    if(Delta_t > t_extra_minimum){
        for(size_t i=1; i< ETS; i++){
                t_next_Vec.push_back( pow( 10.,log10(t_extra_minimum)+double(i)/double(ETS)*( log10(Delta_t)-log10(t_extra_minimum) ) ) );
        }
    }
    for(size_t i=1; i<= saves_number; i++){
        t_next_Vec.push_back(Delta_t*double(i));
    }

    return;
}

/*
        MC slope limiter (3 argument)
        Arguments: a, b, c: values to compare
        Output: minmod(a,b,c)
*/
// double MC(double a, double b, double c)
// {
//
//     double result = 0.;
//
//     if( a>0. && b>0. && c>0. ){
//         result = std::min( a, std::min(b,c) );
//     }
//     else if( a<0. && b<0. && c<0. ){
//         result = std::max( a, std::max(b,c) );
//     }
//
//     return result;
// }

/*
        Superbee slope limiter function
        Arguments: a, b: values to compare
        Output: maxmod(minmod(2*a,b),minmod(a,2*b))
*/
double superbee(double a, double b)
{

    double result = 0.;
    double minmod1, minmod2;

    if( a*b > 0. ){
        if( 2*abs(a) > abs(b) ){
            minmod1 = a;
        }
        else{
            minmod1 = b;
        }

        if( 2*abs(b) > abs(a) ){
            minmod2 = b;
        }
        else{
            minmod2 = a;
        }

        if( abs(minmod1) > abs(minmod2) && minmod1*minmod2 > 0 ){
            result = minmod1;
        }
        else if( abs(minmod1) < abs(minmod2) && minmod1*minmod2 > 0 ){
            result = minmod2;
        }
        else result = 0.;
    }
    else result = 0;

    return result;
}

/*
        van Leer slope limiter function
        Arguments: a, b: values to compare
        Output: van Leer limiter function phi(a/b)
*/
double vanLeer(double a, double b)
{
    double result = 0;

    if( b > 0 ){
        result = ( a + abs(a) )/( b + abs(a) );
    }
    else if( b < 0){
        result = ( a - abs(a) )/( b - abs(a) );
    }

    return result;
}

/*
        Minmod slope limiter function (2 argument)
        Arguments: a, b: values to compare
        Output: van Leer limiter function phi(a/b)
*/
// double minmod2(double a, double b)
// {
//     double result = 0;
//     double sgna, sgnb;
//     if( abs(a) > 1e-20) sgna = a/abs(a);
//     else sgna = 0.;
//     if( abs(b) > 1e-20) sgnb = b/abs(b);
//     else sgnb = 0.;
//
//     result = 0.5*( sgna + sgnb )*std::min( abs(a), abs(b) );
//
//     return result;
// }

/*
        Minmod function for MC flux limiter
        Arguments: r: slope ratio
*/
double MCFlux(double r)
{
    double result = 0.;

    result = std::max(0.,std::max( std::min(2.*r,1.),std::min(r,2.) ));

    return result;
}

/*
        Minmod function for superbee flux limiter
        Arguments: r: slope ratio
*/
double superbeeFlux(double r)
{
    double result = 0.;

    result = std::max(0.,std::min(2.*r,std::min(2.,(1.+r)/2.)));

    return result;
}

/*
        One-sided exponential limiter for outer radial boundary.
        Can be applied to smoothed non-toroidal electric field.
        Inputs: r: radius of cell center in reduced units
                r_max: outer radius of domain in reduced units
                cl: cutoff distance (from outer radius) in reduced units

 */
double OuterEdgeLimiter(double r, double r_max, double cl)
{
    return exp( - (r_max-r)/cl );
}

/*
        Double-sided exponential limiter (bump function) for radial boundaries.
        Can be applied to smooth non-toroidal electric field.
        Inputs: r: radius of cell center in reduced units
                r_max: outer radius of domain in reduced units
                r_min: inner radius of domain in reduced units
                cl: cutoff distance (from outer radius) in reduced units

 */
double EdgeLimiter(double r, double r_max, double r_min, double cl)
{
    double mp = (r_max+r_min)/2.; //midpoint radius of crust
    double hw = (r_max-r_min)/2.; //radial half-width of crust
    return 1. - exp( cl*cl/(hw*hw)*( 1./( (r-mp)*(r-mp)/(hw*hw) - 1. )  + 1. ) );
}

/*
        Double-sided exponential limiter (bump function) for theta=0,pi boundaries.
        Can be applied to smooth non-toroidal electric field.
        Inputs: theta: polar angle of cell center in radians
                cl: cutoff angle (from poles) in reduced units

 */
double PoleEdgeLimiter(double theta, double cl)
{
    double mp = pi/2; //midpoint polar angle of crust
    double hw = pi/2.; //polar half-angle of crust
    return 1. - exp( cl*cl/(hw*hw)*( 1./( (theta-mp)*(theta-mp)/(hw*hw) - 1. )  + 1. ) );
}

/*
        Trapezoid rule integration in 1D
        Inputs: A: vector of integrand sampled in 1D every dx
                dx: vector of spacing between every sample of integrand
        Output: value of integral
*/
double TrapezoidIntegrator(std::vector<double> & A, std::vector<double> & dx)
{
    double integral = 0.;
    for(size_t j=0; j<A.size()-1; j++){
        integral = integral + 0.5*(A[j+1]+A[j])*dx[j];
    }

    return integral;
}

/*
        3D dot product of vector fields defined over 2D coordinates
        Inputs: A, B: VectorFields to take dot product of
                i, j: coordinate indices of location at which dot product is evaluated
*/
double dot(VectorField & A, VectorField & B, size_t i, size_t j)
{
    return A[0][i][j]*B[0][i][j] + A[1][i][j]*B[1][i][j] + A[2][i][j]*B[2][i][j];
}

/*
        3D dot product of vector fields defined over 2D coordinates taken at location of J_r
        Inputs: A, B: VectorFields to take dot product of
                i, j: coordinate indices of location at which dot product is evaluated
*/
double rdot(VectorField & A, VectorField & B, size_t i, size_t j)
{
    return A[0][i][j]*B[0][i][j] + 0.25*(A[1][i][j]+A[1][i+1][j]+A[1][i][j-1]+A[1][i+1][j-1])*B[1][i][j] + 0.5*(A[2][i][j]+A[2][i+1][j])*B[2][i][j];
}

/*
        3D dot product of vector fields defined over 2D coordinates taken at location of J_theta
        Inputs: A, B: VectorFields to take dot product of
                i, j: coordinate indices of location at which dot product is evaluated
*/
double thdot(VectorField & A, VectorField & B, size_t i, size_t j)
{
    return 0.25*(A[0][i][j]+A[0][i][j+1]+A[0][i-1][j]+A[0][i-1][j+1])*B[0][i][j] + A[1][i][j]*B[1][i][j] + 0.5*(A[2][i][j]+A[2][i][j+1])*B[2][i][j];
}

/*
        3D dot product of vector fields defined over 2D coordinates taken at location of J_phi
        Inputs: A, B: VectorFields to take dot product of
                i, j: coordinate indices of location at which dot product is evaluated
*/
double phdot(VectorField & A, VectorField & B, size_t i, size_t j)
{
    return 0.5*(A[0][i][j]+A[0][i-1][j])*B[0][i][j] + 0.5*(A[1][i][j]+A[1][i][j-1])*B[1][i][j] + A[2][i][j]*B[2][i][j];
}

/*
        Radial component of gradient of scalar field A, evaluated with cell-centered forward difference method
        Between neighbouring cells.
        Metric coefficient not included
        Inputs: A: ScalarField to take gradient of
                Deltar: spacing in radial direction of finite difference approximation to gradient
                i, j: coordinate indices of location at which gradient is evaluated
*/
double gradr(ScalarField & A, double Deltar, size_t i, size_t j)
{
    return (A[i][j] - A[i-1][j])/Deltar;
}

/*
        Theta component of gradient of scalar field A, evaluated with cell-centered forward difference method
        Between neighbouring cells.
        Inputs: A: ScalarField to take gradient of
                rDeltatheta: denominator in angular direction of finite difference approximation to gradient. Includes r factor.
                i, j: coordinate indices of location at which gradient is evaluated
*/
double gradth(ScalarField & A, double rDeltatheta, size_t i, size_t j)
{
    return (A[i][j] - A[i][j-1])/rDeltatheta;
}

/*
        Radial component of gradient of scalar field A, evaluated with cell-centered centered difference method
        Computed with values from two cells neigbouring cell i,j in the radial direction
        Metric coefficient not included
        Inputs: A: ScalarField to take gradient of
                Deltar: spacing in radial direction of finite difference approximation to gradient
                i, j: coordinate indices of location at which gradient is evaluated
*/
double gradrCentered(ScalarField & A, double Deltar, size_t i, size_t j)
{
    return (A[i+1][j] - A[i-1][j])/Deltar;
}

/*
        Theta component of gradient of scalar field A, evaluated with cell-centered centered difference method
        Computed with values from two cells neighbouring cell i,j in the polar direction
        Inputs: A: ScalarField to take gradient of
                rDeltatheta: denominator in angular direction of finite difference approximation to gradient. Includes r factor.
                i, j: coordinate indices of location at which gradient is evaluated
*/
double gradthCentered(ScalarField & A, double rDeltatheta, size_t i, size_t j)
{
    return (A[i][j+1] - A[i][j-1])/rDeltatheta;
}

/*
        Computes magnitude of magnetic field and components of unit vector in direction of magnetic field at cell centers
        Also computes largest value of B*eta_H + eta_O across partial domain
        Computes B_r and B_theta components at cell centers by averaging two nearest values at cell edges.
        Inputs: B: VectorField containing magnetic field
                mC: MagCoeffs object containing properties of magnetic field and magnetization
                eta_H: Hall diffusivity in reduced units
                eta_O: Ohmic diffusivity in reduced units. Use eta_O_perp for this in the anisotropic conductivity case
                max_eta_local: maximum value of B*eta_H*(1-4*pi*chi) + eta_O on the local part of the partitioned domain, including correction factor for nonzero magnetic diffusivity (1-4*pi*chi) which is 1 if
                                Landau quantization is disabled. Used to determine the CFL-limited timestep.
                dm: Domain object containing information about the simulation domain
                N_GC: number of ghost cells
                process: Process objects containing information about the simulation domain and the current process
        Output: mC.Bmag: ScalarField containing magnitude of magnetic field
                mC.rBhat, mC.thBhat, mC.phBhat: VectorFields containing unit vector of magnetic field evaluated at locations of J_r, J_theta, J_phi
                max_eta_H_local: maximum value of Bmag*eta_H across partial domain
*/
void Bmag_and_BhatCalc(VectorField & B, MagCoeffs & mC, ScalarField & eta_H, ScalarField & eta_O, double & max_eta_local, const Domain & dm, size_t N_GC, const Process & process)
{
    double Bmag, Br_sq, Btheta_sq;
    double rB_r, rB_theta, rB_phi, rBmag;
    double thB_r, thB_theta, thB_phi, thBmag;
    double phB_r, phB_theta, phB_phi, phBmag;

    max_eta_local = 0.;

    for(size_t i = 1; i < B.shape()[1]-1; i++){
        for(size_t j = 1; j < B.shape()[2]-1; j++){
            if( i == B.shape()[1] - 1 ){
                Br_sq = B[0][i][j]*B[0][i][j];
            }
            else{
                Br_sq = 0.5*(B[0][i][j]*B[0][i][j] + B[0][i+1][j]*B[0][i+1][j]); //0.25*( B[0][i][j]*B[0][i][j] + 2.*B[0][i][j]*B[0][i+1][j] + B[0][i+1][j]*B[0][i+1][j] );
            }
            if( j == B.shape()[2] - 1 ){
                Btheta_sq = B[1][i][j]*B[1][i][j];
            }
            else{
                Btheta_sq = 0.5*(B[1][i][j]*B[1][i][j] + B[1][i][j+1]*B[1][i][j+1]); //0.25*( B[1][i][j]*B[1][i][j] + 2.*B[1][i][j]*B[1][i][j+1] + B[1][i][j+1]*B[1][i][j+1] );
            }
            Bmag = sqrt( Br_sq + Btheta_sq + B[2][i][j]*B[2][i][j] );
            mC.Bmag[i][j] = std::max( Bmag, 1e-10 ); //never let Bmag be zero. Very small Bmag will have the same effect as Bmag = 0, but avoid any divide by zero errors.

            //Compute components of Bhat evaluated at location of B_phi
//            if(Bmag < 1e-30){
//                mC.Bhat[0][i][j] = 0.;
//                mC.Bhat[1][i][j] = 0.;
//                mC.Bhat[2][i][j] = 0.;
//            }
//            else{
//                mC.Bhat[0][i][j] = sqrt(Br_sq)/Bmag;
//                mC.Bhat[1][i][j] = sqrt(Btheta_sq)/Bmag;
//                mC.Bhat[2][i][j] = B[2][i][j]/Bmag;
//            }

            //Compute components of Bhat evaluated at locations of J_r (rBhat), J_theta (thBhat) and J_phi (phBhat)
            rB_r = 0.25*(B[0][i][j]+B[0][i][j-1]+B[0][i+1][j]+B[0][i+1][j-1]);
            rB_theta = B[1][i][j];
            rB_phi = 0.5*(B[2][i][j]+B[2][i][j-1]);

            thB_r = B[0][i][j];
            thB_theta = 0.25*(B[1][i][j]+B[1][i][j+1]+B[1][i-1][j]+B[1][i-1][j+1]);
            thB_phi = 0.5*(B[2][i][j]+B[2][i-1][j]);

            phB_r = 0.5*(B[0][i][j]+B[0][i][j-1]);
            phB_theta = 0.5*(B[1][i][j]+B[1][i-1][j]);
            phB_phi = 0.25*(B[2][i][j]+B[2][i][j-1]+B[2][i-1][j]+B[2][i-1][j-1]);

            rBmag = sqrt( rB_r*rB_r + rB_theta*rB_theta + rB_phi*rB_phi );
            thBmag = sqrt( thB_r*thB_r + thB_theta*thB_theta + thB_phi*thB_phi );
            phBmag = sqrt( phB_r*phB_r + phB_theta*phB_theta + phB_phi*phB_phi );

            //If Bmag is too small (i.e., field is zero), set components of Bhat equal to zero. Else, define them in the usual way.
            if(rBmag < 1e-10){
                mC.rBhat[0][i][j] = 0.;
                mC.rBhat[1][i][j] = 0.;
                mC.rBhat[2][i][j] = 0.;
            }
            else{
                mC.rBhat[0][i][j] = rB_r/rBmag;
                mC.rBhat[1][i][j] = rB_theta/rBmag;
                mC.rBhat[2][i][j] = rB_phi/rBmag;
            }
            if(thBmag < 1e-10){
                mC.thBhat[0][i][j] = 0.;
                mC.thBhat[1][i][j] = 0.;
                mC.thBhat[2][i][j] = 0.;
            }
            else{
                mC.thBhat[0][i][j] = thB_r/thBmag;
                mC.thBhat[1][i][j] = thB_theta/thBmag;
                mC.thBhat[2][i][j] = thB_phi/thBmag;
            }
            if(phBmag < 1e-10){
                mC.phBhat[0][i][j] = 0.;
                mC.phBhat[1][i][j] = 0.;
                mC.phBhat[2][i][j] = 0.;
            }
            else{
                mC.phBhat[0][i][j] = phB_r/phBmag;
                mC.phBhat[1][i][j] = phB_theta/phBmag;
                mC.phBhat[2][i][j] = phB_phi/phBmag;
            }

            if( (1.-4.*pi*mC.chi[i][j])*(mC.Bmag[i][j]*eta_H[i][j] + eta_O[i][j]) > max_eta_local) max_eta_local = (1.-4.*pi*mC.chi[i][j])*( mC.Bmag[i][j]*eta_H[i][j] + eta_O[i][j] );
        }
    }

    //Fill "corner cells" by mirroring across radial boundary
    mC.Bmag[N_GC-1][N_GC-1] = mC.Bmag[N_GC-1][N_GC];
    mC.Bmag[B.shape()[1]-N_GC-1][N_GC-1] = mC.Bmag[B.shape()[1]-N_GC-1][N_GC];
    mC.Bmag[N_GC-1][B.shape()[2]-N_GC] = mC.Bmag[N_GC-1][B.shape()[2]-N_GC-1];
    mC.Bmag[B.shape()[1]-N_GC-1][B.shape()[2]-N_GC] = mC.Bmag[B.shape()[1]-N_GC-1][B.shape()[2]-N_GC-1];

    exchng2Scalar(mC.Bmag, N_GC, process, dm.stridetype_Sca);
    exchng2Vector(mC.rBhat, N_GC, process, dm.stridetype_Sca);
    exchng2Vector(mC.thBhat, N_GC, process, dm.stridetype_Sca);
    exchng2Vector(mC.phBhat, N_GC, process, dm.stridetype_Sca);

    return;
}
