#include <math.h>
#include <vector>

#include <fftw3.h>
#include <shtns.h>

#include "common.h"
#include "initial_conditions.h"
#include "boundary_conditions.h"

void B_Symmetrize(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm)
{
    MPI_Comm comm1D = process.comm1D; //MPI communicator for decomposed domain
    int world_rank = process.world_rank; //rank of current process
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    std::vector<double> Bphi_BC(B.shape()[2]-2*N_GC); //array to hold the Bphi values to be used to symmetrize outer boundary
    std::vector<double> Bphi_BC_sym(B.shape()[2]-2*N_GC); //array to hold the symmetrized Bphi values

    for(size_t i=N_GC-1; i<=B.shape()[1]-N_GC-1; i++){

        // std::cout << "i = " << i << ", process = " << world_rank << std::endl;

        for(size_t j=0; j<B.shape()[2]-2*N_GC; j++){
            Bphi_BC[j] = B[2][i][N_GC+j];
        }

        Bphi_Outer_Sym(Bphi_BC, Bphi_BC_sym, Ntheta, Ntheta_locs, starts, world_rank, comm1D);

        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){
            B[2][i][j] = Bphi_BC_sym[j-N_GC];
        }

    }

    return;
}


/*
        Sets boundary condition on magnetic field
        Inputs: B: magnetic field as a vector field
                Bparams: BParams object containing information about magnetic field and electric field boundary conditions
                Ntheta: extent of combined domain in partitioned direction
                N_GC: number of ghost cells
                process: Process objects containing information about the simulation domain and the current process
                Ntheta_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
                dm: Domain object containing information about the simulation domain
        Output: B with updated boundary values
*/
void B_BoundaryConditions(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm)
{
    MPI_Comm comm1D = process.comm1D; //MPI communicator for decomposed domain
    int world_rank = process.world_rank; //rank of current process
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t Outeri = B.shape()[1]-N_GC-1;

    //Exchange ghost cells in a periodic manner in the theta-direction
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    std::vector<double> Br_BC(B.shape()[2]-2*N_GC); //array to hold the Br values to be Fourier decomposed and then used to compute Btheta
    std::vector<double> Btheta_BC(B.shape()[2]-2*N_GC); //array to hold the Btheta values as computed from Br according to the boundary condition

    //Fill Br_BC with the values existing at the boundary r=Ro
    for(size_t j=0; j<B.shape()[2]-2*N_GC; j++){
        Br_BC[j] = B[0][B.shape()[1]-N_GC-1][N_GC+j];
    }

    //double n_r = 0.5*(dm.n_r[B.shape()[1]-N_GC-1]+dm.n_r[B.shape()[1]-N_GC]);//0.5*(dm.n_r[B.shape()[1]-N_GC-2]+dm.n_r[B.shape()[1]-N_GC-1]);
    Btheta_BC_Calc(Br_BC, Btheta_BC, Ntheta, Ntheta_locs, starts, dm.n_router, world_rank, comm1D);

    //Symmetrize outermost B_phi values. Helps maintain equatorial symmetry
    // if(dm.InitiallyEquatoriallySymmetric == true){
    //     std::vector<double> Bphi_BC(B.shape()[2]-2*N_GC); //array to hold the Bphi values to be used to symmetrize outer boundary
    //     std::vector<double> Bphi_BC_sym(B.shape()[2]-2*N_GC); //array to hold the symmetrized Bphi values
    //
    //     for(size_t j=0; j<B.shape()[2]-2*N_GC; j++){
    //         Bphi_BC[j] = B[2][B.shape()[1]-N_GC-2][N_GC+j];
    //     }
    //
    //     Bphi_Outer_Sym(Bphi_BC, Bphi_BC_sym, Ntheta, Ntheta_locs, starts, world_rank, comm1D);
    //
    //     for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){
    //         B[2][B.shape()[1]-N_GC-2][j] = Bphi_BC_sym[j-N_GC];
    //     }
    // }

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary r=r_min
            */
            if(bparams.B_perp_lower == "perfect_conducting"){
                B[0][N_GC][j] = 0.; //Br = 0
                B[0][N_GC-1-i][j] = -B[0][N_GC+i+1][j]; //Br = 0
            }
            else if(bparams.B_perp_lower == "continuous"){
                B[0][N_GC-1-i][j] = B[0][N_GC+i+1][j]; //Br is continuous
            }

            if(bparams.B_parallel_lower == "perfect_conducting"){
                B[1][N_GC-1-i][j] = B[1][N_GC+i][j]; //Btheta is continuous
                B[2][N_GC-1-i][j] = B[2][N_GC+i][j]; //Bphi is continuous
            }
            else if(bparams.B_parallel_lower == "zero"){
                B[1][N_GC-1-i][j] = -B[1][N_GC+i][j]; //Btheta = 0
                B[2][N_GC-1-i][j] = -B[2][N_GC+i][j]; //Bphi = 0
            }

            /*
                Upper boundary r=r_max
            */
            if(bparams.B_pol_upper == "vacuum"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j]; //Br is continuous
//                B[0][B.shape()[1]-N_GC+i][j] = (dm.r[B.shape()[1]-N_GC+i]-0.5*dm.Deltar[B.shape()[1]-N_GC+i])/(dm.r[B.shape()[1]-N_GC-2-i]-0.5*dm.Deltar[B.shape()[1]-N_GC-2-i])*B[0][B.shape()[1]-N_GC-2-i][j]; //Br is continuous
                //Btheta equals the value related to the Legendre transform of Br such that B is a potential field in the vacuum exterior.
                B[1][B.shape()[1]-N_GC-1+i][j] = 2.*Btheta_BC[j-N_GC] - B[1][B.shape()[1]-N_GC-2-i][j];
                // B_theta_slope = minmod2( (B[1][Outeri-1][j]-B[1][Outeri-2][j])/(0.5*dm.Deltar[Outeri-1]+0.5*dm.Deltar[Outeri-2]), (B[1][Outeri-2][j]-B[1][Outeri-3][j])/(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-3]) );
                // B[1][B.shape()[1]-N_GC-2][j] = B[1][Outeri-2][j] + B_theta_slope*(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-1]);
                // B[1][B.shape()[1]-N_GC-1+i][j] = ( 2.*dm.r_max*Btheta_BC[j-N_GC] - dm.r[B.shape()[1]-N_GC-2]*B[1][B.shape()[1]-N_GC-2-i][j] )/( dm.r_max + 0.5*dm.Deltar[B.shape()[1]-N_GC-1]  );
            }

            // if(i == 0) {
            //     Dataout << Btheta_BC[j-N_GC] << "  ";
            // }

            if(bparams.B_tor_upper == "vacuum"){
                // B_phi_slope = minmod2( (B[2][Outeri-1][j]-B[2][Outeri-2][j])/(0.5*dm.Deltar[Outeri-1]+0.5*dm.Deltar[Outeri-2]), (B[2][Outeri-2][j]-B[2][Outeri-3][j])/(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-3]) );
                // B[2][B.shape()[1]-N_GC-2][j] = B[2][Outeri-2][j] + B_phi_slope*(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-1]);
                B[2][B.shape()[1]-N_GC-1+i][j] = -B[2][B.shape()[1]-N_GC-2-i][j]; //Bphi = 0
//                B[2][B.shape()[1]-N_GC-1+i][j] = ( -2.*dm.r[B.shape()[1]-N_GC-2]*B[2][B.shape()[1]-N_GC-2-i][j] + 1./3.*dm.r[B.shape()[1]-N_GC-3]*B[2][B.shape()[1]-N_GC-3-i][j] )/( dm.r_max + 0.5*dm.Deltar[B.shape()[1]-N_GC-1]  );
            }

        }
        // Dataout << std::endl;
    }
    // Dataout.close();

    //Exchange ghost cells in a periodic manner in the theta-direction after imposing boundary conditions
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at theta=0 and theta=pi: that is, that Btheta and Bphi vanish here and Br is continuous
    // nbrleft and nbrright = -2 if there is no "neighbouring" process
    //Start at i=N_GC-1 and end at i = B.shape()[i]-N_GC to include "corner" ghost cells
    if( nbrleft < 0 ){
        for(size_t i=N_GC-1; i<B.shape()[1]-N_GC+1; i++){
            for(size_t j=0; j<N_GC; j++){
                B[0][i][N_GC-1-j] = B[0][i][N_GC+j];
                B[1][i][N_GC-1-j] = -B[1][i][N_GC+1+j];
                B[2][i][N_GC-1-j] = -B[2][i][N_GC+j];
            }
            B[1][i][N_GC] = 0;
        }
    }
    if( nbrright < 0 ){
        for(size_t i=N_GC-1; i<B.shape()[1]-N_GC+1; i++){
            for(size_t j=0; j<N_GC; j++){
                B[0][i][B.shape()[2]-N_GC+j] = B[0][i][B.shape()[2]-N_GC-1-j];
                B[2][i][B.shape()[2]-N_GC+j] = -B[2][i][B.shape()[2]-N_GC-1-j];
            }
            B[1][i][B.shape()[2]-N_GC+1] = -B[1][i][B.shape()[2]-N_GC-1];
            B[1][i][B.shape()[2]-N_GC] = 0.;
        }
    }

    return;

}

/*
        Sets boundary condition on magnetic field in intermediate RK steps
        Identical to B_BoundaryConditions except does not update the B_theta values in the ghost cells
        Avoiding imposing new B_theta values in the ghost cells at intermediate steps gives a significant increase in speed, especially when running in parallel,
        with almost identical accuracy
        Inputs: B: magnetic field as a vector field
                Bparams: BParams object containing information about magnetic field and electric field boundary conditions
                Ntheta: extent of combined domain in partitioned direction
                N_GC: number of ghost cells
                process: Process objects containing information about the simulation domain and the current process
                Ntheta_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
                dm: Domain object containing information about the simulation domain
        Output: B with updated boundary values
*/
void B_BoundaryConditionsIntermediate(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm)
{
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t Outeri = B.shape()[1]-N_GC-1;

    //Exchange ghost cells in a periodic manner in the theta-direction
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary r=r_min
            */
            if(bparams.B_perp_lower == "perfect_conducting"){
                B[0][N_GC][j] = 0.; //Br = 0
                B[0][N_GC-1-i][j] = -B[0][N_GC+i+1][j]; //Br = 0
            }
            else if(bparams.B_perp_lower == "continuous"){
                B[0][N_GC-1-i][j] = B[0][N_GC+i+1][j]; //Br is continuous
            }

            if(bparams.B_parallel_lower == "perfect_conducting"){
                B[1][N_GC-1-i][j] = B[1][N_GC+i][j]; //Btheta is continuous
                B[2][N_GC-1-i][j] = B[2][N_GC+i][j]; //Bphi is continuous
            }
            else if(bparams.B_parallel_lower == "zero"){
                B[1][N_GC-1-i][j] = -B[1][N_GC+i][j]; //Btheta = 0
                B[2][N_GC-1-i][j] = -B[2][N_GC+i][j]; //Bphi = 0
            }

            /*
                Upper boundary r=r_max
            */
            if(bparams.B_pol_upper == "vacuum"){
                B[0][B.shape()[1]-N_GC+i][j] = B[0][B.shape()[1]-N_GC-2-i][j]; //Br is continuous
                //Btheta equals the value related to the Fourier coefficients of Br such that B is a potential field in the vacuum exterior.
                // B_theta_slope = minmod2( (B[1][Outeri-1][j]-B[1][Outeri-2][j])/(0.5*dm.Deltar[Outeri-1]+0.5*dm.Deltar[Outeri-2]), (B[1][Outeri-2][j]-B[1][Outeri-3][j])/(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-3]) );
                // B[1][B.shape()[1]-N_GC-2][j] = B[1][Outeri-2][j] + B_theta_slope*(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-1]);
                // B[1][B.shape()[1]-N_GC-1+i][j] = B[1][B.shape()[1]-N_GC-2-i][j];
            }

            if(bparams.B_tor_upper == "vacuum"){
                // B_phi_slope = minmod2( (B[2][Outeri-1][j]-B[2][Outeri-2][j])/(0.5*dm.Deltar[Outeri-1]+0.5*dm.Deltar[Outeri-2]), (B[2][Outeri-2][j]-B[2][Outeri-3][j])/(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-3]) );
                // B[2][B.shape()[1]-N_GC-2][j] = B[2][Outeri-2][j] + B_phi_slope*(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-1]);
                B[2][B.shape()[1]-N_GC-1+i][j] = -B[2][B.shape()[1]-N_GC-2-i][j]; //Bphi = 0
            }

        }
    }

    //Exchange ghost cells in a periodic manner in the theta-direction after imposing boundary conditions
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at theta=0 and theta=pi: that is, that Btheta and Bphi vanish here and Br is continuous
    // nbrleft and nbrright = -2 if there is no "neighbouring" process
    if( nbrleft < 0 ){
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            for(size_t j=0; j<N_GC; j++){
                B[0][i][N_GC-1-j] = B[0][i][N_GC+j];
                B[1][i][N_GC-1-j] = -B[1][i][N_GC+1+j];
                B[2][i][N_GC-1-j] = -B[2][i][N_GC+j];
            }
            B[1][i][N_GC] = 0;
        }
    }
    if( nbrright < 0 ){
        for(size_t i=N_GC; i<B.shape()[1]-N_GC; i++){
            for(size_t j=0; j<N_GC; j++){
                B[0][i][B.shape()[2]-N_GC+j] = B[0][i][B.shape()[2]-N_GC-1-j];
                B[2][i][B.shape()[2]-N_GC+j] = -B[2][i][B.shape()[2]-N_GC-1-j];
            }
            B[1][i][B.shape()[2]-N_GC+1] = -B[1][i][B.shape()[2]-N_GC-1];
            B[1][i][B.shape()[2]-N_GC] = 0.;
        }
    }

    return;

}

/*
    Computes Btheta given Br and the requirement that the external field is a potential field using fast Fourier transforms
    Scatters this information to each process
    Inputs: Br_BC: vector of values of Br on the r=Lr boundary
            Btheta_const: mean of Btheta at upper boundary for each process. Combined to give mean across entire domain.
            Ntheta: length of full simulation domain in theta direction excluding ghost cells
            Ntheta_locs, starts: vectors containing the extent and starting indices of the domain in the decomposed direction
            n_r: sqrt(g_rr) = exp(lambda/2) evaluated at r=Ro
            world_rank: rank of current process
            comm1D: MPI communicator for decomposed domain
    Output: Btheta_BC: vector of values of By on the x=Lx boundary
*/
void Btheta_BC_Calc(std::vector<double> & Br_BC, std::vector<double> & Btheta_BC, size_t Ntheta, std::vector<int> & Ntheta_locs, std::vector<int> & starts, double n_r, int world_rank, MPI_Comm comm1D)
{
    std::vector<double> Br;
    std::vector<double> Btheta;
    if(world_rank == 0){
        Br.resize(Ntheta);
        Btheta.resize(Ntheta);
    }

    // Gather the values of Br(r=r_max) from across the different processes
    // Don't take ghost cells
    MPI_Gatherv(&Br_BC.front(), Ntheta_locs[world_rank], MPI_DOUBLE, &Br.front(), Ntheta_locs.data(), starts.data(), MPI_DOUBLE, 0, comm1D);

    //Gather the values of average of Btheta(r=r_max) from across the different processes
    // Don't take ghost cells
    //static int num_procs = std::size(Ntheta_locs);
    //std::vector<double> Btheta_constVec(num_procs);
    //MPI_Allgather(&Btheta_const, 1, MPI_DOUBLE, Btheta_constVec.data(), 1, MPI_DOUBLE, comm1D);

    /*
        On root process, perform Legendre transform of Br(r=r_max), then manipulate the coefficients
        and invert the Legendre transform to obtain Btheta(r=r_max).
    */
    if(world_rank == 0){

        shtns_cfg shtns; // handle to a sht transform configuration
        std::complex<double> *Br_lm, *Btheta_lm; // spherical harmonics coefficients (l,m space): complex numbers.
        std::complex<double> *Br_temp, *Btheta_temp; //temporary Btheta field. Later will input real parts into vector Btheta

        const int nlat = Br.size(); // number of points in the latitude direction  (constraint: nlat >= 2*lmax+1)
        const int lmax = int(floor(double(nlat/2)))-1; // maximum degree of spherical harmonics
        const int mmax = 1; // maximum order of spherical harmonics
        const int nphi = 2*mmax+1; // number of points in the longitude direction (constraint: nphi >= 2*mmax+1)
        const int mres = 1; // periodicity in phi (1 for full-sphere, 2 for half the sphere, 3 for 1/3, etc...)

        const int m_r = 0; //m-value for Legendre polynomials of B_theta representation.
        const int m_theta = 1; //m-value for Legendre polynomials of B_theta representation.

        shtns_verbose(0); // displays informations during initialization.
        shtns_use_threads(0); // enable multi-threaded transforms (if supported).
        shtns = shtns_init( sht_reg_fast, lmax, mmax, mres, nlat, nphi );

        // allocate temporary Br and Btheta fields
        Br_temp = (std::complex<double> *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(std::complex<double>));
        Btheta_temp = (std::complex<double> *) shtns_malloc( NSPAT_ALLOC(shtns) * sizeof(std::complex<double>));

        // allocate Legendre polynomial representations.
        Br_lm = (std::complex<double> *) shtns_malloc( (lmax+1-m_r) * sizeof(std::complex<double>));
        Btheta_lm = (std::complex<double> *) shtns_malloc( (lmax+1-m_theta) * sizeof(std::complex<double>));

        for(int i=0; i<nlat; i++){
            Br_temp[i] = Br[i];
        }

        //Compute Legendre P^0_l series components of Br
        spat_to_SH_ml(shtns, 0, Br_temp, Br_lm, lmax-m_r);

        //Apply transformation from coefficients of Br Legendre series to Btheta Legendre series
        //Include the sqrt(g_{rr}) and other relativistic factors here.
        double l;
        for(int i=0; i<lmax+1-m_theta; i++){
            l = double(i+1);
            Btheta_lm[i] = -Br_lm[i+1]*sqrt(l/(l+1.))*n_r/( 1. - (2.*l+1.)/(2.*l+2.)*(n_r-1.) );
//            Btheta_lm[i] = -Br_lm[i+1]*sqrt(l/(l+1))*n_r/( 1. - (2.*l+1.)/(2.*l+2.)*(n_r-1.) );
            //Shift indices by -1 because there is no P_{l=0}^{m=1} associated Legendre polynomial
        }
        //transform Btheta from Legendre P^1_l series to spatial representation
        SH_to_spat_ml(shtns, m_theta, Btheta_lm, Btheta_temp, lmax+1-m_theta);

        //Shift Btheta values to left cell edge instead of cell center by averaging cell-centered values on either side of each edge
        //At leftmost cell with Btheta vanishing at theta=0, set Btheta = 0 for j=0.
        for(int j=0; j<nlat; j++){
            if( j == 0 ) Btheta[j] = 0.;
            else Btheta[j] = 0.5*( Btheta_temp[j-1].real() + Btheta_temp[j].real() );
        }

        shtns_free(Br_temp);
        shtns_free(Btheta_temp);
        shtns_free(Br_lm);
        shtns_free(Btheta_lm);
        shtns_destroy(shtns);

    }

    //Send the values of Btheta(r=r_max) back to the processes.
    MPI_Scatterv(&Btheta.front(), Ntheta_locs.data(), starts.data(), MPI_DOUBLE, &Btheta_BC.front(), Ntheta_locs.data()[world_rank], MPI_DOUBLE, 0, comm1D);

    return;
}

/*
 *   Computes Btheta given Br and the requirement that the external field is a potential field using fast Fourier transforms
 *   Scatters this information to each process
 *   Inputs: Bphi_outer: vector of values of Bphi on the r=Lr boundary
 *           Ntheta: length of full simulation domain in theta direction excluding ghost cells
 *           Ntheta_locs, starts: vectors containing the extent and starting indices of the domain in the decomposed direction
 *           world_rank: rank of current process
 *           comm1D: MPI communicator for decomposed domain
 *   Output: Bphi_outer_sym: Bphi: vector of values of Bphi on the x=Lx boundary
 */
void Bphi_Outer_Sym(std::vector<double> & Bphi_outer, std::vector<double> & Bphi_outer_sym, size_t Ntheta, std::vector<int> & Ntheta_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D)
{
    std::vector<double> Bphi;
    std::vector<double> Bphi_sym;
    if(world_rank == 0){
        Bphi.resize(Ntheta);
        Bphi_sym.resize(Ntheta);
    }

    // Gather the values of Br(r=r_max) from across the different processes
    // Don't take ghost cells
    MPI_Gatherv(&Bphi_outer.front(), Ntheta_locs[world_rank], MPI_DOUBLE, &Bphi.front(), Ntheta_locs.data(), starts.data(), MPI_DOUBLE, 0, comm1D);

    /*
     *       On root process, symmetrize B_phi
     */
    if(world_rank == 0){

        for(int j=0; j<Ntheta/2; j++){
            Bphi_sym[j] = 0.5*( Bphi[j] - Bphi[Ntheta-1-j] );
            Bphi_sym[Ntheta-1-j] = -Bphi_sym[j];
        }

    }

    //Send the values of Btheta(r=r_max) back to the processes.
    MPI_Scatterv(&Bphi_sym.front(), Ntheta_locs.data(), starts.data(), MPI_DOUBLE, &Bphi_outer_sym.front(), Ntheta_locs.data()[world_rank], MPI_DOUBLE, 0, comm1D);

    return;
}

/*
        Sets boundary condition on phi component of magnetic field
        Used after Bphi update but before Br and Btheta update
        Inputs: B: magnetic field after Bphi update as a vector field
                Bparams: BParams object containing information about magnetic field and electric field boundary conditions
                Ntheta: extent of combined domain in partitioned direction
                N_GC: number of ghost cells
                process: Process objects containing information about the simulation domain and the current process
                Ntheta_locs and starts: vectors containing the extent and starting indices of the domain in the decomposed direction
                dm: Domain object containing information about the simulation domain
        Output: B with updated boundary values
*/
void Bphi_BoundaryConditions(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm)
{
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    size_t Outeri = B.shape()[1]-N_GC-1;

    //Exchange ghost cells in a periodic manner in the theta-direction
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<B.shape()[2]-N_GC; j++){

            /*
                Lower boundary r=r_min
            */
            if(bparams.B_parallel_lower == "perfect_conducting"){
                B[2][N_GC-1-i][j] = B[2][N_GC+i][j]; //Bphi is continuous
            }
            else if(bparams.B_parallel_lower == "zero"){
                B[2][N_GC-1-i][j] = -B[2][N_GC+i][j]; //Bphi = 0
            }

            /*
                Upper boundary r=r_max
            */
            if(bparams.B_tor_upper == "vacuum"){
                // B_phi_slope = minmod2( (B[2][Outeri-1][j]-B[2][Outeri-2][j])/(0.5*dm.Deltar[Outeri-1]+0.5*dm.Deltar[Outeri-2]), (B[2][Outeri-2][j]-B[2][Outeri-3][j])/(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-3]) );
                // B[2][B.shape()[1]-N_GC-2][j] = B[2][Outeri-2][j] + B_phi_slope*(0.5*dm.Deltar[Outeri-2]+0.5*dm.Deltar[Outeri-1]);
                B[2][B.shape()[1]-N_GC-1+i][j] = -B[2][B.shape()[1]-N_GC-2-i][j]; //Bphi = 0
//                B[2][B.shape()[1]-N_GC-1+i][j] = ( -2.*dm.r[B.shape()[1]-N_GC-2]*B[2][B.shape()[1]-N_GC-2-i][j] + 1./3.*dm.r[B.shape()[1]-N_GC-3]*B[2][B.shape()[1]-N_GC-3-i][j] )/( dm.r_max + 0.5*dm.Deltar[B.shape()[1]-N_GC-1]  );
            }
        }
    }

    //Exchange ghost cells in a periodic manner in the theta-direction after imposing boundary conditions
    exchng2Vector(B, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at theta=0 and theta=pi: that is, that Btheta and Bphi vanish here and Br is continuous
    // nbrleft and nbrright = -2 if there is no "neighbouring" process
    //Start at i=N_GC-1 and end at i = B.shape()[i]-N_GC to include "corner" ghost cells
    if( nbrleft < 0 ){
        for(size_t i=N_GC-1; i<B.shape()[1]-N_GC+1; i++){
            for(size_t j=0; j<N_GC; j++){
                B[2][i][N_GC-1-j] = -B[2][i][N_GC+j];
            }
        }
    }
    if( nbrright < 0 ){
        for(size_t i=N_GC-1; i<B.shape()[1]-N_GC+1; i++){
            for(size_t j=0; j<N_GC; j++){
                B[2][i][B.shape()[2]-N_GC+j] = -B[2][i][B.shape()[2]-N_GC-1-j];
            }
        }
    }

    return;

}

/*
    Sets boundary values of the electric field

    Inputs: phE, cE: redshifted electric field and conjugate electric field in reduced units at cell edges (location of J_phi)
            Bparams: BParams object containing information about magnetic field and electric field boundary conditions
            N_GC: number of ghost cells
            process: Process objects containing information about the simulation domain and the current process
    Output: rE, thE, phE: redshifted electric field with boundary values set
*/
void E_BoundaryConditions(VectorField & phE, VectorField & cE, const BParams & bparams, size_t N_GC, const Process & process, const Domain & dm)
{
    int nbrleft = process.nbrleft; //the rank of the process to the left of the current process
    int nbrright = process.nbrright; //the rank of the process to the right of the current process

    // size_t Outeri = phE.shape()[1]-N_GC-1;
    // double phE_r_slope, phE_th_slope, phE_phi_slope;

    //Impose boundary conditions at non-periodic boundaries
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<phE.shape()[2]-N_GC; j++){
            /*
                Lower boundary r=r_min
            */
            if(bparams.E_perp_lower == "perfect_conducting"){
                phE[0][N_GC-1-i][j] = phE[0][N_GC+i][j]; //Er is continuous
                cE[0][N_GC-1-i][j] = cE[0][N_GC+i][j]; //Er is continuous
            }
            if(bparams.E_parallel_lower == "perfect_conducting"){
                // cE[1][N_GC-1-i][j] = -cE[1][N_GC+i][j]; //Etheta = 0
                // cE[2][N_GC-1-i][j] = -cE[2][N_GC+i][j]; //Ephi = 0
                phE[1][N_GC-1-i][j] = -phE[1][N_GC+1+i][j]; //Etheta = 0
                phE[2][N_GC-1-i][j] = -phE[2][N_GC+1+i][j]; //Ephi = 0
                phE[1][N_GC][j] = 0.;
                phE[2][N_GC][j] = 0.;
                cE[1][N_GC-1-i][j] = -cE[1][N_GC+1+i][j]; //Etheta = 0
                cE[2][N_GC-1-i][j] = -cE[2][N_GC+1+i][j]; //Ephi = 0
                cE[1][N_GC][j] = 0.;
                cE[2][N_GC][j] = 0.;
            }
            else if(bparams.E_parallel_lower == "continuous"){
                // thE[2][N_GC-1-i][j] = thE[2][N_GC+1+i][j]; //Ephi is continuous
                phE[2][N_GC-1-i][j] = phE[2][N_GC+1+i][j]; //Ephi is continuous
                cE[2][N_GC-1-i][j] = cE[2][N_GC+1+i][j]; //Ephi is continuous
            }
            /*
                Upper boundary r=r_max
            */
            // phE at outer boundary i is constructed by linear interpolation from phE at i-1 and i-2
            // phE[0][Outeri][j] = 0.5*(phE[0][Outeri][j]+phE[0][Outeri-1][j]); //
            // phE[1][Outeri][j] = 0.5*(phE[1][Outeri][j]+phE[1][Outeri-1][j]); //
            // phE[2][Outeri][j] = 0.5*(phE[2][Outeri][j]+phE[2][Outeri-1][j]); //
            // phE_r_slope = minmod2( (phE[0][Outeri][j]-phE[0][Outeri-1][j])/dm.Deltar[Outeri-1], (phE[0][Outeri-1][j]-phE[0][Outeri-2][j])/dm.Deltar[Outeri-2] );
            // phE_th_slope = minmod2( (phE[1][Outeri][j]-phE[1][Outeri-1][j])/dm.Deltar[Outeri-1], (phE[1][Outeri-1][j]-phE[1][Outeri-2][j])/dm.Deltar[Outeri-2] );
            // phE_phi_slope = minmod2( (phE[2][Outeri][j]-phE[2][Outeri-1][j])/dm.Deltar[Outeri-1], (phE[2][Outeri-1][j]-phE[2][Outeri-2][j])/dm.Deltar[Outeri-2] );
            // phE[0][Outeri][j] = phE[0][Outeri-1][j] + phE_r_slope*dm.Deltar[Outeri-1]; //
            // phE[1][Outeri][j] = phE[1][Outeri-1][j] + phE_th_slope*dm.Deltar[Outeri-1]; //
            // phE[2][Outeri][j] = phE[2][Outeri-1][j] + phE_phi_slope*dm.Deltar[Outeri-1]; //
            // phE[0][Outeri][j] = phE[0][Outeri-1][j] + (phE[0][Outeri-1][j]-phE[0][Outeri-2][j])/dm.Deltar[Outeri-2]*( dm.r[Outeri-1]-dm.r[Outeri-2] + 0.5*(dm.Deltar[Outeri-1]-dm.Deltar[Outeri-2]) ); //
            // phE[1][Outeri][j] = phE[1][Outeri-1][j] + (phE[1][Outeri-1][j]-phE[1][Outeri-2][j])/dm.Deltar[Outeri-2]*( dm.r[Outeri-1]-dm.r[Outeri-2] + 0.5*(dm.Deltar[Outeri-1]-dm.Deltar[Outeri-2]) ); //
            // phE[2][Outeri][j] = phE[2][Outeri-1][j] + (phE[2][Outeri-1][j]-phE[2][Outeri-2][j])/dm.Deltar[Outeri-2]*( dm.r[Outeri-1]-dm.r[Outeri-2] + 0.5*(dm.Deltar[Outeri-1]-dm.Deltar[Outeri-2]) ); //

//            if(bparams.E_perp_upper == "vacuum"){
//                phE[0][phE.shape()[1]-N_GC-1][j] = 0; //Er vanishes
//            }
//            else if(bparams.E_perp_upper == "continuous"){
//                phE[0][phE.shape()[1]-N_GC+i][j] = phE[0][phE.shape()[1]-N_GC-2-i][j]; //Er is continuous
//            }
//            if(bparams.E_parallel_upper == "vacuum"){
//                phE[1][phE.shape()[1]-N_GC+i][j] = -phE[1][phE.shape()[1]-N_GC-2-i][j]; //Etheta = 0
//                phE[2][phE.shape()[1]-N_GC+i][j] = -phE[2][phE.shape()[1]-N_GC-2-i][j]; //Ephi = 0
////                phE[1][phE.shape()[1]-N_GC-1][j] = 0.; //Etheta = 0
////                phE[2][phE.shape()[1]-N_GC-1][j] = 0.; //Ephi = 0
//                //E[1][E.shape()[1]-N_GC-1][j] = 0.; //E_theta does not need to vanish at outer boundary
//                //E[2][E.shape()[1]-N_GC-1][j] = 0.; //E_phi does not need to vanish at outer boundary
//            }
//            else if(bparams.E_parallel_upper == "continuous"){
//                phE[1][phE.shape()[1]-N_GC+i][j] = phE[1][phE.shape()[1]-N_GC-2-i][j]; //Etheta is continuous
//                phE[2][phE.shape()[1]-N_GC+i][j] = phE[2][phE.shape()[1]-N_GC-2-i][j]; //Ephi is continuous
//            }
        }
    }

    //Exchange ghost cells in a periodic manner in the y-direction
    exchng2Vector(cE, N_GC, process, dm.stridetype_Vec);
    exchng2Vector(phE, N_GC, process, dm.stridetype_Vec);

    //Impose boundary conditions at theta=0 and theta=pi: that is, that Etheta and Ephi vanish here and Er is continuous
    // nbrleft and nbrright = -2 if there is no "neighbouring" process
    if( nbrleft < 0 ){
        for(size_t i=N_GC; i<phE.shape()[1]-N_GC; i++){
            for(size_t j=0; j<N_GC; j++){
                // cE[0][i][N_GC-1-j] = cE[0][i][N_GC+j];
                // cE[1][i][N_GC-1-j] = -cE[1][i][N_GC+j];
                // cE[2][i][N_GC-1-j] = -cE[2][i][N_GC+j];
                phE[0][i][N_GC-1-j] = phE[0][i][N_GC+1+j];
                phE[1][i][N_GC-1-j] = -phE[1][i][N_GC+1+j];
                phE[2][i][N_GC-1-j] = -phE[2][i][N_GC+1+j];
                cE[0][i][N_GC-1-j] = cE[0][i][N_GC+1+j];
                cE[1][i][N_GC-1-j] = -cE[1][i][N_GC+1+j];
                cE[2][i][N_GC-1-j] = -cE[2][i][N_GC+1+j];
            }
            phE[1][i][N_GC] = 0.;
            phE[2][i][N_GC] = 0.;
            cE[1][i][N_GC] = 0.;
            cE[2][i][N_GC] = 0.;
        }
    }
    if( nbrright < 0 ){
        for(size_t i=N_GC; i<phE.shape()[1]-N_GC; i++){
            // cE[0][i][cE.shape()[2]-N_GC] = cE[0][i][cE.shape()[2]-N_GC-1];
            // cE[1][i][cE.shape()[2]-N_GC] = -cE[1][i][cE.shape()[2]-N_GC-1];
            // cE[2][i][cE.shape()[2]-N_GC] = -cE[2][i][cE.shape()[2]-N_GC-1];
            phE[0][i][phE.shape()[2]-N_GC+1] = phE[0][i][phE.shape()[2]-N_GC-1];
            phE[1][i][phE.shape()[2]-N_GC+1] = -phE[1][i][phE.shape()[2]-N_GC-1];
            phE[2][i][phE.shape()[2]-N_GC+1] = -phE[2][i][phE.shape()[2]-N_GC-1];
            cE[0][i][phE.shape()[2]-N_GC+1] = cE[0][i][phE.shape()[2]-N_GC-1];
            cE[1][i][phE.shape()[2]-N_GC+1] = -cE[1][i][phE.shape()[2]-N_GC-1];
            cE[2][i][phE.shape()[2]-N_GC+1] = -cE[2][i][phE.shape()[2]-N_GC-1];
            //phE[0][i][phE.shape()[2]-N_GC] = phE[0][i][phE.shape()[2]-N_GC-1];
            // cE[1][i][cE.shape()[2]-N_GC+1] = -cE[1][i][cE.shape()[2]-N_GC-2];
            // cE[2][i][cE.shape()[2]-N_GC+1] = -cE[2][i][cE.shape()[2]-N_GC-2];
            phE[1][i][phE.shape()[2]-N_GC] = 0.;
            phE[2][i][phE.shape()[2]-N_GC] = 0.;
            cE[1][i][phE.shape()[2]-N_GC] = 0.;
            cE[2][i][phE.shape()[2]-N_GC] = 0.;
        }
    }

    return;

}

/*
        Sets boundary condition on heat flux density using temperature
        Inputs: T, q: redshifted temperature and heat flux density times lapse function squared in reduced units
                Tcore: core temperature in reduced units
                tC: TransCoeffs object containing thermal conductivity kappa
                mC: MagCoeffs object containing properties of magnetic field and magnetization
                tparams: TParams object containing information about initial temperature and thermal boundary conditions
                dm: Domain object containing information about the simulation domain
                N_GC: number of ghost cells. Note: this function is written assuming N_GC = 2
                process: Process objects containing information about the simulation domain and the current process
        Output: T with updated boundary values
*/
void T_BoundaryConditions(ScalarField & T, VectorField & q, double Tcore, TransCoeffs & tC, MagCoeffs & mC, TParams & tparams, const Domain & dm, size_t N_GC, const Process & process)
{
    static double sigma_SBbar = sigma_SB*t_0*pow(T_0,3.)/(s_0*L_0); //dimensionless Stefan-Boltzmann constant
    static double g_s14 = tparams.g_s14; //gravitational acceleration at stellar surface in units of 1e14 cm/s^2
    static double r_s = tparams.r_s; //outer (surface) radius of star in reduced units
    static double r_b = dm.r[T.shape()[0]-N_GC-1]+0.5*dm.Deltar[T.shape()[0]-N_GC-1]; //inner radius of heat-blanketing envelope = outer radius of simulation domain. Note that we add Deltar to r to get this because of how r[T.shape()[0]-N_GC-1] is defined
    double Ts, dTsdTb; //local surface temperature and its derivative w.r.t. redshifted envelope bottom temperature in reduced units
    double Ts_star; //local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s
    double Tb8, Tb9; //local temperature at bottom of heat blanketing envelope in units of 10^8 K or 10^9 K
    double theta_mc, theta_B; //magnetic colatitude in radians and magnetic field inclination angle in radians

    double a, b, RadialTerm, TangentialTerm, kappa_r, kappa_th;

    size_t O_index = T.shape()[0]-N_GC-1; //radial index of exterior cell (whose value we are setting)
    double DeltarOuter = 0.5*( dm.Deltar[O_index-1] + dm.Deltar[O_index] );
    double B12 = tparams.B_pol*B_0/1e12; //initial magnetic field at poles in 10^12 G

    static double g_f = dm.n_t[O_index-1]*dm.n_t[O_index-1]*pow(r_s/r_b,2.); //factor of exp(nu/2)*(R_out/R_b)^2
    static double x_g = (1.-1./(dm.n_r[O_index-1]*dm.n_r[O_index-1]))*dm.r[O_index-1]/r_s;
    // static double a_g; //general relativistic correction parameter for magnetic field inclination angle. Equals 1/2 in nonrelativistic limit x_g = 0. Only used for T_upper = "magnetic_P01"
    // if( dm.n_r[O_index-1] <= 1+1e-10 ) a_g = 0.5;
    // else a_g = -( (1.-x_g)*log(1.-x_g)+x_g-0.5*x_g*x_g )/( (log(1.-x_g)+x_g+0.5*x_g*x_g)*sqrt(1.-x_g) );

    double Ts4_sum = 0.; //running sum of (local surface temperature)^4 values in reduced units.

    exchng2Scalar(T, N_GC, process, dm.stridetype_Sca); //needed to compute TangentialTerm correctly when running on multiple processes

    //Impose boundary conditions at non-periodic boundaries
    for(size_t j=N_GC; j<T.shape()[1]-N_GC; j++){
        /*
            Lower boundary r=r_min
        */
        if(tparams.T_lower == "fake_core"){
            T[N_GC-1][j] = 2.*Tcore - T[N_GC][j]; //set core temperature based on fake core temperature evolution
            T[N_GC-2][j] = T[N_GC-1][j]; //dT/dr = 0 across fake core of width N ghost cells
        }
        /*
            Upper boundary r=r_max
        */
        if(tparams.T_upper == "non-magnetic_G83"){
            Tb8 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e8; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^8 K
            T_sGudmundsson1983(Tb8, g_s14, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Gudmundsson, Pethick and Epstein, ApJ 272, 286 (1983)
        }
        else if(tparams.T_upper == "non-magnetic_P01"){
            Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^8 K
            T_sPotekhin2001(Tb9, g_s14, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
        }
        else if(tparams.T_upper == "magnetic_P01"){
            Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
            // theta_B = atan(a_g*tan(dm.theta[j])); //magnetic field inclination angle in radians. Equals the local angle between the magnetic field and the radial direction. Includes general relativistic correction
            theta_B = acos(mC.phBhat[0][O_index][j]); //magnetic field inclination angle in radians
            T_sPotekhinMag2001(Tb9, B12, g_s14, theta_B, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
        }
        else if(tparams.T_upper == "magnetic_P15"){
            Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
            theta_mc = dm.theta[j];  //magnetic colatitude in radians
            T_sPotekhin2015(Tb9, B12, g_s14, theta_mc, Ts, Ts_star, dTsdTb);// T_s - T_b relation for heat blanketing envelope is from Potekhin, Pons & Page, Space Sci. Rev. 191, 239 (2015), Appendix B.
        }
        else {
            std::cout << "Outer thermal boundary condition not set" << std::endl;
        }
        a = g_f*sigma_SBbar*Ts_star*Ts_star*Ts_star*( Ts_star - 2.*dTsdTb/dm.n_t[O_index-1]*(T[O_index][j] + T[O_index-1][j]) );
        b = 2.*g_f*sigma_SBbar*Ts_star*Ts_star*Ts_star*dTsdTb/dm.n_t[O_index-1];

        if(tparams.conductivity_anisotropy == true){
            kappa_r = 0.5*( (tC.kappa_delta[O_index-1][j]+tC.kappa_delta[O_index][j])*mC.thBhat[0][O_index][j]*mC.thBhat[0][O_index][j] + (tC.kappa_perp[O_index-1][j]+tC.kappa_perp[O_index][j]) );
            kappa_th = 0.5*( (tC.kappa_delta[O_index-1][j]+tC.kappa_delta[O_index][j])*mC.thBhat[0][O_index][j]*mC.thBhat[1][O_index][j] - (tC.kappa_H[O_index-1][j]+tC.kappa_H[O_index][j])*mC.thBhat[2][O_index][j] );
            RadialTerm = dm.n_t[O_index-1]*kappa_r/(DeltarOuter*dm.n_r[O_index-1]); //K_r
            TangentialTerm = 0*dm.n_t[O_index-1]*kappa_th/r_b*0.5/tan(dm.Deltatheta/2.)*( T[O_index-1][j+1] - T[O_index-1][j-1] )/2.; //K_theta
            // std::cout << "j = " << j << ", thBhat_r = " << mC.thBhat[0][O_index][j] << std::endl;
        }
        else{
            kappa_r = 0.5*(tC.kappa[O_index-1][j]+tC.kappa[O_index][j]); //if only one isotropic thermal conductivity, do not need to take an average
            RadialTerm = dm.n_t[O_index-1]*kappa_r/(DeltarOuter*dm.n_r[O_index-1]); //K_r
            TangentialTerm = 0.;
        }

        q[0][O_index][j] = g_f*sigma_SBbar*Ts_star*Ts_star*Ts_star*Ts_star;
        T[O_index][j] = T[O_index-1][j]; //( RadialTerm - b )/( RadialTerm + b )*T[O_index-1][j] - ( TangentialTerm + a )/( RadialTerm + b );
        T[O_index+1][j] = (T[O_index][j]-T[O_index-1][j])/(dm.n_r[O_index-1]*DeltarOuter)*( dm.n_r[O_index-1]*DeltarOuter + dm.n_r[O_index]*0.5*(dm.Deltar[O_index] + dm.Deltar[O_index+1]) ) + T[O_index-1][j];

        Ts4_sum += 2.*pi*sin(dm.theta[j])*2.*sin(dm.Deltatheta/2.)*Ts*Ts*Ts*Ts; //local surface temperature^4 integrated over area element.

    }

    tparams.Ts4_av = Ts4_sum;

    //Exchange ghost cells in a periodic manner in the y-direction
    exchng2Scalar(T, N_GC, process, dm.stridetype_Sca);

    //Impose boundary conditions at theta=0 and theta=pi: that is, that the heat flux in the angular direction dT/dtheta=0 vanishes
    // nbrleft and nbrright = -2 if there is no "neighbouring" process
    if( process.nbrleft < 0 ){
        for(size_t i=0; i<T.shape()[0]; i++){
            for(size_t j=0; j<N_GC; j++){
                T[i][N_GC-1-j] = T[i][N_GC+j];
            }
        }
    }
    if( process.nbrright < 0 ){
        for(size_t i=0; i<T.shape()[0]; i++){
            for(size_t j=0; j<N_GC; j++){
                T[i][T.shape()[1]-N_GC+j] = T[i][T.shape()[1]-N_GC-1-j];
            }
        }
    }

    return;
}

/*
        Non-magnetic heat-blanketing envelope model from Gudmundsson, Pethick and Epstein, ApJ 272, 286 (1983)
        Calculates both T_s(T_b) and dT_s/dT_b(T_b).
        Inputs: Tb8: local temperature of inner edge of heat-blanketing envelope in units of 10^8 K
                g_s14: surface gravitational field strength in units of 10^14 cm/s^2
        Outputs: Ts: local surface temperature in reduced units
                 Ts_star: local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s.
                        Equals Ts for this envelope model.
                 dTsdTb: derivative of Ts with respect to Tb (dimensionless)
*/
void T_sGudmundsson1983(double Tb8, double g_s14, double & Ts, double & Ts_star, double & dTsdTb){

    static double a_GPE = 1./1.82; //exponent of T_b8 in Gudmundsson et al. 1983 model
    static double g_s14pp455 = pow( g_s14, 0.455 ); //gs,14^0.455

    Ts = 1.e6*pow( g_s14pp455*Tb8/1.288,a_GPE )/T_0; //surface temperature in reduced units
    Ts_star = Ts;
    dTsdTb = a_GPE*Ts*T_0/(1.e8*Tb8); //derivative of surface temperature w.r.t T_b. Dimensionless

    return;
}

/*
        Non-magnetic heat-blanketing envelope model from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
        Calculates both T_s(T_b) and dT_s/dT_b(T_b). Generally similar results to using T_sGudmundsson1983.
        Inputs: Tb9: local temperature of inner edge of heat-blanketing envelope in units of 10^9 K
                g_s14: surface gravitational field strength in units of 10^14 cm/s^2
        Outputs: Ts: local surface temperature in reduced units
                 Ts_star: local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s.
                        Equals Ts for this envelope model.
                 dTsdTb: derivative of Ts with respect to Tb (dimensionless)
*/
void T_sPotekhin2001(double Tb9, double g_s14, double & Ts, double & Ts_star, double & dTsdTb){

    double zeta = Tb9 - 0.001*pow(g_s14,0.25)*sqrt( 7.*Tb9 ); //effective temperature at inner edge of heat-blanketing envelope in units of 10^9 K. Dimensionless
    Ts = 1.e6*pow( g_s14*( pow(7.*zeta,2.25) + pow(zeta/3.,1.25) ),0.25 )/T_0; //local surface temperature in reduced units
    Ts_star = Ts;
    double dzetadTb9 = 1. - 0.001*pow(g_s14,0.25)*sqrt( 7./Tb9 ); //derivative of zeta w.r.t. Tb9. Dimensionless.
    dTsdTb = 0.25*g_s14*(1e15/(Ts*T_0*Ts*T_0*Ts*T_0))*( 2.25*pow(7.*zeta,2.25) + 1.25*pow(zeta/3.,1.25) )/zeta*dzetadTb9; //derivative of surface temperature w.r.t T_b. dimensionless

    return;
}

/*
        Magnetic heat-blanketing envelope model from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
        Uses modified coefficients taken from Pons, Miralles and Geppert A&A 496, 207–216 (2009)
        Calculates both T_s(T_b) and dT_s/dT_b(T_b).
        Inputs: Tb9: local temperature of inner edge of heat-blanketing envelope in units of 10^9 K
                B12: magnetic field strength at magnetic pole in units of 10^12 G
                g_s14: surface gravitational field strength in units of 10^14 cm/s^2
                theta_B: magnetic field inclination angle with general relativistic correction (radians)
        Outputs: Ts: local surface temperature in reduced units
                 Ts_star: local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s.
                        Equals Ts for this envelope model.
                 dTsdTb: derivative of Ts with respect to Tb (dimensionless)
*/
void T_sPotekhinMag2001(double Tb9, double B12, double g_s14, double theta_B, double & Ts, double & Ts_star, double & dTsdTb){

    double zeta = Tb9 - 0.001*pow(g_s14,0.25)*sqrt( 7.*Tb9 ); //effective temperature at inner edge of heat-blanketing envelope in units of 10^9 K. Dimensionless
    double chi_par = 1. + 0.05*pow( B12,0.25 )/pow( Tb9,0.240 );
    double chi_perp = sqrt( 1. + 0.07*B12/pow( 0.03+Tb9,0.559 ) )/pow( 1. + 0.9*B12/(0.03+Tb9),0.4 );
    double chi_92 = pow( chi_par,9./2. )*cos(theta_B)*cos(theta_B) + pow( chi_perp,9./2. )*sin(theta_B)*sin(theta_B);
    double chi = pow( chi_92,2./9. );
    double Ts0 = 1.e6*pow( g_s14*( pow(7.*zeta,2.25) + pow(zeta/3.,1.25) ),0.25 )/T_0; //local surface temperature without magnetic field effect in reduced units
    Ts = chi*Ts0; //local surface temperature with magnetic field effect in reduced units
    Ts_star = Ts;
    double dzetadTb9 = 1. - 0.001*pow(g_s14,0.25)*sqrt( 7./Tb9 ); //derivative of zeta w.r.t. Tb9. Dimensionless.
    double dTs0dTb = 0.25*g_s14*(1e15/(Ts0*T_0*Ts0*T_0*Ts0*T_0))*( 2.25*pow(7.*zeta,2.25) + 1.25*pow(zeta/3.,1.25) )/zeta*dzetadTb9; //derivative of surface temperature w.r.t T_b. dimensionless
    double dchi_pardTb = -0.012*pow( B12,0.25 )/pow( Tb9,1.240 );
    double dchi_perpdTb = -0.019565*B12/pow( 0.03+Tb9,1.559 )/sqrt( 1. + 0.07*B12/pow( 0.03+Tb9,0.559 ) )/pow( 1. + 0.9*B12/(0.03+Tb9),0.4 )
                    + 0.36*B12/pow( 0.03+Tb9,2. )*sqrt( 1. + 0.07*B12/pow( 0.03+Tb9,0.559 ) )/pow( 1. + 0.9*B12/(0.03+Tb9),1.4 );
    double dchidTb = T_0/1e9*chi/chi_92*( dchi_pardTb *pow( chi_par,7./2. )*cos(theta_B)*cos(theta_B) + dchi_perpdTb*pow( chi_perp,7./2. )*sin(theta_B)*sin(theta_B) );
    dTsdTb = dTs0dTb*chi + Ts*dchidTb;

    return;
}

/*
        Magnetic heat-blanketing envelope model from Potekhin, Pons & Page, Space Sci. Rev. 191, 239 (2015), Appendix B. Includes envelope neutrino emission correction.
        Calculates both T_s(T_b) and dT_s/dT_b(T_b) assuming magnetic dipole field threading envelope.
        Inputs: Tb9: local temperature of inner edge of heat-blanketing envelope in units of 10^9 K
                B12: magnetic field strength at magnetic pole in units of 10^12 G
                g_s14: surface gravitational field strength in units of 10^14 cm/s^2
                theta: magnetic colatitude (radians). For dipole field aligned with polar axis, this is just the polar angle coordinate theta for theta<pi/2 and pi-theta for theta>pi/2.
        Outputs: Ts: local surface temperature in reduced units
                 Ts_star: local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s
                 dTsdTb: derivative of Ts with respect to Tb (dimensionless)
*/
void T_sPotekhin2015(double Tb9, double B12, double g_s14, double theta, double & Ts, double & Ts_star, double & dTsdTb){

    static double g_s14pp65 = pow( g_s14,0.65 ); //gs,14^0.65
    double a, T0, T1, Tp0, Tpmax, Tp, denom, f1, f2, f, Teq, Teq0, a2, a1, a_factor_num, a_factor_denom;
    double dT1dTb9, dT0dTb9, da1dTb, da2dTb, dTp0dTb9, dTpdTb, dfTb9, dTeqdTb;

    a = 0.337/( 1. + 0.02*sqrt(B12) ); //dimensionless
    T0 = pow( 15.7*Tb9*sqrt(Tb9) + 1.36*Tb9,0.3796 ); //dimensionless. Not to be confused with T_0, the characteristic temperature!
    T1 = 1.13*pow( B12,0.119 )*pow( Tb9,a ); //dimensionless
    Tp0 = 1e6*pow( g_s14*( T1*T1*T1*T1 + ( 1. + 0.15*sqrt(B12) )*T0*T0*T0*T0 ),0.25 ); //in K
    Tpmax = 1e6*( 5.2*g_s14pp65 + 0.093*sqrt(g_s14*B12) ); // in K
    Tp = Tp0/pow( 1. + Tp0*Tp0*Tp0*Tp0/(Tpmax*Tpmax*Tpmax*Tpmax),0.25 ); //in K
    denom = ( B12 + 450.*Tb9 + 119.*B12*Tb9 )*( B12 + 450.*Tb9 + 119.*B12*Tb9 )*( B12 + 450.*Tb9 + 119.*B12*Tb9 )*( B12 + 450.*Tb9 + 119.*B12*Tb9 ); //dimensionless
    f1 = pow( 1230.*Tb9,3.35 )*B12*sqrt(1. + 2.*B12*B12)/denom; //dimensionless
    f2 = 0.0066*B12*B12*sqrt(B12)/( sqrt(Tb9) + 0.00258*B12*B12*sqrt(B12) ); //dimensionless
    f = 1. + f1 + f2; //dimensionless
    Teq = Tp/f; //in K
    Teq0 = Tp0/f; //in K
    a2 = 10.*B12/( sqrt(Tb9) + 0.1*B12*pow( Tb9,-0.25 ) ); //dimensionless
    a1 = a2*sqrt(Tb9)/3.; //dimensionless
    a_factor_num = ( 1. + a1 + a2 )*cos(theta)*cos(theta); //dimensionless
    a_factor_denom = ( 1. + a1*abs(cos(theta)) + a2*cos(theta)*cos(theta) ); //dimensionless
    Ts = ( Teq + a_factor_num/a_factor_denom*(Tp - Teq) )/T_0; //local surface temperature in reduced units
    Ts_star = ( Teq0 + a_factor_num/a_factor_denom*(Tp0 - Teq0) )/T_0;//local surface temperature without corrections due to neutrino emission in crust: the temperature corresponding to the actual blackbody emission required at r=r_b to produce a surface temperature Ts at r=r_s

    dT1dTb9 = a*T1/Tb9; //dimensionless
    dT0dTb9 = 0.3796*T0*( 23.55*sqrt(Tb9) + 1.36 )/( 15.7*Tb9*sqrt(Tb9) + 1.36*Tb9 ); //dimensionless
    da2dTb = -a2/1e9*( 0.5/sqrt(Tb9) - 0.025*B12*pow(Tb9,-1.25) )/( sqrt(Tb9) + 0.1*B12*pow(Tb9,-0.25) ); //in K^-1
    da1dTb = a1/(2.*Tb9*1e9) + da2dTb*sqrt(Tb9)/3.; //in K^-1
    dTp0dTb9 = 1e24*g_s14/(Tp0*Tp0*Tp0)*( T1*T1*T1*dT1dTb9 + ( 1. + 0.15*sqrt(B12) )*T0*T0*T0*dT0dTb9 ); //in K
    dTpdTb = dTp0dTb9/1e9*Tp/Tp0*( 1. - Tp*Tp*Tp*Tp/(Tpmax*Tpmax*Tpmax*Tpmax) ); //dimensionless
    dfTb9 = f1*3.35/Tb9 - 4.*f1*(450.+119.*B12)/( B12 + 450.*Tb9 + 119.*B12*Tb9 ) - 0.5*f2/( sqrt(Tb9) + 0.00258*B12*B12*sqrt(B12) )/sqrt(Tb9); //dimensionless
    dTeqdTb = Teq/Tp*( dTpdTb - Teq*dfTb9/1e9 ); //dimensionless
    dTsdTb = dTeqdTb + a_factor_num/a_factor_denom*(dTpdTb - dTeqdTb) + ( da1dTb + da2dTb )*cos(theta)*cos(theta)/a_factor_denom*(Tp - Teq)
            - a_factor_num/(a_factor_denom*a_factor_denom)*( da1dTb*abs(cos(theta)) + da2dTb*cos(theta)*cos(theta) )*(Tp - Teq); //dimensionless

    return;
}

/*
        Sets boundary condition on temperature in inner radial boundary. Used in IMEX timestepping methods.
        Inputs: Tcore: current core redshifted temperature in reduced units
                tC: TransCoeffs object containing thermal conductivity kappa
                N_GC: number of ghost cells
                dm, process: Domain and Process objects containing information about the simulation domain and the current process
                corethermalparams: CoreThermalParams object containing information about thermodynamic properties of the core
                a_ii: diagonal entry in Butcher tableau for this implicit RK step
        Output: b_0_term, d_0_term: contributions to the tridiagonal matrix b_0 and column vector d_0
*/
void T_InnerBC_r(double & b_0_term, double & d_0_term, double Tcore, const Domain & dm, const CoreThermalParams & corethermalparams, double a_ii)
{
    double C_v = gsl_spline_eval( corethermalparams.C_vCore, Tcore, corethermalparams.C_vCore_acc );
    double Q_nu = gsl_spline_eval( corethermalparams.Q_nuCore, Tcore, corethermalparams.Q_nuCore_acc );

    double DeltaT = 1e-4*Tcore; //temperature step used in computing dQ_nuCoredT, the finite difference derivative of Q_nuCore with respect to temperature
    double dQ_nuCoredT = ( gsl_spline_eval( corethermalparams.Q_nuCore, Tcore+0.5*DeltaT, corethermalparams.Q_nuCore_acc ) - gsl_spline_eval( corethermalparams.Q_nuCore, Tcore-0.5*DeltaT, corethermalparams.Q_nuCore_acc ) )/DeltaT;
    // double denom = 1. - dm.Deltat*a_ii*( dQ_nuCoredT/C_v + corethermalparams.HeatFlux_pref/C_v ); //denominator of implicit calculation of T_core
    double denom = 1. - 0*dm.Deltat*a_ii*( dQ_nuCoredT/C_v ); //denominator of implicit calculation of T_core

    //T_{i=N_GC-1} = b_0_term*T_{i=N_GC} + d_0_term
    // b_0_term = -dm.Deltat*a_ii*corethermalparams.HeatFlux_pref/C_v/denom;
    // d_0_term = ( Tcore + dm.Deltat*a_ii*( Q_nu/C_v + corethermalparams.HeatFlux_H/C_v - dQ_nuCoredT*Tcore/C_v ) )/denom;
    b_0_term = 0.;//-dm.Deltat*a_ii*corethermalparams.HeatFlux_pref/C_v/denom;
    d_0_term = ( Tcore + 0*dm.Deltat*a_ii*( Q_nu/C_v - dQ_nuCoredT*Tcore/C_v - corethermalparams.HeatFlux/C_v ) )/denom;

    return;
}

/*
        Sets boundary condition on temperature at outer radial boundary. Used in IMEX timestepping methods.
        Inputs: j: Index in theta direction
                T: temperature as a scalar field
                Tcore: core temperature in reduced units
                tC: TransCoeffs object containing thermal conductivity kappa
                mC: MagCoeffs object containing properties of magnetic field and magnetization
                tparams: TParams object containing information about initial temperature and thermal boundary conditions
                dm: Domain object containing information about the simulation domain
                N_GC: number of ghost cells. Note: this function is written assuming N_GC = 2
        Output: b_n_term, d_n_term: contributions to the tridiagonal matrix b_n and column vector d_n
*/
void T_OuterBC_r(size_t j, double & b_n_term, double & d_n_term, ScalarField & T, TransCoeffs & tC, MagCoeffs & mC, const TParams & tparams, const Domain & dm, size_t N_GC)
{
    static double sigma_SBbar = sigma_SB*t_0*pow(T_0,3.)/(s_0*L_0); //dimensionless Stefan-Boltzmann constant
    static double g_s14 = tparams.g_s14; //gravitational acceleration at stellar surface in units of 1e14 cm/s^2
    static double r_s = tparams.r_s; //outer (surface) radius of star in reduced units
    static double r_b = dm.r[T.shape()[0]-N_GC-1]-0.5*dm.Deltar[T.shape()[0]-N_GC-1]; //inner radius of heat-blanketing envelope = outer radius of simulation domain
    double Ts, Ts_star= 0., dTsdTb = 0.; //local surface temperature and its derivative w.r.t. redshifted envelope bottom temperature in reduced units
    double kappa_r, kappa_th; //surface-averaged radial and tangential thermal conductivities in reduced units
    double Tb8, Tb9; //local temperature at bottom of heat blanketing envelope in units of 10^8 K or 10^9 K
    double theta_mc, theta_B; //magnetic colatitude in radians and magnetic field inclination angle in radians

    double RadialTerm, TangentialTerm, a = 0., b = 0.;

    size_t O_index = T.shape()[0]-N_GC-1; //radial index of exterior cell (whose value we are setting)
    double DeltarOuter = 0.5*( dm.Deltar[O_index-1] + dm.Deltar[O_index] );
    double B12 = tparams.B_pol*B_0/1e12; //initial magnetic field at poles in 10^12 G

    static double g_f = dm.n_t[O_index-1]*dm.n_t[O_index-1]*pow(r_s/r_b,2.); //factor of exp(nu)*(R_out/R_b)^2
    static double x_g = (1.-1./(dm.n_r[O_index-1]*dm.n_r[O_index-1]))*dm.r[O_index-1]/r_s;
    static double a_g; //general relativistic correction parameter for magnetic field inclination angle. Equals 1/2 in nonrelativistic limit x_g = 0. Only used for T_upper = "magnetic_P01"
    // if( dm.n_r[O_index-1] <= 1+1e-10 ) a_g = 0.5;
    // else a_g = -( (1.-x_g)*log(1.-x_g)+x_g-0.5*x_g*x_g )/( (log(1.-x_g)+x_g+0.5*x_g*x_g)*sqrt(1.-x_g) );

    if(tparams.T_upper == "non-magnetic_G83"){
        Tb8 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e8; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
        T_sGudmundsson1983(Tb8, g_s14, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Gudmundsson, Pethick and Epstein, ApJ 272, 286 (1983)
    }
    else if(tparams.T_upper == "non-magnetic_P01"){
        Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
        T_sPotekhin2001(Tb9, g_s14, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
    }
    else if(tparams.T_upper == "magnetic_P01"){
        Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
        // theta_B = atan(a_g*tan(dm.theta[j])); //magnetic field inclination angle in radians. Equals the local angle between the magnetic field and the radial direction. Includes general relativistic correction
        theta_B = acos(mC.phBhat[0][O_index][j]); //magnetic field inclination angle in radians
        T_sPotekhinMag2001(Tb9, B12, g_s14, theta_B, Ts, Ts_star, dTsdTb); // T_s - T_b relation for non-magnetic heat blanketing envelope from Potekhin & Yakovlev, Astron. Astrophys. 374, 213 (2001)
    }
    else if(tparams.T_upper == "magnetic_P15"){
        Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
        Tb9 = T[O_index-1][j]/dm.n_t[O_index-1]*T_0/1.e9; //local temperature at inner edge of heat-blanketing envelope (outer edge of simulation domain) in units of 10^9 K
        theta_mc = dm.theta[j];  //magnetic colatitude in radians
        T_sPotekhin2015(Tb9, B12, g_s14, theta_mc, Ts, Ts_star, dTsdTb);// T_s - T_b relation for heat blanketing envelope is from Potekhin, Pons & Page, Space Sci. Rev. 191, 239 (2015), Appendix B.
    }
    else {
        std::cout << "Outer thermal boundary condition not set" << std::endl;
    }

    a = g_f*sigma_SBbar*Ts_star*Ts_star*Ts_star*( Ts_star - 2.*dTsdTb/dm.n_t[O_index-1]*(T[O_index][j] + T[O_index-1][j]) );
    b = 2.*g_f*sigma_SBbar*Ts_star*Ts_star*Ts_star*dTsdTb/dm.n_t[O_index-1];

    if(tparams.conductivity_anisotropy == true){
        kappa_r = 0.5*( (tC.kappa_delta[O_index-1][j]+tC.kappa_delta[O_index][j])*mC.thBhat[0][O_index][j]*mC.thBhat[0][O_index][j] + (tC.kappa_perp[O_index-1][j]+tC.kappa_perp[O_index][j]) );
        kappa_th = 0.5*( (tC.kappa_delta[O_index-1][j]+tC.kappa_delta[O_index][j])*mC.thBhat[0][O_index][j]*mC.thBhat[1][O_index][j] - (tC.kappa_H[O_index-1][j]+tC.kappa_H[O_index][j])*mC.thBhat[2][O_index][j] );
        RadialTerm = dm.n_t[O_index-1]*kappa_r/(DeltarOuter*dm.n_r[O_index-1]); //K_r
        TangentialTerm = 0*dm.n_t[O_index-1]*kappa_th/r_b*0.5/tan(dm.Deltatheta/2.)*( T[O_index-1][j+1] - T[O_index-1][j-1] )/2.; //K_theta
    }
    else{
        kappa_r = 0.5*(tC.kappa[O_index-1][j]+tC.kappa[O_index][j]); //if only one isotropic thermal conductivity, do not need to take an average
        RadialTerm = dm.n_t[O_index-1]*kappa_r/(DeltarOuter*dm.n_r[O_index-1]); //K_r
        TangentialTerm = 0.;
    }

    //T_{i=O_index} = b_n_term*T_{i=O_index-1} + d_n_term
    b_n_term = ( RadialTerm - b )/( RadialTerm + b );
    d_n_term = -( TangentialTerm + a )/( RadialTerm + b );

    return;
}

/*
        Computes predicted temperature in first ghost cell at boundary between decomposed domains.
        Evolves T in theta direction using explicit method. Only includes unmixed d^2/dtheta^2 and spatial curvature term, as does
        the function "Compute_T_IMEX" which uses "T_InteriorBC_theta".
        Inputs: j_in: Indices in theta direction of set of ghost cells whose temperature to approximate
                T: redshifted temperature in reduced units. This is intermediate temperature after ADI evaluation in r-direction.
                thC: ThermCoeffs object containing specific heat capacity c_v and neutrino emissivity q_nu
                trC: TransCoeffs object containing thermal conductivity kappa
                mC: MagCoeffs object containing properties of magnetic field and magnetization
                dm, process: Domain and Process objects containing information about the simulation domain and the current process
                N_GC: number of ghost cells
                a_ii: diagonal entry in Butcher tableau for this implicit RK step
        Output: T_f: updated redshifted temperature in reduced units, with first ghost cell values updated explicitly (forward Euler) as a first prediction
*/
void T_InteriorBC_theta(size_t j_in, ScalarField & T, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const Domain & dm, size_t N_GC, double a_ii, ScalarField & T_f) {

    double Ftheta, QT; //flux difference function in r and theta directions and cross term
    size_t Nr = T.shape()[0]; //radial extent of T_f, including ghost cells
    ScalarField q_th(boost::extents[Nr][2]); //temporary VectorField object to hold necessary heat fluxes into/out of ghost cell under consideration. Only the part of q_th containing dT/dtheta

    if(trC.conductivity_anisotropy == true){
        double kappa_delta_jmh, kappa_perp_jmh;
        for(size_t i=N_GC; i<Nr-N_GC-1; i++){
            for(size_t j=j_in; j<j_in+2; j++){
                kappa_delta_jmh = 2./( 1./trC.kappa_delta[i][j-1] + 1./trC.kappa_delta[i][j] );
                kappa_perp_jmh = 2./( 1./trC.kappa_perp[i][j-1] + 1./trC.kappa_perp[i][j] );

                q_th[i][j-j_in] = -dm.n_t[i]*( kappa_delta_jmh*mC.rBhat[1][i][j]*mC.rBhat[1][i][j] + kappa_perp_jmh )*gradth(T,dm.Deltatheta*dm.r[i],i,j);
            }
        }
    }
    else{
        double kappa_jmh;
        for(size_t i=N_GC; i<Nr-N_GC-1; i++){
            for(size_t j=j_in; j<j_in+2; j++){
                kappa_jmh = 2./( 1./trC.kappa[i][j-1] + 1./trC.kappa[i][j] );
                q_th[i][j-j_in] = -dm.n_t[i]*kappa_jmh*gradth(T,dm.Deltatheta*dm.r[i],i,j);
            }
        }
    }

    //Divide by overall factor of c_v[i][j]
    for(size_t i=N_GC; i<Nr-N_GC-1; i++){
        Ftheta = sin(dm.theta[j_in]+0.5*dm.Deltatheta)*q_th[i][1] - sin(dm.theta[j_in]-0.5*dm.Deltatheta)*q_th[i][0];
        QT = ( -Ftheta/(dm.r[i]*2.*sin(dm.theta[j_in])*sin(dm.Deltatheta/2.)) )/thC.c_v[i][j_in];

        T_f[i][j_in] = T[i][j_in] + dm.Deltat*a_ii*QT;
    }

}
