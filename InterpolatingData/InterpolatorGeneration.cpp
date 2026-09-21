#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_erf.h>
#include <typeinfo>

/*
        Generates interpolating functions required for finite-temperature formulae for Landau quantization-induced
        magnetization and its derivatives

        Uses GSL spline library to generate the interpolating functions, then serializes them so they can be loaded
        and stored for repeated use within the ElectronMHDAxisymmetricSpherical

*/

const double pi = 3.141592653589793238463;
const double c = 29979245800; //speed of light in cm/s
const double hbarc = 197.3269804; //hbar times c in MeV*fm
const double unit_e = 4.80320425e-10; //elementary charge in units of statcoulomb
const double k_B = 8.61733034e-11; //Boltzmann constant in MeV/K
const double k_Bcgs = 1.380649e-16; //Boltzmann constant in cgs units (erg/K)
const double M_e = 0.51099895000; //electron mass in MeV
const double M_N = 938.92; //mean nucleon mass in MeV (938.91875434 MeV)
const double M_n = 939.56542052; //(bare) neutron mass in MeV
const double M_p = 938.27208943; //(bare) proton mass in MeV
const double MeVtoErg = 1/6.2415e5; //conversion factor from MeV to erg
const double alpha_e = 0.0072973525693; //electromagnetic fine structure constant (dimensionless)
const double B_crit = 4.41400564e13; //quantum critical magnetic field (Schwinger field) in G
const double eB_crit = M_e*M_e; //critical magnetic field times elementary charge in MeV^2

const double GaussConverter = 6.241509074e-8/hbarc; //conversion factor between 1 statCoulomb*Gauss to __ MeV/fm/hbarc = ___ fm^-2
const double GaussConverter2 = 6.241509074e-8*hbarc; //conversion factor between 1 statCoulomb*Gauss to __ MeV/fm*hbarc = ___ MeV^2
const double GaussConverter3 = 4.002719868e16; //conversion factor: 1 sqrt(MeV/fm^3) = 4.002719868e16 G
const double GaussConverter4 = GaussConverter3/pow(hbarc,1.5); //conversion factor: 1 MeV^2 = 4.002719868e16/(197.3269804)^(3/2) G = 1.444027592e13 G (hbarc=c=1 all energy units)

struct LowTInterp_params{
    double a, b; // a = sqrt(M_e^2+2*eB*n)/T = M_n/T, b = mu_e/T
};

struct HighTInterp_params{
    double a, b, A; // a = mu_e/M_e, b = eB/M_e^2, A = (mu_e^2-M_e^2)/eB
};

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

double G_1Integrand(double x, void * params);
double G_2Integrand(double x, void * params);
double NH_2Integrand(double x, void * params);
double I_1Integrand(double x, void * params);
double NI_2Integrand(double x, void * params);
double h_1Integrand(double y, void * params);
double h_2Integrand(double y, void * params);

int main(){

    /*
        Part 1: Interpolating functions for polylogarithms of half-integer order.

        Li_{-1/2}(-exp(z)), Li_{1/2}(-exp(z)), Li_{3/2}(-exp(z)), Li_{5/2}(-exp(z)), Li_{7/2}(-exp(z))

        Because of a lack of C++ libraries for computing these, we perform the initial computation of these functions for
        -50 < z < 50 using Python's mpmath library, generating a data table Polylogarithms.dat, and then interpolate using GSL based on this
        data table

        Now included within "Compute_OmegaB_Interpolators" function defined in microphysics_Bdep.cpp

        Used in low-temperature approximations to magnetization-related thermodynamic functions

    */

//    std::ostringstream filename;
//	filename << "Polylogarithms.dat";
//	std::fstream Datafile;
//    Datafile.open(filename.str().c_str(),std::fstream::in);
//	std::string line, dummyLine;
//    size_t N_dat = 0; //number of rows in data table. Initialize to zero
//	getline(Datafile, dummyLine); //skips header (1 line)
//	while (std::getline(Datafile, line))
//		++N_dat;
//
//	Datafile.close();
//
//    //vectors to store values from data table
//	std::vector<double> z_InterpVec; //z, related to arguments of polylogarithms x by z = ln(-x)
//	std::vector<double> PolyLogN1_2_InterpVec; //Li_{-1/2}(-exp(z))
//	std::vector<double> PolyLog1_2_InterpVec; //Li_{-1/2}(-exp(z))
//	std::vector<double> PolyLog3_2_InterpVec; //Li_{-1/2}(-exp(z))
//    std::vector<double> PolyLog5_2_InterpVec; //Li_{-1/2}(-exp(z))
//	std::vector<double> PolyLog7_2_InterpVec; //Li_{-1/2}(-exp(z))
//
//    //temporary variables to hold data table values
//    double z_dt, PL1_dt, PL2_dt, PL3_dt, PL4_dt, PL5_dt;
//    double dummy; //dummy variable to avoid saving all quantities from loaded data table
//
//	Datafile.open(filename.str().c_str(),std::fstream::in);
//	std::getline(Datafile, dummyLine); //skips header (1 line)
//	//Gets background TOV solution values from data table
//	for(size_t i=0; i<N_dat; i++)
//	{
//		Datafile >> z_dt >> PL1_dt >> PL2_dt >> PL3_dt >> PL4_dt >> PL5_dt;
//
//		z_InterpVec.push_back(z_dt);
//		PolyLogN1_2_InterpVec.push_back(PL1_dt);
//		PolyLog1_2_InterpVec.push_back(PL2_dt);
//		PolyLog3_2_InterpVec.push_back(PL3_dt);
//		PolyLog5_2_InterpVec.push_back(PL4_dt);
//		PolyLog7_2_InterpVec.push_back(PL5_dt);
//
//	}
//	Datafile.close();
//
//	//converts vectors of data to arrays to use in interpolation
//	double* z_Interp = &z_InterpVec[0];
//	double* PolyLogN1_2_Interp = &PolyLogN1_2_InterpVec[0];
//	double* PolyLog1_2_Interp = &PolyLog1_2_InterpVec[0];
//	double* PolyLog3_2_Interp = &PolyLog3_2_InterpVec[0];
//    double* PolyLog5_2_Interp = &PolyLog5_2_InterpVec[0];
//    double* PolyLog7_2_Interp = &PolyLog7_2_InterpVec[0];
//
//	//Generates spline and acc objects for each quantity to interpolate them and assigns them to DataInterpolation object
//	gsl_interp_accel *PolyLogN1_2_acc = gsl_interp_accel_alloc ();
//	gsl_spline *PolyLogN1_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
//	gsl_spline_init(PolyLogN1_2, z_Interp, PolyLogN1_2_Interp, N_dat);
//
//	gsl_interp_accel *acc_PolyLog1_2 = gsl_interp_accel_alloc ();
//	gsl_spline *PolyLog1_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
//	gsl_spline_init(PolyLog1_2, z_Interp, PolyLog1_2_Interp, N_dat);
//
//	gsl_interp_accel *acc_PolyLog3_2 = gsl_interp_accel_alloc ();
//	gsl_spline *PolyLog3_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
//	gsl_spline_init(PolyLog3_2, z_Interp, PolyLog3_2_Interp, N_dat);
//
//    gsl_interp_accel *acc_PolyLog5_2 = gsl_interp_accel_alloc ();
//	gsl_spline *PolyLog5_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
//	gsl_spline_init(PolyLog5_2, z_Interp, PolyLog5_2_Interp, N_dat);
//
//    gsl_interp_accel *acc_PolyLog7_2 = gsl_interp_accel_alloc ();
//	gsl_spline *PolyLog7_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
//	gsl_spline_init(PolyLog7_2, z_Interp, PolyLog7_2_Interp, N_dat);

	/*
        Part 2: 2D interpolating function data for the two-argument functions

        G_1(a,b), G_2(a,b), -H_2(a,b), I_1(a,b), -I_2(a,b)

        where a = sqrt(M_e^2+2*eB*n)/T = M_n/T, b = mu_e/T

        Use approximate forms for a > 50, b > 50, so this is upper limit for each.
        At lower limit, lowest values of mu_e and sqrt(M_e^2+2*eB*n) is M_e (in strongly quantizing limit n=0 for latter)
        So set lower limit of range to be m_e/(2 MeV) ~ 0.2 where 2 MeV ~ 2.33e10 K (beyond highest temperature we consider)

        Used in low-temperature approximations to magnetization-related thermodynamic functions

	*/

    std::vector<double> a;
    std::vector<double> b;

    size_t num_interp = 200; //number of values of a and b to use to generate interpolation functions

    double abmin = 0.5;//0.2;
    double abmax = 51; //include a bit above maximum value that will be consider to avoid possible interpolation errors.

    a = linspace(abmin,abmax,num_interp);
    b = linspace(abmin,abmax,num_interp);

    //Functions to integrate to determine data to be interpolated
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(3000);
    double result, error;

    struct LowTInterp_params paramsLT;

    gsl_function G_1f;
    G_1f.function = &G_1Integrand;
    G_1f.params = &paramsLT;

    gsl_function G_2f;
    G_2f.function = &G_2Integrand;
    G_2f.params = &paramsLT;

    gsl_function NH_2f;
    NH_2f.function = &NH_2Integrand;
    NH_2f.params = &paramsLT;

    gsl_function I_1f;
    I_1f.function = &I_1Integrand;
    I_1f.params = &paramsLT;

    gsl_function NI_2f;
    NI_2f.function = &NI_2Integrand;
    NI_2f.params = &paramsLT;

    //Create data files for integrated values
    std::fstream DataG_1;
    DataG_1.open("G_1.dat",std::fstream::out);
    std::fstream DataG_2;
    DataG_2.open("G_2.dat",std::fstream::out);
    std::fstream DataNH_2;
    DataNH_2.open("NH_2.dat",std::fstream::out);
    std::fstream DataI_1;
    DataI_1.open("I_1.dat",std::fstream::out);
    std::fstream DataNI_2;
    DataNI_2.open("NI_2.dat",std::fstream::out);

    //Extra space at top-left corner of data files
    DataG_1 << "             ";
    DataG_2 << "             ";
    DataNH_2 << "             ";
    DataI_1 << "             ";
    DataNI_2 << "             ";

    for(size_t j=0; j<b.size(); j++){
        DataG_1 << std::setprecision(10) << b[j] << "   ";
        DataG_2 << std::setprecision(10) << b[j] << "   ";
        DataNH_2 << std::setprecision(10) << b[j] << "   ";
        DataI_1 << std::setprecision(10) << b[j] << "   ";
        DataNI_2 << std::setprecision(10) << b[j] << "   ";
    }
    DataG_1 << std::endl;
    DataG_2 << std::endl;
    DataNH_2 << std::endl;
    DataI_1 << std::endl;
    DataNI_2 << std::endl;

    //Integrate functions. Interpolation performed on integrals
    std::cout << "Performing integrals for G_1 through -I_2" << std::endl;
    for(size_t i=0; i<a.size(); i++){

        DataG_1 << a[i] << "   ";
        DataG_2 << a[i] << "   ";
        DataNH_2 << a[i] << "   ";
        DataI_1 << a[i] << "   ";
        DataNI_2 << a[i] << "   ";

        for(size_t j=0; j<b.size(); j++){

            std::cout << "a = " << a[i] << ", b = " << b[j] << std::endl;

            paramsLT = {a[i],b[j]};

            gsl_integration_qagiu(&G_1f, 1., 0, 1e-7, 1000, w, &result, &error);
            DataG_1 << result << "   ";

            gsl_integration_qagiu(&G_2f, 1., 0, 1e-7, 1000, w, &result, &error);
            DataG_2 << result << "   ";

            gsl_integration_qagiu(&NH_2f, 1., 0, 1e-7, 1000, w, &result, &error);
            DataNH_2 << result << "   ";

            gsl_integration_qagiu(&I_1f, 1., 0, 1e-7, 1000, w, &result, &error);
            DataI_1 << result << "   ";

            gsl_integration_qagiu(&NI_2f, 1., 0, 1e-7, 1000, w, &result, &error);
            DataNI_2 << result << "   ";

        }

        DataG_1 << std::endl;
        DataG_2 << std::endl;
        DataNH_2 << std::endl;
        DataI_1 << std::endl;
        DataNI_2 << std::endl;

    }
    DataG_1.close();
    DataG_2.close();
    DataNH_2.close();
    DataI_1.close();
    DataNI_2.close();

    std::cout << "Integrals for G_1 through -I_2 complete" << std::endl;

    //Delete contents of vectors a and b, since will reuse them with different values
    a.clear();
    b.clear();

    /*
        Part 3: 2D interpolating functions for the two-argument functions

        h_1(a,b), h_2(a,b)

        where a = mu_e/M_e, b = eB/M_e^2

        Used in high-temperature approximations to magnetization-related thermodynamic functions

    */

    num_interp = 60;

    double a_min = 10., a_max = 170.;
    double b_min = 1., b_max = 5000.;

    //logarithmically-spaced values of a and b
    for(size_t i = 0; i<num_interp; i++){
        a.push_back( pow( 10., log10(a_min)+(log10(a_max)-log10(a_min))/(double(num_interp-1))*double(i) ) );
        b.push_back( pow( 10., log10(b_min)+(log10(b_max)-log10(b_min))/(double(num_interp-1))*double(i) )  );
    }

    const double zeta3_2 = 2.61237534869; //zeta(3/2)

    //Functions to integrate to determine data to be interpolated
    struct HighTInterp_params paramsHT;

    gsl_function h_1f;
    h_1f.function = &h_1Integrand;
    h_1f.params = &paramsHT;

    gsl_function h_2f;
    h_2f.function = &h_2Integrand;
    h_2f.params = &paramsHT;

    //Create data files for integrated values
    std::fstream Datah_1;
    Datah_1.open("h_1.dat",std::fstream::out);
    std::fstream Datah_2;
    Datah_2.open("h_2.dat",std::fstream::out);

    //Extra space at top-left corner of data files
    Datah_1 << "             ";
    Datah_2 << "             ";

    //Adds b = 0 entry, where the integrals are both zero
    Datah_1 << std::setprecision(10) << 0. << "   ";
    Datah_2 << 0. << "   ";

    for(size_t j=0; j<b.size(); j++){
        Datah_1 << b[j] << "   ";
        Datah_2 << b[j] << "   ";
    }
    Datah_1 << std::endl;
    Datah_2 << std::endl;
    double result2, error2;

    //Integrate functions. Interpolation performed on integrals.
    //Split the integral into two parts due to singularity at y=0.

    double pts[3] = {0.,0.,1.};
    size_t numpts = 3;
    std::cout << "Performing integrals for h_1 and h_2" << std::endl;
    //Adds row of a=1 entries (both integrals evaluate to zero for a=1)
    Datah_1 << 1. << "   ";
    Datah_2 << 1. << "   ";
    for(size_t i=0; i<b.size()+1; i++){
        Datah_1 << 0. << "   ";
        Datah_2 << 0. << "   ";
    }
    Datah_1 << std::endl;
    Datah_2 << std::endl;

    for(size_t i=0; i<a.size(); i++){

        Datah_1 << a[i] << "   " << 0. << "   ";
        Datah_2 << a[i] << "   " << 0. << "   ";

        for(size_t j=0; j<b.size(); j++){

            std::cout << "a = " << a[i] << ", b = " << b[j] << std::endl;

            paramsHT = {a[i],b[j]};

            gsl_integration_qagp(&h_1f, pts, numpts, 0, 1e-7, 3000, w, &result, &error);
            gsl_integration_qagiu(&h_1f, 1., 0, 1e-7, 3000, w, &result2, &error2);
            Datah_1 << M_e*( result + result2 + (a[i] - 1.)*sqrt(2./pi)*zeta3_2 ) << "   ";

            gsl_integration_qagp(&h_2f, pts, numpts, 0, 1e-7, 3000, w, &result, &error);
            gsl_integration_qagiu(&h_2f, 1., 0, 1e-7, 3000, w, &result2, &error2);
            Datah_2 << M_e*( result + result2 ) << "   ";

        }

        Datah_1 << std::endl;
        Datah_2 << std::endl;
    }
    Datah_1.close();
    Datah_2.close();

    std::cout << "Integrals for h_1 and h_2 complete" << std::endl;

    return 0;
}

double G_1Integrand(double x, void * params){

    struct LowTInterp_params *p = (struct LowTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( a*x-b > 300. ) result = x*x/sqrt(x*x-1.)*exp(-(a*x-b));
    else result = x*x/sqrt(x*x-1.)/( exp(a*x-b) + 1. );

    return result;
}

double G_2Integrand(double x, void * params){

    struct LowTInterp_params *p = (struct LowTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( a*x-b > 300. ) result = 1./sqrt(x*x-1.)*exp(-(a*x-b));
    else result = 1./sqrt(x*x-1.)/( exp(a*x-b) + 1. );

    return result;
}


double NH_2Integrand(double x, void * params){

    struct LowTInterp_params *p = (struct LowTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( a*x-b > 300. ) result = a*x/sqrt(x*x-1.)*exp(-(a*x-b));
    else result = a*x/sqrt(x*x-1.)*exp(a*x-b)/pow( exp(a*x-b) + 1.,2. );

    return result;
}


double I_1Integrand(double x, void * params){

    struct LowTInterp_params *p = (struct LowTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( a*x-b > 300. ) result = x/sqrt(x*x-1.)*exp(-(a*x-b));
    else result = x/sqrt(x*x-1.)/( exp(a*x-b) + 1. );

    return result;
}


double NI_2Integrand(double x, void * params){

    struct LowTInterp_params *p = (struct LowTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( a*x-b > 300. ) result = a/sqrt(x*x-1.)*exp(-(a*x-b));
    else result = a/sqrt(x*x-1.)*exp(a*x-b)/pow( exp(a*x-b) + 1.,2. );

    return result;
}

double h_1Integrand(double y, void * params){

    struct HighTInterp_params *p = (struct HighTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( y < 1e-4 ){
        result = -( pow(a,3.) - 3.*a + 2. )/(9.*b)*sqrt(y) + ( 3.*pow(a,5.) - 10.*pow(a,3.) + 15.*a - 8. )/(90.*b*b)*pow(y,1.5);
    }
    else if( y/b > 2000. ){
        result = (y-1)/pow(y,2.5)*( 1. - a + b/(2.*y)*( 1. - exp( (1.-a*a)*y/b)/a ) );
    }
    else{
        result = (y/tanh(y) - 1.)/pow(y,2.5)*( 1. - a + sqrt(pi*b/y)/2.*( exp( y/b + gsl_sf_log_erfc(sqrt(y/b)) ) - exp( y/b + gsl_sf_log_erfc(a*sqrt(y/b)) )  ) );
    }

    return result;
}

double h_2Integrand(double y, void * params){

    struct HighTInterp_params *p = (struct HighTInterp_params *) params;
    double a = p -> a;
    double b = p -> b;

    double result;

    if( y < 1e-4 ){
        result = ( a - 1. )/3.*sqrt(y) - 1./(9.*b)*( pow(a,3.) - 3.*a + 2. )*pow(y,1.5);
    }
    else if( y/b > 20000. ){
        result = b*(y - 1.)/(2.*pow(y,2.5))*( 1. - exp( (1.-a*a)*y/b )/a );
    }
    else{
        result = (y/tanh(y) - 1.)*sqrt(pi*b)/(2.*y*y)*( exp( y/b + gsl_sf_log_erfc(sqrt(y/b)) ) - exp( y/b + gsl_sf_log_erfc(a*sqrt(y/b)) )  );
    }

    return result;
}
