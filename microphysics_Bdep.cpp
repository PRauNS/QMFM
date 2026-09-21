#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "common.h"
#include "microphysics.h"
#include "microphysics_Bdep.h"

/*
        Generate interpolating functions using data from "Polylogarithm_DataTable.py" and "InterpolatorGeneration.cpp"
        These are used in computing the finite temperature magnetization and its derivatives
        Input: OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
*/
void Compute_OmegaB_Interpolators(OmegaB_deriv_Interpolation & OmB_Interpolators)
{
    /*
        Part 1: Interpolating functions for polylogarithms of half-integer order.

        Li_{-1/2}(-exp(z)), Li_{1/2}(-exp(z)), Li_{3/2}(-exp(z)), Li_{5/2}(-exp(z)), Li_{7/2}(-exp(z))

        Because of a lack of C++ libraries for computing these, we perform the initial computation of these functions for
        -50 < z < 50 using Python's mpmath library, generating a data table Polylogarithms.dat, and then interpolate
        using GSL based on this data table

        Used in low-temperature approximations to magnetization-related thermodynamic functions
    */

    std::ostringstream filename;
	filename << "InterpolatingData/Polylogarithms.dat";
	std::fstream Datafile;
    Datafile.open(filename.str().c_str(),std::fstream::in);
	std::string line, dummyLine;
    size_t N_dat = 0; //number of rows in data table. Initialize to zero
	getline(Datafile, dummyLine); //skips header (1 line)
	while (std::getline(Datafile, line))
		++N_dat;

	Datafile.close();

    //vectors to store values from data table
	std::vector<double> z_InterpVec; //z, related to arguments of polylogarithms x by z = ln(-x)
	std::vector<double> PolyLogN1_2_InterpVec; //Li_{-1/2}(-exp(z))
	std::vector<double> PolyLog1_2_InterpVec; //Li_{-1/2}(-exp(z))
	std::vector<double> PolyLog3_2_InterpVec; //Li_{-1/2}(-exp(z))
    std::vector<double> PolyLog5_2_InterpVec; //Li_{-1/2}(-exp(z))
	std::vector<double> PolyLog7_2_InterpVec; //Li_{-1/2}(-exp(z))

    //temporary variables to hold data table values
    double z_dt, PL1_dt, PL2_dt, PL3_dt, PL4_dt, PL5_dt;

	Datafile.open(filename.str().c_str(),std::fstream::in);
	std::getline(Datafile, dummyLine); //skips header (1 line)
	for(size_t i=0; i<N_dat; i++){
		Datafile >> z_dt >> PL1_dt >> PL2_dt >> PL3_dt >> PL4_dt >> PL5_dt;

		z_InterpVec.push_back(z_dt);
		PolyLogN1_2_InterpVec.push_back(PL1_dt);
		PolyLog1_2_InterpVec.push_back(PL2_dt);
		PolyLog3_2_InterpVec.push_back(PL3_dt);
		PolyLog5_2_InterpVec.push_back(PL4_dt);
		PolyLog7_2_InterpVec.push_back(PL5_dt);
	}
	Datafile.close();

	//converts vectors of data to arrays to use in interpolation
	double* z_Interp = &z_InterpVec[0];
	double* PolyLogN1_2_Interp = &PolyLogN1_2_InterpVec[0];
	double* PolyLog1_2_Interp = &PolyLog1_2_InterpVec[0];
	double* PolyLog3_2_Interp = &PolyLog3_2_InterpVec[0];
    double* PolyLog5_2_Interp = &PolyLog5_2_InterpVec[0];
    double* PolyLog7_2_Interp = &PolyLog7_2_InterpVec[0];

	//Generates spline and acc objects for each quantity to interpolate them and assigns them to DataInterpolation object
	gsl_interp_accel *PolyLogN1_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *PolyLogN1_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(PolyLogN1_2, z_Interp, PolyLogN1_2_Interp, N_dat);

	gsl_interp_accel *PolyLog1_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *PolyLog1_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(PolyLog1_2, z_Interp, PolyLog1_2_Interp, N_dat);

	gsl_interp_accel *PolyLog3_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *PolyLog3_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(PolyLog3_2, z_Interp, PolyLog3_2_Interp, N_dat);

    gsl_interp_accel *PolyLog5_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *PolyLog5_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(PolyLog5_2, z_Interp, PolyLog5_2_Interp, N_dat);

    gsl_interp_accel *PolyLog7_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *PolyLog7_2 = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(PolyLog7_2, z_Interp, PolyLog7_2_Interp, N_dat);

    //Stores interpolation functions in instance of class "OmegaB_deriv_interpolators"
    //Also creates accelerator objects for each interpolator and stores them in this class instance.

    OmB_Interpolators.PolyLogN1_2 = PolyLogN1_2;
    OmB_Interpolators.PolyLogN1_2_acc = PolyLogN1_2_acc;
    OmB_Interpolators.PolyLog1_2 = PolyLog1_2;
    OmB_Interpolators.PolyLog1_2_acc = PolyLog1_2_acc;
    OmB_Interpolators.PolyLog3_2 = PolyLog3_2;
    OmB_Interpolators.PolyLog3_2_acc = PolyLog3_2_acc;
    OmB_Interpolators.PolyLog5_2 = PolyLog5_2;
    OmB_Interpolators.PolyLog5_2_acc = PolyLog5_2_acc;
    OmB_Interpolators.PolyLog7_2 = PolyLog7_2;
    OmB_Interpolators.PolyLog7_2_acc = PolyLog7_2_acc;

	/*
        Part 2: 2D interpolating functions for the two-argument functions

        G_1(a,b), G_2(a,b), -H_2(a,b), I_1(a,b), -I_2(a,b)

        where a = sqrt(M_e^2+2*eB*n)/T = M_n/T, b = mu_e/T

        G_1.dat through NI_2.dat are organized by row values 0.1 < a < 50 and column values 0.1 < b < 50

        Used in low-temperature approximations to magnetization-related thermodynamic functions.

        Load in 2D data tables G_1.dat, G_2.dat, NH_2.dat, I_1.dat and NI_2.dat and then interpolate
        using GSL based on these data tables. For numerical accuracy, better to interpolate the logarithm
        of these functions, and then compute them using G_1 = exp(LogG_1Interp), H_2 = -exp(LogNH_2Interp), etc.
	*/

    std::ostringstream filenameG_1;
	filenameG_1 << "InterpolatingData/G_1.dat";
	std::fstream DatafileG_1;
    DatafileG_1.open(filenameG_1.str().c_str(),std::fstream::in);
    N_dat = 0; //number of rows in data table. Initialize to zero
	while (std::getline(DatafileG_1, line))
		++N_dat;
	DatafileG_1.close();

    //vectors to store values from data table
	std::vector<double> a; //a values. Same for each data set in Part 2
	std::vector<double> b; //b values. Same for each data set in Part 2
    double logG_1[N_dat-1][N_dat-1];
    double logG_2[N_dat-1][N_dat-1];
    double logNH_2[N_dat-1][N_dat-1];
    double logI_1[N_dat-1][N_dat-1];
    double logNI_2[N_dat-1][N_dat-1];

    //temporary variables to hold data table values
    double a_dt, b_dt, G_1_dt, G_2_dt, NH_2_dt, I_1_dt, NI_2_dt;

	DatafileG_1.open(filenameG_1.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                DatafileG_1 >> b_dt;
                b.push_back(b_dt);
            }
        }
        else{
            DatafileG_1 >> a_dt;
            a.push_back(a_dt);
            for(size_t j=0; j<N_dat-1; j++){
                DatafileG_1 >> G_1_dt;
                logG_1[i-1][j] = log(G_1_dt);
            }
        }
	}
	DatafileG_1.close();

    std::ostringstream filenameG_2;
	filenameG_2 << "InterpolatingData/G_2.dat";
	std::fstream DatafileG_2;
    DatafileG_2.open(filenameG_2.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                DatafileG_2 >> b_dt;
            }
        }
        else{
            DatafileG_2 >> a_dt;
            for(size_t j=0; j<N_dat-1; j++){
                DatafileG_2 >> G_2_dt;
                logG_2[i-1][j] = log(G_2_dt);
            }
        }
	}
	DatafileG_2.close();

    std::ostringstream filenameNH_2;
	filenameNH_2 << "InterpolatingData/NH_2.dat";
	std::fstream DatafileNH_2;
    DatafileNH_2.open(filenameNH_2.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                DatafileNH_2 >> b_dt;
            }
        }
        else{
            DatafileNH_2 >> a_dt;
            for(size_t j=0; j<N_dat-1; j++){
                DatafileNH_2 >> NH_2_dt;
                logNH_2[i-1][j] = log(NH_2_dt);
            }
        }
	}
	DatafileNH_2.close();

    std::ostringstream filenameI_1;
	filenameI_1 << "InterpolatingData/I_1.dat";
	std::fstream DatafileI_1;
    DatafileI_1.open(filenameI_1.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                DatafileI_1 >> b_dt;
            }
        }
        else{
            DatafileI_1 >> a_dt;
            for(size_t j=0; j<N_dat-1; j++){
                DatafileI_1 >> I_1_dt;
                logI_1[i-1][j] = log(I_1_dt);
            }
        }
	}
	DatafileI_1.close();

    std::ostringstream filenameNI_2;
	filenameNI_2 << "InterpolatingData/NI_2.dat";
	std::fstream DatafileNI_2;
    DatafileNI_2.open(filenameNI_2.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                DatafileNI_2 >> b_dt;
            }
        }
        else{
            DatafileNI_2 >> a_dt;
            for(size_t j=0; j<N_dat-1; j++){
                DatafileNI_2 >> NI_2_dt;
                logNI_2[i-1][j] = log(NI_2_dt);
            }
        }
	}
	DatafileNI_2.close();

    //Declare interpolating functions
    const gsl_interp2d_type *T_G_1 = gsl_interp2d_bicubic;
    const gsl_interp2d_type *T_G_2 = gsl_interp2d_bicubic;
    const gsl_interp2d_type *T_NH_2 = gsl_interp2d_bicubic;
    const gsl_interp2d_type *T_I_1 = gsl_interp2d_bicubic;
    const gsl_interp2d_type *T_NI_2 = gsl_interp2d_bicubic;

    double *z_G_1 = (double *) malloc(a.size() * b.size() * sizeof(double));
    double *z_G_2 = (double *) malloc(a.size() * b.size() * sizeof(double));
    double *z_NH_2 = (double *) malloc(a.size() * b.size() * sizeof(double));
    double *z_I_1 = (double *) malloc(a.size() * b.size() * sizeof(double));
    double *z_NI_2 = (double *) malloc(a.size() * b.size() * sizeof(double));

    gsl_interp_accel *LogG_1_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *LogG_1_acc_b = gsl_interp_accel_alloc();
    gsl_interp_accel *LogG_2_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *LogG_2_acc_b = gsl_interp_accel_alloc();
    gsl_interp_accel *LogNH_2_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *LogNH_2_acc_b = gsl_interp_accel_alloc();
    gsl_interp_accel *LogI_1_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *LogI_1_acc_b = gsl_interp_accel_alloc();
    gsl_interp_accel *LogNI_2_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *LogNI_2_acc_b = gsl_interp_accel_alloc();

    gsl_spline2d *spline_LogG_1 = gsl_spline2d_alloc(T_G_1, a.size(), b.size());
    gsl_spline2d *spline_LogG_2 = gsl_spline2d_alloc(T_G_2, a.size(), b.size());
    gsl_spline2d *spline_LogNH_2 = gsl_spline2d_alloc(T_NH_2, a.size(), b.size());
    gsl_spline2d *spline_LogI_1 = gsl_spline2d_alloc(T_I_1, a.size(), b.size());
    gsl_spline2d *spline_LogNI_2 = gsl_spline2d_alloc(T_NI_2, a.size(), b.size());

    //Integrate functions. Interpolation performed on integrals
    for(size_t i=0; i<a.size(); i++){
        for(size_t j=0; j<b.size(); j++){
            gsl_spline2d_set(spline_LogG_1, z_G_1, i, j, logG_1[i][j]);
            gsl_spline2d_set(spline_LogG_2, z_G_2, i, j, logG_2[i][j]);
            gsl_spline2d_set(spline_LogNH_2, z_NH_2, i, j, logNH_2[i][j]);
            gsl_spline2d_set(spline_LogI_1, z_I_1, i, j, logI_1[i][j]);
            gsl_spline2d_set(spline_LogNI_2, z_NI_2, i, j, logNI_2[i][j]);
        }
    }

    //Interpolate integrated functions
    gsl_spline2d_init(spline_LogG_1, a.data(), b.data(), z_G_1, a.size(), b.size());
    gsl_spline2d_init(spline_LogG_2, a.data(), b.data(), z_G_2, a.size(), b.size());
    gsl_spline2d_init(spline_LogNH_2, a.data(), b.data(), z_NH_2, a.size(), b.size());
    gsl_spline2d_init(spline_LogI_1, a.data(), b.data(), z_I_1, a.size(), b.size());
    gsl_spline2d_init(spline_LogNI_2, a.data(), b.data(), z_NI_2, a.size(), b.size());

    OmB_Interpolators.LogG_1 = spline_LogG_1;
    OmB_Interpolators.LogG_1_acc_a = LogG_1_acc_a;
    OmB_Interpolators.LogG_1_acc_b = LogG_1_acc_b;
    OmB_Interpolators.LogG_2 = spline_LogG_2;
    OmB_Interpolators.LogG_2_acc_a = LogG_2_acc_a;
    OmB_Interpolators.LogG_2_acc_b = LogG_2_acc_b;
    OmB_Interpolators.LogNH_2 = spline_LogNH_2;
    OmB_Interpolators.LogNH_2_acc_a = LogNH_2_acc_a;
    OmB_Interpolators.LogNH_2_acc_b = LogNH_2_acc_b;
    OmB_Interpolators.LogI_1 = spline_LogI_1;
    OmB_Interpolators.LogI_1_acc_a = LogI_1_acc_a;
    OmB_Interpolators.LogI_1_acc_b = LogI_1_acc_b;
    OmB_Interpolators.LogNI_2 = spline_LogNI_2;
    OmB_Interpolators.LogNI_2_acc_a = LogNI_2_acc_a;
    OmB_Interpolators.LogNI_2_acc_b = LogNI_2_acc_b;

    free(z_G_1);
    free(z_G_2);
    free(z_NH_2);
    free(z_I_1);
    free(z_NI_2);

    /*
        Part 3: 2D interpolating functions for the two-argument functions

        h_1(a,b), h_2(a,b)

        where a = mu_e/M_e, b = eB/M_e^2

        h_1.dat and h_2.dat data are organized by row values 1 < a < 170 and column values 0 < b < 5000

        Used in high-temperature approximations to magnetization-related thermodynamic functions

        Load in 2D data tables h_1.dat and h_2.dat and then interpolate using GSL based on these data tables
    */

    std::ostringstream filenameh_1;
	filenameh_1 << "InterpolatingData/h_1.dat";
	std::fstream Datafileh_1;
    Datafileh_1.open(filenameh_1.str().c_str(),std::fstream::in);
    N_dat = 0; //number of rows in data table. Initialize to zero
	while (std::getline(Datafileh_1, line))
		++N_dat;
	Datafileh_1.close();

    //vectors to store values from data table
	a.clear(); //a and values. Same for both data sets. Use same vectors as before.
	b.clear();
    double h_1[N_dat-1][N_dat-1];
    double h_2[N_dat-1][N_dat-1];

    //temporary variables to hold data table values
    double h_1_dt, h_2_dt;
    Datafileh_1.open(filenameh_1.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                Datafileh_1 >> b_dt;
                b.push_back(b_dt);
            }
        }
        else{
            Datafileh_1 >> a_dt;
            a.push_back(a_dt);
            for(size_t j=0; j<N_dat-1; j++){
                Datafileh_1 >> h_1_dt;
                h_1[i-1][j] = h_1_dt;
            }
        }
	}
	Datafileh_1.close();

    std::ostringstream filenameh_2;
	filenameh_2 << "InterpolatingData/h_2.dat";
	std::fstream Datafileh_2;
    Datafileh_2.open(filenameh_2.str().c_str(),std::fstream::in);
	for(size_t i=0; i<N_dat; i++){
        if(i == 0){ //for first row, which consists of b values, simply fill b_Vec
            for(size_t j=0; j<N_dat-1; j++){
                Datafileh_2 >> b_dt;
            }
        }
        else{
            Datafileh_2 >> a_dt;
            for(size_t j=0; j<N_dat-1; j++){
                Datafileh_2 >> h_2_dt;
                h_2[i-1][j] = h_2_dt;
            }
        }
	}
	Datafileh_2.close();

    //Declare interpolating functions
    const gsl_interp2d_type *T_h_1 = gsl_interp2d_bicubic;
    const gsl_interp2d_type *T_h_2 = gsl_interp2d_bicubic;

    double *z_h_1 = (double *) malloc(a.size() * b.size() * sizeof(double));
    double *z_h_2 = (double *) malloc(a.size() * b.size() * sizeof(double));

    gsl_interp_accel *h_1_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *h_1_acc_b = gsl_interp_accel_alloc();
    gsl_interp_accel *h_2_acc_a = gsl_interp_accel_alloc();
    gsl_interp_accel *h_2_acc_b = gsl_interp_accel_alloc();

    gsl_spline2d *spline_h_1 = gsl_spline2d_alloc(T_h_1, a.size(), b.size());
    gsl_spline2d *spline_h_2 = gsl_spline2d_alloc(T_h_2, a.size(), b.size());

    //Integrate functions. Interpolation performed on integrals
    for(size_t i=0; i<a.size(); i++){
        for(size_t j=0; j<b.size(); j++){

            gsl_spline2d_set(spline_h_1, z_h_1, i, j, h_1[i][j] );
            gsl_spline2d_set(spline_h_2, z_h_2, i, j, h_2[i][j] );

        }
    }

    //Interpolate integrated functions
    gsl_spline2d_init(spline_h_1, a.data(), b.data(), z_h_1, a.size(), b.size());
    gsl_spline2d_init(spline_h_2, a.data(), b.data(), z_h_2, a.size(), b.size());

    OmB_Interpolators.h_1 = spline_h_1;
    OmB_Interpolators.h_1_acc_a = h_1_acc_a;
    OmB_Interpolators.h_1_acc_b = h_1_acc_b;
    OmB_Interpolators.h_2 = spline_h_2;
    OmB_Interpolators.h_2_acc_a = h_2_acc_a;
    OmB_Interpolators.h_2_acc_b = h_2_acc_b;

    free(z_h_1);
    free(z_h_2);

    /*
        Part 4: 1D interpolating functions for the two-argument functions

        i_1(x), i_2(x)

        where x = (mu_e^2-M_e^2)/eB

        Used in high-temperature approximations to magnetization-related thermodynamic functions

        Perform the required integrals in this function and then interpolate using GSL based on these data tables
    */
        //logarithmically-spaced values of a and b

    size_t num_interp = 200;
    double x_min = 1e-3, x_max = 1e4;
    double x_array[num_interp];
    double i_1_ar[num_interp], i_2_ar[num_interp], i_3_ar[num_interp]; //arrays to hold values of integrals and later to initialize interpolators

    struct HighTInterp_params x_param;

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (2000);

    double result, result2, error;

    gsl_function i_1;
    i_1.function = &i_1Integrand;
    i_1.params = &x_param;

    gsl_function i_2;
    i_2.function = &i_2Integrand;
    i_2.params = &x_param;

    gsl_function i_3;
    i_3.function = &i_3Integrand;
    i_3.params = &x_param;

    double y_low = 1e-8, y_high = 0.2;
    double pts[3] = {0.,0.,y_low};
    size_t numpts = 3;
    double zeta3_2 = 2.61237534869; //zeta(3/2)

    for(size_t i = 0; i<num_interp; i++){

        x_array[i] = pow( 10., log10(x_min)+(log10(x_max)-log10(x_min))/(double(num_interp-1))*double(i) );
        x_param = { x_array[i] };

        gsl_integration_qagp(&i_1, pts, numpts, 0, 1e-8, 2000, w, &result, &error);
        gsl_integration_qagiu(&i_1, y_low, 0, 1e-8, 2000, w, &result2, &error);
        i_1_ar[i] = result + result2 + sqrt(2./pi)*zeta3_2;

        gsl_integration_qag(&i_2, 0., y_high, 0, 1e-8, 2000, 3, w, &result, &error);
        gsl_integration_qagiu(&i_2, y_high, 0, 1e-8, 2000, w, &result2, &error);
        i_2_ar[i] = result + result2;

        gsl_integration_qag(&i_3, 0., y_high, 0, 1e-8, 2000, 3, w, &result, &error);
        gsl_integration_qagiu(&i_3, y_high, 0, 1e-8, 2000, w, &result2, &error);
        i_3_ar[i] = result + result2;
    }

    gsl_interp_accel *i_1_acc = gsl_interp_accel_alloc ();
	gsl_spline *spline_i_1 = gsl_spline_alloc (gsl_interp_steffen, num_interp);
	gsl_spline_init(spline_i_1, x_array, i_1_ar, num_interp);
    gsl_interp_accel *i_2_acc = gsl_interp_accel_alloc ();
	gsl_spline *spline_i_2 = gsl_spline_alloc (gsl_interp_steffen, num_interp);
	gsl_spline_init(spline_i_2, x_array, i_2_ar, num_interp);
	gsl_interp_accel *i_3_acc = gsl_interp_accel_alloc ();
	gsl_spline *spline_i_3 = gsl_spline_alloc (gsl_interp_steffen, num_interp);
	gsl_spline_init(spline_i_3, x_array, i_3_ar, num_interp);

    OmB_Interpolators.i_1 = spline_i_1;
    OmB_Interpolators.i_1_acc = i_1_acc;
    OmB_Interpolators.i_2 = spline_i_2;
    OmB_Interpolators.i_2_acc = i_2_acc;
    OmB_Interpolators.i_3 = spline_i_3;
    OmB_Interpolators.i_3_acc = i_3_acc;

    gsl_integration_workspace_free (w);

    return;
}

/*
    Integrand for function i_1(x) used in high-T approximation for magnetization and its derivatives
    Input: y: quantity being integrated over
           params: parameters for integrand (x)
*/
double i_1Integrand(double y, void * params)
{

    struct HighTInterp_params *p = (struct HighTInterp_params *) params;
    double x = p -> x;

    double result;

    if( y < 1e-4 ){
        result = sqrt(y)/3.*( -x + x*x*y/2. + x*( 1./15. - x*x/6. )*y*y );
    }
    else if( y > 40. ){
        result = ( y - 1. )/pow(y,2.5)*( exp(-y*x) - 1. );
    }
    else{
        result = ( y/tanh(y) - 1. )/pow(y,2.5)*( exp(-y*x) - 1. );
    }

    return result;
}

/*
    Integrand for function i_2(x) used in high-T approximation for magnetization and its derivatives
    Input: y: quantity being integrated over
           params: parameters for integrand (x)
*/
double i_2Integrand(double y, void * params)
{

    struct HighTInterp_params *p = (struct HighTInterp_params *) params;
    double x = p -> x;

    double result;

    if( y < 1e-4 ){
        result = sqrt(y)/3.*( 1. - y*y/15. + 2.*y*y*y*y/315. )*exp(-y*x);
    }
    else if( y > 40. ){
        result = ( y - 1. )/pow(y,1.5)*exp(-y*x);
    }
    else{
        result = ( y/tanh(y) - 1. )/pow(y,1.5)*exp(-y*x);
    }

    return result;
}

/*
    Integrand for function i_3(x) used in high-T approximation for magnetization and its derivatives (only for varying T case)
    Input: y: quantity being integrated over
           params: parameters for integrand (x)
*/
double i_3Integrand(double y, void * params)
{
    struct HighTInterp_params *p = (struct HighTInterp_params *) params;
    double x = p -> x;

    double result;

    if( y < 1e-6 ){
        result = pow(y,1.5)/3.*( 1. - y*y/15. + 2.*y*y*y*y/315. - y*y*y*y*y*y/1575. )*exp(-y*x);
    }
    else if( y > 40. ){
        result = ( y - 1. )/sqrt(y)*exp(-y*x);
    }
    else{
        result = ( y/tanh(y) - 1. )/sqrt(y)*exp(-y*x);
    }

    return result;
}

/*
        Computes electron chemical potential mu_e profile, assuming B = 0.
        We assume that this remains fixed throughout the simulation, and that n_e adjusts
        due to the evolution of B and T
        Inputs: n_e: electron number density in fm^{-3}. A ScalarField object that can vary in two dimensions.
        Output: mu_e: electron chemical potential in MeV. A RadialScalarField object that only varies in the radial direction
*/
void Compute_mu_e(ScalarField & n_e, RadialScalarField & mu_e)
{
    for(size_t i=0; i<n_e.shape()[0]; i++){
        mu_e[i] = sqrt( pow( 3.*pi*pi*n_e[i][0],2./3.)*hbarc*hbarc + M_e*M_e );
    }

    return;
}

/*
        Anisotropic (magnetic field dependent) crust electrical and thermal conductivity in s^{-1}. Returns Ohmic diffusivity in reduced units and/or thermal conductivity in reduced units.
        Returns result with Landau quantization effects (Shubnikov-de Haas oscillations) if input parameter quantization == true. If this is true, then eta_O_perp and
        eta_O_parallel = eta_O_delta + eta_O_perp, the Ohmic diffusivities perpendicular and parallel to the magnetic field, will be different. Otherwise they will equal each other.
        The thermal conductivities kappa_perp, kappa_parallel = kappa_delta + kappa_perp and kappa_H will always be distinct for nonzero B, in both quantizing and non-quantizing cases as long
        as anisotropy conductivity is enabled.

        Uses results from Potekhin A&A, 351, 787 (1999), Potekhin, Baiko, Haensel and Yakovlev, A&A 346, 345 (1999),
        Gnedin et al. MNRAS 324, 725 (2001), and Schmitt and Shternin, "Reaction rates and transport in neutron stars"
        in The Physics and Astrophysics of Neutron Stars. Springer, Heidelberg, pg. 455-574 (2018)
        Use Pearson et al. MNRAS 481, 2994 (2018), Eq. (C21) and Table C10 for computing x_nuc=sqrt(3/5)*x_p for BSk24, fit to Douchin and Haensel, A&A 380, 151 (2001) data for SLy4.
        Note: Eq. (C21) has a typo in the denominator of the first term: 1 + p_4... should be 1 - p_4... (checked via their EOS code). Also replaced p_6 with 0.325 to improve fit
        For electron-electron scattering contribution (thermal conductivity) uses Shternin and Yakovlev PRD 74, 043004 (2006) with non-degenerate electron correction
        from Cassisi et al. ApJ, 661, 1094 (2007)
        For ion (phonon) contribution (thermal conductivity) uses Chugunov and Haensel MNRAS 381, 1143 (2007) with typo corrections based on Potekhin's conductivity code

        Inputs: T: temperature in reduced units
                Bmag: magnetic field magnitude in reduced units
                n_e: electron number density in fm^{-3}
                A: mass number of nuclei (excludes dripped neutrons and protons)
                Z: atomic number of nuclei (excludes dripped protons)
                n_i: ion number density in fm^{-3}
                n_b: baryon number density in fm^{-3}
                cparams: ConductParams object which conductivities to compute, which impurity model to use, whether to include quantization effects, neutron drip density
                tparams: TParams object containing information about core temperature.
                N_GC: number of ghost cells
                n_t: lapse function
                c_v: specific heat capacity in reduced units
        Output: eta_O_delta, eta_O_perp: ScalarField object of Ohmic diffusivity parallel to magnetic field minus that perpendicular to it/perpendicular to it in reduced units. eta_O_delta < 0, and tends to 0 as B decreases
                kappa_delta, kappa_perp, kappa_H: ScalarField object of thermal conductivity parallel to magnetic field minus that perpendicular to it/perpendicular to it/Hall thermal conductivity in reduced units
                max_eta_T: maximum value of thermal diffusivity across simulation domain (reduced units)
*/
void ConductivityCalcB(ScalarField & T, ScalarField & Bmag, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC,
    std::vector<double> & n_t, ScalarField & c_v, ScalarField & eta_O_delta, ScalarField & eta_O_perp, ScalarField & kappa_delta, ScalarField & kappa_perp, ScalarField & kappa_H, double & max_eta_T, const Process & process, const Domain & dm)
{
    static double un1 = 2.78; //n = -1 frequency moment of bcc Coulombic lattice
    static double un2 = 12.973; //n = -2 frequency moment of bcc Coulombic lattice
    static double hbar = 6.582119569e-22; //reduced Planck constant in MeV*s
    static double kappa_pref = 1e39*pi*pi*k_Bcgs*c*c; //constant prefactor for computing thermal conductivities. 1e39 factor converts n_e to cm^{-3} and MeVtoErg/(c*c) converts mu_e to erg/c^2 = g
    static double T_plasma_pref = 4.*pi*alpha_e*hbarc*hbarc*hbarc; //constant prefactor used in computing electron plasma temperature

    static double eB_pref = unit_e*B_0*statCGtoMeV2; //constant prefactor for computing e*B in MeV^2

    static double n_drip = cparams.n_b_nd; //neutron drip density in fm^{-3} from equation of state
    double T_MeV, k_F, x_r, mu_e, a_i, q_BZ, Gamma, T_plasma, eta, eta_0, Z_cell;
    double v_F, beta, q_D, q_i2, k_TF2, w0, x_nuc, x_nuct, tp, q_s, s, w1, D;
    double G0, G2, G_s_sigma, G_s_kappa, w, sw, sw1, svF2, expi_sw, expi_sww, expi_sw1, expi_sw1w1, Lambda_1 = 0., Lambda_2 = 0., Lambda0_1 = 0., Lambda0_2 = 0.;
    double eB, nu, n_max, RegEps, omega_g, b, x, x2, E, Etilde, atilde, L, Ltilde;
    double A1, B1, C1, D1, A2, B2, C2, D2;
    double CoulombLog0, CoulombLogPar, CoulombLogPerp, PSCorrection, tau_par_ei, tau_perp_ei, tau_par_e, tau_perp_e;
    double a_m, xi, xi_s, u_0, zeta, Qparallel, Qperpxi, Qperpxi0, Phi_eps, Psi_eps;
    double tau_ee, T_plasma_e, th, I_l, I_t, I_lt, A_fact, Bpow, C, C_1, C_2, T_TF, t, NonDegenerateFactor_ee; //for electron-electron scattering thermal conductivity
    double kappa_i, kappa_ii, kappa_ie, LogFact, kappa_0, c_s, F_th, w_DW, w_form, y, exp_nw, exp_nwy2, Lambda_phe, Arho6Ap, L_ph; //for ionic thermal conductivity
    double sigma_par, sigma_perp; //electrical conductivity parallel/perpendicular to B in s^{-1}
    double lgnmax, sqrtx, sqrtb;

    double Q = 0., s_imp, w1_imp, s_impvF2, Lambda_1_imp, Lambda_2_imp, Lambda0_1_imp, Lambda0_2_imp, CoulombLog_imp, tau_ei_imp;

    // parameters for specific heat capacity of a bcc lattice
    static double alpha[5] = {0.932446, 0.334547, 0.265764, 4.757014e-3, 4.7770935e-3};
    static double an[9] = {1., 0.1839, 0.593586, 5.4814e-3, 5.01813e-4, 0, 3.2947e-7, 0, 5.8356e-11};
    static double bn[8] = {261.66, 0, 7.07997, 0, 0.0409484, 3.97355e-4, 5.11148e-5, 2.19749e-6};
    double C_vi, th2, th3, th4, th5, th6, th7, th8, th9, th10, th11, Ath, Bth, Athp, Athpp, Bthp, Bthpp;

    for(size_t i=1; i<eta_O_delta.shape()[0]-N_GC+1; i++){

        T_plasma = sqrt(756.40518016759142*Z[i]*Z[i]*n_i[i]/A[i]); //plasma temperature in MeV. 756.40518016759142 = 4*pi*alpha_e*(hbar*c)^3/(1 amu*c^2)
        eta_0 = 0.19/pow(Z[i],1./6.);
        a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm
        q_BZ = pow( 6.*pi*pi*n_i[i],1./3. ); //Brillouin zone boundary wavenumber in fm^{-1}

        //Set Q = Z_imp^2 in A. Potekhin's code
        if(cparams.Z_impurity == "Carreau2020BSk24"){
            Q = Q_impCarreau2020BSk24(n_b[i],n_drip);
        }
        else if(cparams.Z_impurity == "Vigano2013Q100"){
            Q = Q_impVigano2013Q100(n_b[i]);
        }
        else{
            Q = 0.;
        }

        for(size_t j=0; j<eta_O_delta.shape()[1]; j++){

            if(i < N_GC){
                T_MeV = T[N_GC+1-i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            else{
                T_MeV = T[i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            Gamma = alpha_e*hbarc*Z[i]*Z[i]/(T_MeV*a_i); //plasma coupling parameter (dimensionless)
            eta = T_MeV/T_plasma; //ratio of temperature to plasma temperature
            k_F = pow(3.*pi*pi*n_e[i][j],1./3.)*hbarc; //Fermi wavenumber = Fermi momentum in MeV
            x_r = k_F/M_e; //also used for p_0bar in some formulae
            mu_e = M_e*sqrt(x_r*x_r + 1.); //electron chemical potential in MeV
            v_F = sqrt( 1. - M_e*M_e/(mu_e*mu_e) ); //Fermi velocity in units of c
            beta = pi*alpha_e*Z[i]*v_F; //dimensionless beta parameter. beta_Z in Gnedin et al. 2001
            q_D = sqrt(3.*Gamma)*hbarc/a_i; //inverse Debye screening length in MeV
            q_i2 = q_D*q_D*( 1. + 0.06*Gamma )*exp(-sqrt(Gamma)); //q_i^2 in MeV^2
            k_TF2 = 4.*alpha_e/pi*k_F*mu_e; //Thomas-Fermi wave number squared in MeV^2

            Z_cell = n_e[i][j]/n_i[i];
            if(n_b[i] < n_drip) x_nuc = 1.15*pow(A[i],1./3.)/a_i; //Ion radius divided by Wigner-Seitz cell radius. From Kaminker et al., Astron. Astrophys. 343, 1009 (1999)
            //From Pearson et al. MNRAS 481, 2994 (2018) for BSk24 EOS. Eq. (C21) and Table C10. Use Zeq=Zcell, not Z.
            else if(cparams.EOS == "BSk24") x_nuc = ( 0.1035 + 1.944*pow(n_b[i],0.5717) )/( 1. - 608.*pow(n_b[i],3.143) ) + 0.0225*Z_cell*pow(n_b[i],1.26);
            else if(cparams.EOS == "SLy4") x_nuc = ( 0.08633 + 1.0173*pow(n_b[i],0.4857) )/( 1. - 237.6*pow(n_b[i],2.615) ) + 0.2732*Z_cell*pow(n_b[i],0.7222); //Pearson et al. 2018 formula fitted to Douchin and Haensel 2001 data
            else x_nuc = ( 0.1035 + 1.944*pow(n_b[i],0.5717) )/( 1. - 608.*pow(n_b[i],3.143) ) + 0.0225*Z_cell*pow(n_b[i],1.26); //Uses BSk24 formula by default
            tp = 1./( 0.1558 + 0.225*n_b[i] + 9.452*n_b[i]*n_b[i] ); //from CONDBSk subroutine of "condBSk.f" code of A. Potekhin (https://www.ioffe.ru/astro/conduct/condBSk.f)
            x_nuct = x_nuc*tp/( 0.6 + tp ); //x_nuc with thermal correction from

            q_s = sqrt( q_i2 + k_TF2 )*exp(-0.5*beta); //screening length in MeV
            s = q_s*q_s/(4.*k_F*k_F); // s = q_s^2/(2*k_F)^2. Also called the Coulomb screening parameter a_s (dimensionless)
            w0 = un2*4.*k_F*k_F/(q_D*q_D)*( 1. + beta/3. ); //w(q=2*k_F); called "w" in e.g. Potekhin et al. 1999. Also called a_DW, the Debye-Waller factor, in Potekhin 1999.
            w1 = 14.7327*x_nuc*x_nuc*( 1. + sqrt(x_nuc)*Z[i]/13. )*( 1. + beta/3. ); //finite ion size correction factor to w0
            D = exp(-0.42075*sqrt(x_r/(A[i]*Z[i]))*un1*exp(-9.1*eta));

            //Ionic quantum correction factors
            G0 = eta/sqrt(eta*eta + eta_0*eta_0)*( 1. + 0.122*beta*beta );
            G2 = 0.0105*( 1. - 1./Z[i] )*( 1. + v_F*v_F*v_F*beta )*eta/pow( eta*eta + 0.0081,1.5 )*( 1. + x_nuct*x_nuct*sqrt(2.*Z[i]) ); //1. + pow(x_nuct,2)*sqrt(2.*Z[i]) factor is from Gnedin et al. 2001.

            w = w0 + w1; //finite size ion-corrected Debye-Waller factor
            sw = s*w; //s*w
            sw1 = s*w1; //s*w1
            svF2 = v_F*v_F*s; //s times (v_F/c)^2 = s*(1-(M_e/mu_e^2)

            //For impurity scattering
            if(cparams.Z_impurity != "zero"){
                s_imp = k_TF2/(4.*k_F*k_F);
                w1_imp = 14.7327*x_nuc*x_nuc;
                s_impvF2 = v_F*v_F*s_imp;
                Lambda_1_imp = 0.5*( std::log1p(1./s_imp) - 1./(1.+s_imp ) );
                Lambda_2_imp = v_F*v_F*( (2.*s_imp +1.)/(2.*s_imp +2.) - s_imp*std::log1p(1./s_imp ) );
                Lambda0_1_imp = ( std::log1p(1./s_imp) + 1./(1.+1./s_imp)*(1.-exp(-w1_imp)) - (1.+s_imp*w1_imp)*( expE1_rp(s_imp*w1_imp)-exp(-w1_imp)*expE1_rp(s_imp*w1_imp+w1_imp) ) )/2.;
                Lambda0_2_imp = ( v_F*v_F*(exp(-w1_imp)-1.+w1_imp)/w1_imp - s_impvF2/(1.+1./s_imp)*(1.-exp(-w1_imp)) - 2.*s_impvF2*std::log1p(1./s_imp) + s_impvF2*(2.+s_imp*w1_imp)*( expE1_rp(s_imp*w1_imp)-exp(-w1_imp)*expE1_rp(s_imp*w1_imp+w1_imp) ) )/2.;
                CoulombLog_imp = (Lambda_1_imp - Lambda_2_imp) - (Lambda0_1_imp - Lambda0_2_imp); //Coulomb logarithm for impurity scattering
                tau_ei_imp = 5.699361922e-17/( Q/Z[i]*CoulombLog_imp*sqrt(1. + x_r*x_r) ); //Impurity scattering relaxation time in s
            }
            else tau_ei_imp = 1e100; //if Q = 0 (no impurities), set impurity scattering relaxation time to be extremely large so it does not affect result

            //Note: gsl_sf_expint_Ei = -int_{-x}^{inf}dt*exp(-t)/t
            if( w < 1e-15 ){
                Lambda_1 = 0.5*w*( 2. - 1./(1.+s) - 2.*s*std::log1p(1./s) );
                Lambda_2 = v_F*v_F*0.5*w*( 1.5 - 3.*s - 1./(1.+s) + 3.*s*s*std::log1p(1./s) );
            }
            else if( w > 100. ){
                Lambda_1 = 0.5*( std::log1p(1./s) - 1./(1.+s) - 1./(sw*sw) );
                Lambda_2 = v_F*v_F*( (2.*s+1.)/(2.*s+2.) - s*std::log1p(1./s) );
            }
            else if(sw < 10.){
                expi_sw = -gsl_sf_expint_Ei(-sw);
                expi_sww = -gsl_sf_expint_Ei(-sw-w);
                Lambda_1 = ( std::log1p(1./s) + 1./(1.+1./s)*-std::expm1(-w) - (1.+sw)*exp(sw)*( expi_sw - expi_sww ) )/2.;
                Lambda_2 = ( v_F*v_F*(std::expm1(-w)+w)/w - svF2/(1.+1./s)*-std::expm1(-w)
                            - 2.*svF2*std::log1p(1./s) + svF2*(2.+sw)*exp(sw)*( expi_sw - expi_sww ) )/2.;
            }
            else{
                Lambda_1 = ( std::log1p(1./s) + 1./(1.+1./s)*-std::expm1(-w) - (1.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
                Lambda_2 = ( v_F*v_F*(std::expm1(-w)+w)/w - svF2/(1.+1./s)*-std::expm1(-w)
                            - 2.*svF2*std::log1p(1./s) + svF2*(2.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
            }
            if( w1 < 1e-15 ){
                Lambda0_1 = 0.5*w1*( 2. - 1./(1.+s) - 2.*s*std::log1p(1./s) );
                Lambda0_2 = v_F*v_F*0.5*w1*( 1.5 - 3.*s - 1./(1.+s) + 3.*s*s*std::log1p(1./s) );
            }
            else if( w1 > 100. ){
                Lambda0_1 = 0.5*( std::log1p(1./s) - 1./(1.+s) - 1./(sw1*sw1) );
                Lambda0_2 = v_F*v_F*( (2.*s+1.)/(2.*s+2.) - s*std::log1p(1./s) );
            }
            else if(sw1 < 10.){
                expi_sw1 = -gsl_sf_expint_Ei(-sw1);
                expi_sw1w1 = -gsl_sf_expint_Ei(-sw1-w1);
                Lambda0_1 = ( std::log1p(1./s) + 1./(1.+1./s)*-std::expm1(-w1) - (1.+sw1)*exp(sw1)*( expi_sw1 - expi_sw1w1 ) )/2.;
                Lambda0_2 = ( v_F*v_F*(std::expm1(-w1)+w1)/w1 - svF2/(1.+1./s)*-std::expm1(-w1) - 2.*svF2*std::log1p(1./s) + svF2*(2.+sw1)*exp(sw1)*( expi_sw1 - expi_sw1w1 ) )/2.;
            }
            else{
                Lambda0_1 = ( std::log1p(1./s) + 1./(1.+1./s)*-std::expm1(-w1) - (1.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
                Lambda0_2 = ( v_F*v_F*(std::expm1(-w1)+w1)/w1 - svF2/(1.+1./s)*-std::expm1(-w1) - 2.*svF2*std::log1p(1./s) + svF2*(2.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
            }

            //Quantities required to compute the anisotropic, B-dependent transport coefficients
            eB = eB_pref*Bmag[i][j]; //e*B in MeV^2
            nu = k_F*k_F/(2.*eB); //highest occupied Landau level is integer part of this
            // n_max = std::floor(nu); //highest occupied Landau level at T=0
            double eps = 0.1;
            n_max = std::max( 0., std::floor(nu) - 1. + 0.5*( 2. + tanh( (nu-std::floor(nu))/eps ) - tanh( (std::ceil(nu)-nu)/eps ) ) ); //smoothed occupied Landau level at T=0
            if(n_max < 1.) RegEps = 0.; //regulates divergence in PSCorrection when nu is exactly an integer
            else RegEps = 0.5;
            omega_g = eB/(mu_e*hbar); //electron gyrofrequency in s^{-1}

            b = Bmag[i][j]*B_0/B_crit; //magnetic field magnitude in units of the quantum critical field
            double f = (nu - n_max);
            x2 = 2.*( f - 0.5*( tanh( f/eps ) - tanh( (1.-f)/eps ) ) ); // x^2 = sqrt(nu-n_max). Vanishes at exact filling of a Landau level. To minimize amplification of roundoff error, we regularize this with a smoothing function
            x = sqrt(x2);

            E = -std::expm1(-w0)/w0; //(1. - exp(-w0))/w0;
            Etilde = 1./(10. + 5./b);
            atilde = ( sqrt(s) + 1./(2. + 0.5*w0) )*( sqrt(s) + 1./(2. + 0.5*w0) );
            L = std::log1p(1./atilde); //log(1.+1./atilde)
            Ltilde = atilde*L;
            A1 = (30. - 15.*E - (15. - 6.*E)*v_F*v_F)/(30. - 10.*E - (20. - 5.*E)*v_F*v_F);
            B1 = 1.5 - 0.5*E + 0.25*v_F*v_F/(1. - 2./3.*v_F*v_F);
            C1 = (1. - E + 0.75*v_F*v_F)/(1. + v_F*v_F);
            D1 = 1. + 0.06*L*L/(n_max*n_max+1e-15); //regularized
            A2 = 0.8*( 1. + Ltilde ) + 0.2*L;
            B2 = (0.68 - 1.3*Etilde)*pow( Ltilde,1./6. );
            C2 = 1.42 - Etilde + sqrt(Ltilde)/3.;
            D2 = (0.52 - Etilde)*pow( Ltilde,0.25 );

            lgnmax = log(n_max);
            sqrtb = sqrt(b);
            sqrtx = sqrt(x);
            CoulombLogPar = 1./x/sqrt( D1/pow( x + sqrtb/x_r*( A1 - B1*sqrtx*x + C1*x*(x - sqrtx)/n_max ),2. )
                                       + L*L*pow( (3.*x2 - 1.)/( 2.*n_max + 1.5*x2/((1. + 2.*b)*(1. + 2.*b)) ) + 0.07 + E/5.,2. ) ); //ratio of Coulomb logarithm for conductivity perpendicular to magnetic field divided by zero field Coulomb logarithm (CoulombLogarithm0)
            CoulombLogPerp = 1. + sqrtb/(x_r*x2)*( sqrtb*A2/x_r + B2*lgnmax*x - (C2 + D2*lgnmax)*sqrtx*x2 ); //ratio of Coulomb logarithm for conductivity parallel to magnetic field divided by zero field Coulomb logarithm (CoulombLogarithm0)

            // Thermal conductivity for electron-electron scattering. Note: u in Shternin and Yakovlev 2006 = v_F
            T_plasma_e = sqrt(T_plasma_pref*n_e[i][j]/mu_e); //electron plasma temperature in MeV
            th = 1.7320508075688772*T_plasma_e/T_MeV; //1.7320508075688772 = sqrt(3)
            I_l = 1./v_F*( 0.1587 - 0.02538/( 1. + 0.0435*th ) )*log( 1. + 128.56/( 37.1*th + 10.83*th*th + th*th*th ) );
            A_fact = 20. + 450.*v_F*v_F*v_F;
            C_1 = 0.05607 + 0.03216*v_F*v_F;
            C_2 = 0.0254 + 0.04127*v_F*v_F*v_F*v_F;
            C = A_fact*exp( C_1/C_2 );
            I_t = v_F*v_F*v_F*( 2.404/C + ( C_2 - 2.404/C )/( 1. + 0.1*th*v_F ) )*log( 1. + C/( A_fact*th*v_F + th*th*v_F*v_F ) );
            A_fact = 12.2 + 25.2*v_F*v_F*v_F;
            Bpow = 1. - 0.75*v_F;
            C_1 = 0.123636 + 0.016234*v_F*v_F;
            C_2 = 0.0762 + 0.05714*v_F*v_F*v_F*v_F;
            C = A_fact*exp( C_1/C_2 );
            I_lt = v_F*( 18.52*v_F*v_F/C + ( C_2 - 18.2*v_F*v_F/C )/( 1. + 0.1558*pow( th,Bpow ) ) )*log( 1. + C/( A_fact*th + 10.83*th*th*v_F*v_F + pow( th*v_F,8./3. ) ) );
            T_TF = T_MeV/( mu_e-M_e ); //temperature to Fermi temperature ratio (dimensionless). Called theta in Cassisi et al. 2007
            t = 25.*T_TF;
            NonDegenerateFactor_ee = ( 1. + t*t )/( 1. + t + 0.4342*t*t*sqrt(T_TF) ); //multiply nu_ee by this factor to correctly interpolate to non-degenerate electrons
            tau_ee = pi*T_MeV*mu_e/( 36.*alpha_e*alpha_e*( I_l + I_t + I_lt )*NonDegenerateFactor_ee*n_e[i][j]*hbarc*hbarc*hbarc )*hbar; //electron-electron scattering time in s.

            // Ionic thermal conductivity. First computes the specific heat capacity per ion in units of k_B for the harmonic Coulomb lattice.
            th = T_plasma/T_MeV; // = eta^-1
            th2 = th*th;
            th3 = th2*th;
            th4 = th3*th;
            th5 = th4*th;
            th6 = th5*th;
            th7 = th6*th;
            th8 = th7*th;
            th9 = th8*th;
            th10 = th9*th;
            th11 = th10*th;
            Ath = an[0] + an[1]*th + an[2]*th2 + an[3]*th3 + an[4]*th4 + an[6]*th6 + an[8]*th8;
            Bth = bn[0] + bn[2]*th2 + bn[4]*th4 + bn[5]*th5 + bn[6]*th6 + bn[7]*th7 + alpha[3]*an[5]*th9 + alpha[4]*an[7]*th11;
            Athp = an[1] + 2.*an[2]*th + 3.*an[3]*th2 + 4.*an[4]*th3 + 6.*an[6]*th5 + 8.*an[8]*th7;
            Athpp = 2.*an[2] + 6.*an[3]*th + 12.*an[4]*th2 + 30.*an[6]*th4 + 56.*an[8]*th6;
            Bthp = 2.*bn[2]*th + 4.*bn[4]*th3 + 5.*bn[5]*th4 + 6.*bn[6]*th5 + 7.*bn[7]*th6 + 9.*alpha[3]*an[5]*th8 + 11.*alpha[4]*an[7]*th10;
            Bthpp = 2.*bn[2] + 12.*bn[4]*th2 + 20.*bn[5]*th3 + 30.*bn[6]*th4 + 42.*bn[7]*th5 + 72.*alpha[3]*an[5]*th7 + 110.*alpha[4]*an[7]*th9;
            C_vi = th2*( pow( alpha[0]/(exp(0.5*alpha[0]*th)-exp(-0.5*alpha[0]*th)),2. ) + pow( alpha[1]/(exp(0.5*alpha[1]*th)-exp(-0.5*alpha[1]*th)),2. )
                                     + pow( alpha[2]/(exp(0.5*alpha[2]*th)-exp(-0.5*alpha[2]*th)),2. ) + (Athpp*Bth*Bth-2.*Athp*Bthp*Bth+2.*Ath*Bthp*Bthp-Ath*Bth*Bthpp)/(Bth*Bth*Bth) ); //specific heat per ion divided by k_B (dimensionless)

            LogFact = log( 2. + 1./(sqrt(3.*Gamma)*Gamma) );
            kappa_0 = k_Bcgs*T_plasma*n_i[i]*a_i*a_i/hbarc*c*1e26; //normalization constant in units of erg/(cm*K*s). Numerical factor converts fm^-2 to cm^-2
            kappa_ii = kappa_0*sqrt( 16./(Gamma*Gamma*Gamma*Gamma*Gamma*LogFact*LogFact) + 0.16 + Gamma*Gamma/5929.*exp( 2./3.*th )  ); //ion-ion scattering part of ionic thermal conductivity in erg/(cm*K*s)

            c_s = T_plasma/( 3.*q_BZ*hbarc )*c; //sound speed in cm/s
            F_th = 0.014 + 0.03/( exp( th/5. ) + 1. );
            w_DW = 1.683*sqrt( x_r/(A[i]*Z[i]) )*( 0.5*un1*exp( -9.1/th ) + un2/th );
            w_form = 43.*x_nuc*x_nuc;
            w = w_DW + w_form;
            y = q_BZ*hbarc/( 2.*k_F ); //dimensionless
            exp_nw = exp(-w);
            exp_nwy2 = exp(-w*y*y);
            if( y > 1e-12 ){
                Lambda_phe = 0.5*( exp_nwy2*expE1_rp(w*y*y) - exp_nw*expE1_rp(w) - v_F*v_F*( exp_nwy2 - exp_nw )/w );
            }
            else{
                Lambda_phe = log( 1./y ) - v_F*v_F*( 1. - y*y )/2.;
            }
            Arho6Ap = 1.66053906892e-24*A[i]*n_i[i]*1e39/1e6;// A*rho_6/A'. Approximate rho = m_u*n_b where m_u = 1.66053906892e-24 g is the atomic mass unit and A' = n_b/n_i
            L_ph = 320.*a_i/( 1. + x_r*x_r )/Lambda_phe*26./Z[i]*F_th/0.01*sqrt( Arho6Ap ); //in fm
            kappa_ie = 1./3.*k_Bcgs*C_vi*n_i[i]*c_s*L_ph*1e26; //ion-electron scattering part of ionic thermal conductivity in erg/(cm*K*s). Numerical factor converts fm^-2 to cm^-2

            kappa_i = 1./( 1./kappa_ie + 1./kappa_ii );

            //Different values of G_s used for electrical and/or thermal conductivity.
            if(cparams.EorT == "Electrical"){
                G_s_sigma = G0; //G_sigma
                CoulombLog0 = ( (Lambda_1 - Lambda_2) - (Lambda0_1-Lambda0_2) )*G_s_sigma*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                PSCorrection = 2.*k_F*k_F*k_F/(3.*eB)/( sqrt( k_F*k_F - 2.*eB*n_max ) - 2./(3.*eB)*( pow( k_F*k_F - 2.*eB*n_max,1.5 ) - k_F*k_F*k_F ) - eB/6.*( 1./sqrt( k_F*k_F - 2.*eB*n_max+ RegEps ) - 1./k_F ) ); //phase space correction factor: (N_0(mu_e)/N_B(mu_e))^2

                if( Bmag[i][j] < 1e-10 || cparams.quantization == false ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLog0*sqrt(1. + x_r*x_r)); //electron-ion relaxation time in s. If magnetic field is extremely weak, simply use B=0 result
                    tau_par_ei = tau_perp_ei;
                }
                else if( nu > 1. ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLogPerp*CoulombLog0*sqrt(1. + x_r*x_r)); //perpendicular to magnetic field relaxation time in s for nu > 1
                    tau_par_ei = 5.699361922e-17*PSCorrection/(Z[i]*CoulombLogPar*CoulombLog0*sqrt(1. + x_r*x_r)); //parallel to magnetic field relaxation time in s for nu > 1
                }
                else{
                    //Uses Eq. (A8-A10) of Potekhin A&A, 351, 787 (1999)
                    a_m = 1./sqrt(eB); //magnetic quantum length in MeV^{-1}
                    xi = 2.*x_r*x_r/b;
                    xi_s = 0.5*a_m*a_m*q_s*q_s;
                    u_0 = xi + xi_s;
                    zeta = w*0.5/(a_m*a_m*k_F*k_F);
                    Qparallel = (1.-exp(-zeta*xi))/u_0 - expE1_rp(u_0) + (1. + zeta)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi = (1. + u_0)*expE1_rp(u_0) - 1. + exp(-zeta*xi) - (1. + u_0 + zeta*u_0)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi0 = (1. + xi_s)*expE1_rp(xi_s) - (1. + xi_s + zeta*xi_s)*expE1_rp( xi_s*(1. + zeta) );
                    Phi_eps = x_r*x_r/(2.*Qparallel);
                    Psi_eps = b/(x_r*x_r)*( mu_e*mu_e/(M_e*M_e)*Qperpxi + Qperpxi0 );

                    tau_perp_ei = (eB*mu_e*k_F)/(pi*(Z[i]*Z[i]*alpha_e*alpha_e)*hbarc*hbarc*hbarc*n_i[i]*M_e*M_e*Psi_eps)*hbar; //perpendicular to magnetic field relaxation time in s for nu<1
                    tau_par_ei = Phi_eps*mu_e*eB/(2.*pi*k_F*(Z[i]*Z[i]*alpha_e*alpha_e)*hbarc*hbarc*hbarc*n_i[i])*hbar; //parallel to magnetic field relaxation time in s for nu<1
                }

                //Includes impurity scattering in electron-ion relaxation times
                tau_perp_ei = tau_perp_ei/( 1. + tau_perp_ei/tau_ei_imp );
                tau_par_ei = tau_par_ei/( 1. + tau_par_ei/tau_ei_imp );

                sigma_par = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_par_ei; //returns parallel to B electrical conductivity in s^{-1}. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                sigma_perp = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_perp_ei; //returns perpendicular to B electrical conductivity in s^{-1}, excluding magnetic field-dependent factor, which cancels when inverting conductivity tensor.
                eta_O_perp[i][j] = c*c*t_0/(4.*pi*sigma_perp*L_0*L_0); //Ohmic diffusivity perpendicular to B in reduced units
                eta_O_delta[i][j] = c*c*t_0/(4.*pi*L_0*L_0)*(sigma_perp-sigma_par)/(sigma_par*sigma_perp); //Ohmic diffusivity parallel to B minus that perpendicular to B in reduced units

            }
            else if(cparams.EorT == "Thermal"){
                G_s_kappa = G0 + G2; //G_kappa.
                CoulombLog0 = ( (Lambda_1 - Lambda_2) - (Lambda0_1-Lambda0_2) )*G_s_kappa*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                PSCorrection = 2.*k_F*k_F*k_F/(3.*eB)/( sqrt( k_F*k_F - 2.*eB*n_max ) - 2./(3.*eB)*( pow( k_F*k_F - 2.*eB*n_max,1.5 ) - k_F*k_F*k_F ) - eB/6.*( 1./sqrt( k_F*k_F - 2.*eB*n_max+ RegEps ) - 1./k_F ) ); //phase space correction factor: (N_0(mu_e)/N_B(mu_e))^2

                if( Bmag[i][j] < 1e-10 || cparams.quantization == false ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLog0*sqrt(1. + x_r*x_r)); //if magnetic field is extremely weak or we ignore Landau quantization, simply use B=0 result
                    tau_par_ei = tau_perp_ei;
                }
                if( nu > 1. ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLogPerp*CoulombLog0*sqrt(1. + x_r*x_r)); //perpendicular to magnetic field relaxation time in s for nu > 1
                    tau_par_ei = 5.699361922e-17*PSCorrection/(Z[i]*CoulombLogPar*CoulombLog0*sqrt(1. + x_r*x_r)); //parallel to magnetic field relaxation time in s for nu > 1
                }
                else{
                    //Uses Eq. (A8-A10) of Potekhin A&A, 351, 787 (1999)
                    a_m = 1./sqrt(eB); //magnetic quantum length in MeV^{-1}
                    xi = 2.*x_r*x_r/b;
                    xi_s = 0.5*a_m*a_m*q_s*q_s;
                    u_0 = xi + xi_s;
                    zeta = w*0.5/(a_m*a_m*k_F*k_F);
                    Qparallel = (1.-exp(-zeta*xi))/u_0 - expE1_rp(u_0) + (1. + zeta)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi = (1. + u_0)*expE1_rp(u_0) - 1. + exp(-zeta*xi) - (1. + u_0 + zeta*u_0)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi0 = (1. + xi_s)*expE1_rp(xi_s) - (1. + xi_s + zeta*xi_s)*expE1_rp( xi_s*(1. + zeta) );
                    Phi_eps = x_r*x_r/(2.*Qparallel);
                    Psi_eps = b/(x_r*x_r)*( mu_e*mu_e/(M_e*M_e)*Qperpxi + Qperpxi0 );

                    tau_perp_ei = (eB*mu_e*k_F)/(pi*pow(Z[i]*alpha_e,2.)*hbarc*hbarc*hbarc*n_i[i]*M_e*M_e*Psi_eps)*hbar; //perpendicular to magnetic field relaxation time in s for nu<1
                    tau_par_ei = Phi_eps*mu_e*eB/(2.*pi*k_F*pow(Z[i]*alpha_e,2.)*hbarc*hbarc*hbarc*n_i[i])*hbar; //parallel to magnetic field relaxation time in s for nu<1
                }

                //Includes impurity scattering in electron-ion relaxation times
                tau_perp_ei = tau_perp_ei/( 1. + tau_perp_ei/tau_ei_imp );
                tau_par_ei = tau_par_ei/( 1. + tau_par_ei/tau_ei_imp );

                //Includes electron-electron scattering in relaxation time for thermal conductivity
                if( eB/mu_e < T_MeV && eB/mu_e < M_e ){
                    tau_perp_e = tau_perp_ei*tau_ee/( tau_perp_ei + tau_ee );
                }
                else tau_perp_e = tau_perp_ei;
                if( eB/mu_e < T_MeV || eB/mu_e < M_e ){
                    tau_par_e = tau_par_ei*tau_ee/( tau_par_ei + tau_ee );
                }
                else tau_par_e = tau_par_ei;

                tau_perp_e = tau_perp_ei*tau_ee/( tau_perp_ei + tau_ee );
                tau_par_e = tau_par_ei*tau_ee/( tau_par_ei + tau_ee );
                kappa_perp[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*( tau_perp_e/( 1. + pow(omega_g*tau_perp_e,2.) ) ) )*t_0/(s_0*L_0*L_0); //returns thermal conductivity perpendicular to B in reduced units.
                kappa_delta[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*tau_par_e )*t_0/(s_0*L_0*L_0) - kappa_perp[i][j]; //returns thermal conductivity parallel to B minus perpendicular to B in reduced units.
                kappa_H[i][j] = kappa_perp[i][j]*(omega_g*tau_perp_e); //returns thermal Hall conductivity (Leduc-Righi effect) in reduced units. Equals kappa_perp times the dimensionless omega_g*tau_perp
                kappa_perp[i][j] += kappa_i*t_0/(s_0*L_0*L_0); //adds ionic thermal conductivity to kappa_perp. Its contribution to kappa_delta will cancel, since it adds equally to kappa_perp and kappa_par. It does not contribute to kappa_H, so we avoid adding it until after computing kappa_H.

                if( (kappa_delta[i][j] + kappa_perp[i][j])/c_v[i][j] > max_eta_T ) max_eta_T = (kappa_delta[i][j] + kappa_perp[i][j])/c_v[i][j]; //computes maximum value of thermal diffusivity
            }
            else if(cparams.EorT == "Both"){
                G_s_sigma = G0; //G_sigma
                CoulombLog0 = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_sigma*D; //Fitted formula for Coulomb logarithm from Gnedin et al. 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                PSCorrection = 2.*k_F*k_F*k_F/(3.*eB)/( sqrt( k_F*k_F - 2.*eB*n_max ) - 2./(3.*eB)*( pow( k_F*k_F - 2.*eB*n_max,1.5 ) - k_F*k_F*k_F ) - eB/6.*( 1./sqrt( k_F*k_F - 2.*eB*n_max + RegEps ) - 1./k_F ) ); //phase space correction factor: (N_0(mu_e)/N_B(mu_e))^2

                if( Bmag[i][j] < 1e-6 || cparams.quantization == false ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLog0*sqrt(1. + x_r*x_r)); //electron-ion relaxation time in s. If magnetic field is extremely weak, simply use B=0 result
                    tau_par_ei = tau_perp_ei;
                }
                else if( nu > 1. ){
                    //Uses Eq. (42-45) of Potekhin A&A, 351, 787 (1999)
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLogPerp*CoulombLog0*sqrt(1. + x_r*x_r)); //perpendicular to B relaxation time in s for nu > 1
                    tau_par_ei = 5.699361922e-17*PSCorrection/(Z[i]*CoulombLogPar*CoulombLog0*sqrt(1. + x_r*x_r)); //parallel to B relaxation time in s for nu > 1
                }
                else{
                    //Uses Eq. (A8-A10) of Potekhin A&A, 351, 787 (1999)
                    a_m = 1./sqrt(eB); //magnetic quantum length in MeV^{-1}
                    xi = 2.*x_r*x_r/b;
                    xi_s = 0.5*a_m*a_m*q_s*q_s;
                    u_0 = xi + xi_s;
                    zeta = w*0.5/(a_m*a_m*k_F*k_F);
                    Qparallel = (1.-exp(-zeta*xi))/u_0 - expE1_rp(u_0) + (1. + zeta)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi = (1. + u_0)*expE1_rp(u_0) - 1. + exp(-zeta*xi) - (1. + u_0 + zeta*u_0)*expE1_rp( u_0*(1. + zeta) )*exp( zeta*(xi_s - u_0) );
                    Qperpxi0 = (1. + xi_s)*expE1_rp(xi_s) - (1. + xi_s + zeta*xi_s)*expE1_rp( xi_s*(1. + zeta) );
                    Phi_eps = x_r*x_r/(2.*Qparallel);
                    Psi_eps = b/(x_r*x_r)*( mu_e*mu_e/(M_e*M_e)*Qperpxi + Qperpxi0 );

                    tau_perp_ei = (eB*mu_e*k_F)/(pi*pow(Z[i]*alpha_e,2.)*hbarc*hbarc*hbarc*n_i[i]*M_e*M_e*Psi_eps)*hbar; //perpendicular to B relaxation time in s for nu<1
                    tau_par_ei = Phi_eps*mu_e*eB/(2.*pi*k_F*pow(Z[i]*alpha_e,2.)*hbarc*hbarc*hbarc*n_i[i])*hbar; //parallel to B relaxation time in s for nu<1
                }

                //Includes impurity scattering in electron-ion relaxation times
                tau_perp_ei = tau_perp_ei/( 1. + tau_perp_ei/tau_ei_imp );
                tau_par_ei = tau_par_ei/( 1. + tau_par_ei/tau_ei_imp );

                sigma_par = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_par_ei; //returns parallel to B electrical conductivity in s^{-1}. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                sigma_perp = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_perp_ei; //returns perpendicular to B electrical conductivity in s^{-1}, excluding magnetic field-dependent factor, which cancels when inverting conductivity tensor.
                eta_O_perp[i][j] = c*c*t_0/(4.*pi*sigma_perp*L_0*L_0); //Ohmic diffusivity perpendicular to B in reduced units
                // eta_O_delta[i][j] = c*c*t_0/(4.*pi*sigma_par*L_0*L_0) - eta_O_perp[i][j]; //Ohmic diffusivity parallel to B minus that perpendicular to B in reduced units
                eta_O_delta[i][j] = c*c*t_0/(4.*pi*L_0*L_0)*(sigma_perp-sigma_par)/(sigma_par*sigma_perp); //Ohmic diffusivity parallel to B minus that perpendicular to B in reduced units

                G_s_kappa = G0 + G2; //G_kappa. Add correction factor to G_sigma.
                CoulombLog0 = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_kappa*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2

                if( Bmag[i][j] < 1e-6 || cparams.quantization == false ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLog0*sqrt(1. + x_r*x_r)); //electron-ion relaxation time in s. If magnetic field is extremely weak, simply use B=0 result
                    tau_par_ei = tau_perp_ei;
                }
                else if( nu > 1. ){
                    tau_perp_ei = 5.699361922e-17/(Z[i]*CoulombLogPerp*CoulombLog0*sqrt(1. + x_r*x_r)); //relaxation time for scattering perpendicular to B in s for nu > 1
                    tau_par_ei = 5.699361922e-17*PSCorrection/(Z[i]*CoulombLogPar*CoulombLog0*sqrt(1. + x_r*x_r)); //relaxation time for scattering parallel to B in s for nu > 1
                }

                //Includes impurity scattering in electron-ion relaxation times
                tau_perp_ei = tau_perp_ei/( 1. + tau_perp_ei/tau_ei_imp );
                tau_par_ei = tau_par_ei/( 1. + tau_par_ei/tau_ei_imp );

                //Includes electron-electron scattering in relaxation time for thermal conductivity
                if( eB/mu_e < T_MeV && eB/mu_e < M_e ){
                    tau_perp_e = tau_perp_ei*tau_ee/( tau_perp_ei + tau_ee );
                }
                else tau_perp_e = tau_perp_ei;
                if( eB/mu_e < T_MeV || eB/mu_e < M_e ){
                    tau_par_e = tau_par_ei*tau_ee/( tau_par_ei + tau_ee );
                }
                else tau_par_e = tau_par_ei;

                kappa_perp[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*( tau_perp_e/( 1. + pow(omega_g*tau_perp_e,2.) ) ) )*t_0/(s_0*L_0*L_0); //returns thermal conductivity perpendicular to B in reduced units. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                kappa_delta[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*tau_par_e )*t_0/(s_0*L_0*L_0) - kappa_perp[i][j]; //returns thermal conductivity parallel to B minus perpendicular to B in reduced units. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                kappa_H[i][j] = kappa_perp[i][j]*(omega_g*tau_perp_e); //returns thermal Hall conductivity (Leduc-Righi effect) in reduced units. Equals kappa_perp times the dimensionless omega_g*tau_perp
                kappa_perp[i][j] += kappa_i*t_0/(s_0*L_0*L_0); //adds ionic thermal conductivity to kappa_perp. Its contribution to kappa_delta will cancel, since it adds equally to kappa_perp and kappa_par. It does not contribute to kappa_H, so we avoid adding it until after computing kappa_H.
                if( (kappa_delta[i][j] + kappa_perp[i][j])/c_v[i][j] > max_eta_T && i < eta_O_delta.shape()[0]-N_GC+1 ) max_eta_T = (kappa_delta[i][j] + kappa_perp[i][j])/c_v[i][j]; //computes maximum value of thermal diffusivity

            }
            else{
                std::cout << "Invalid value for EorT in sigmaCalcB" << std::endl;
            }
        }
    }

    return;
}

/*
    Computes specific heat capacity at nonzero B in reduced units. Includes electron, ion and free neutron contributions

    Thermal harmonic part of ionic heat capacity c_vi_h from Baiko, Potekhin and Yakovlev PRE 64, 057402 (2001), assuming bcc lattice
    Thermal anharmonic part of ionic heat capacity c_vi_ah from Baiko and Chugunov, MNRAS 510, 2628–2643 (2022)
    Electron polarization term c_vei from Potekhin and Chabrier Contrib. Plasma Phys. 50, 82 (2010) and Potekhin and Chabrier A&A, 550, A43 (2013) Appendix C.2
    Neutron heat capacity c_vn from Pastore, Chamel and Margueron, MNRAS 448, 1887 (2015)
    Neutron 1S0 pairing gap in crust parameterization taken from model of Kaminker, Haensel and Yakovlev A&A 373, L17 (2001)
        and Andersson, Comer and Glampedakis Nucl. Phys. A 763, 212 (2005)
    Specific pairing gap model taken from Schwenk, Friman and Brown Nucl. Phys. A 713, 191 (2003), SFB parametrization in Ho et al. PRC 91, 015806 (2015) Table II.
    Neutron 1S0 superfluidity control function R_A(T/T_cn) taken from Eq. (18), R_A of Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999)

    Inputs: T: temperature in reduced units
            n_e: electron number density in fm^{-3}
            A: mass number of nuclei (excludes dripped neutrons and protons)
            Z: atomic number of nuclei (excludes dripped protons)
            n_i: ion number density in fm^{-3}
            n_b: baryon number density in fm^{-3}
            n_nf: free neutron number density in fm^{-3}
            mu_nf: free neutron chemical potential in MeV
            tparams: TParams object containing information about nucleon superfluid gaps.
            N_GC: number of ghost cells
            n_t: lapse function
            c_v: ScalarField object containing specific heat capacity of magnetized electrons (B-dependent) in reduced units, computed by "Omega_xyVaryingT"
    Output: c_v: ScalarField object containing specific heat capacity of electrons, ions and free neutrons in reduced units
*/
void c_vfuncCalcB(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, RadialScalarField & n_nf, RadialScalarField & mu_nf, const TParams & tparams, size_t N_GC, std::vector<double> & n_t, ScalarField & c_v)
{
    // parameters for bcc lattice
    static double T_pconst = hbarc*hbarc*hbarc*4.*pi*alpha_e/931.49410372;// constant appearing in plasma temperature calculation where m_u = 931.49410372 MeV is the atomic mass unit
    static double alpha[5] = {0.932446, 0.334547, 0.265764, 4.757014e-3, 4.7770935e-3};
    static double an[9] = {1., 0.1839, 0.593586, 5.4814e-3, 5.01813e-4, 0, 3.9247e-7, 0, 5.8356e-11};
    static double bn[8] = {261.66, 0, 7.07997, 0, 0.0409484, 3.97355e-4, 5.11148e-5, 2.19749e-6};
    static double Acl[3] = {10.2, 248., 2.03e5};
    static double q = 0.205, ec = exp(1.);
    static double A11 = -10., A12 = 6e-3;
    static double A1q = -0.62/6., A2q = -0.56;
    static double A13 = -Acl[0] - A11;
    static double A14 = ( A1q - A11*A12 )/A13;
    static double A21 = pow( -2.*A2q/Acl[1], 4./3. );
    double eta, eta2, Q, dQdeta, d2Qdeta2, g, h, dgdeta, dhdeta, d2gdeta2, d2hdeta2;
    double s, b_ei1, b_ei2, b_ei3, b_ei4, f_inf, Ax;
    double A1, A2, A3;
    double th2, th3, th4, th5, th6, th7, th8, th9, th10, th11;
    double T_MeV, T_p, a_i, Gamma, Gamma_q, Zeq, th, Ath, Bth, Athp, Athpp, Bthp, Bthpp;
    double k_Fn, E_Fn, T_cl, pref, x_r, c_vi_h, c_vi_ah, c_vei, c_vn;
    //parametrization of 1S0 neutron pairing gap. Delta0 in MeV, k0, k2 in fm^{-1}, k1, k3 in fm^{-2}
    double Delta0 = tparams.Delta0_nCrust, k0 = tparams.k0_nCrust, k1 = tparams.k1_nCrust, k2 = tparams.k2_nCrust, k3 = tparams.k3_nCrust;

    double kF0, kF2, Delta_n, T_cn, tau, v_A;

    double R_A = 1.; //superfluidity control function for neutron specific heat capacity. Taken to be 1 if no superfluidity.

    for(size_t i=0; i<c_v.shape()[0]-N_GC+1; i++){

        a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm
        s = 1./( 1. + 0.01*pow( log(Z[i]),1.5 ) + 0.097/(Z[i]*Z[i]) );
        b_ei1 = 1. - 1.1866*pow( Z[i],-0.267 ) + 0.27/Z[i];
        b_ei2 = 1. + 2.25/pow( Z[i],1./3. )*( 1. + 0.684*pow( Z[i],5. ) + 0.222*pow( Z[i],6. ) )/( 1. + 0.222*pow( Z[i],6. ) );
        b_ei3 = 41.5/( 1. + log(Z[i]) );
        b_ei4 = 0.395*log(Z[i]) + 0.347/pow( Z[i],1.5 );

        for(size_t j=0; j<c_v.shape()[1]; j++){

            if(i < N_GC){
                T_MeV = T[N_GC+1-i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            else{
                T_MeV = T[i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            x_r = pow(3.*pi*pi*n_e[i][j],1./3.)*hbarc/M_e;

            Zeq = n_e[i][j]/n_i[i]; //proton number in unit cell. Only differs from Z near crust-core transition
            T_p = sqrt(T_pconst*Z[i]*Zeq*n_i[i]/A[i]); //plasma temperature in MeV
            th = T_p/T_MeV;
            Gamma = alpha_e*hbarc*Z[i]*Zeq/(T_MeV*a_i); //plasma coupling parameter (dimensionless)
            Gamma_q = Gamma/th; //plasma coupling parameter at T=T_p (dimensionless)

            th2 = th*th;
            th3 = th2*th;
            th4 = th3*th;
            th5 = th4*th;
            th6 = th5*th;
            th7 = th6*th;
            th8 = th7*th;
            th9 = th8*th;
            th10 = th9*th;
            th11 = th10*th;
            Ath = an[0] + an[1]*th + an[2]*th2 + an[3]*th3 + an[4]*th4 + an[6]*th6+ an[8]*th8;
            Bth = bn[0] + bn[2]*th2 + bn[4]*th4 + bn[5]*th5 + bn[6]*th6 + bn[7]*th7 + alpha[3]*an[6]*th9 + alpha[4]*an[8]*th11;
            Athp = an[1] + 2.*an[2]*th + 3.*an[3]*th2 + 4.*an[4]*th3 + 6.*an[6]*th5 + 8.*an[8]*th7;
            Athpp = 2.*an[2] + 6.*an[3]*th + 12.*an[4]*th2 + 30.*an[6]*th4 + 56.*an[8]*th6;
            Bthp = 2.*bn[2]*th + 4.*bn[4]*th3 + 5.*bn[5]*th4 + 6.*bn[6]*th5 + 7.*bn[7]*th6 + 9.*alpha[3]*an[6]*th8 + 11.*alpha[4]*an[8]*th10;
            Bthpp = 2.*bn[2] + 12.*bn[4]*th2 + 20.*bn[5]*th3+ 30.*bn[6]*th4 + 42.*bn[7]*th5 + 72.*alpha[3]*an[6]*th7 + 110.*alpha[4]*an[8]*th9;
            c_vi_h = n_i[i]*1e39*k_Bcgs*th2*( pow(alpha[0]/(exp(0.5*alpha[0]*th)-exp(-0.5*alpha[0]*th)),2.) + pow(alpha[1]/(exp(0.5*alpha[1]*th)-exp(-0.5*alpha[1]*th)),2.)
                                     + pow(alpha[2]/(exp(0.5*alpha[2]*th)-exp(-0.5*alpha[2]*th)),2.) + (Athpp*Bth*Bth-2.*Athp*Bthp*Bth+2.*Ath*Bthp*Bthp-Ath*Bth*Bthpp)/(Bth*Bth*Bth) ); //in erg/K/cm^3

            A1 = -2.*A11/th*( 6.*A12*A12*th4 + 3.*A12*th2 + 1. )/pow( 1. + A12*th2,3. ) - 2.*A13/th*( 6.*A14*A14*th4 + 3.*A14*th2 + 1. )/pow( 1. + A14*th2,3. );
            A2 = 1.5*Acl[1]/th2*( 2. + 3.*A21*th4 )/pow( 1. + A21*th4,1.25 );
            A3 = 4.*Acl[2]/th3;
            c_vi_ah = n_i[i]*1e39*k_Bcgs*( A1/Gamma_q + A2/(Gamma_q*Gamma_q) + A3/(Gamma_q*Gamma_q*Gamma_q) );

            eta = q*th;
            if(eta < 6.){
                eta2 = eta*eta;
                g = log( 1. + exp(eta2) );
                h = log( ec - (ec-2.)*exp(-eta2) );
                dgdeta = 2.*eta*exp(eta2)/( 1. + exp(eta2) );
                dhdeta = 2.*(ec-2.)*eta/( 2. - ec + exp(eta2+1.) );
                d2gdeta2 = 2.*exp(eta2)*( 2.*eta2 + exp(eta2) + 1. )/pow( 1. + exp(eta2),2. );
                d2hdeta2 = -2.*( ec - 2. )*( exp(eta2+1.)*(2.*eta2 - 1.) + ec - 2. )/pow( 2. - ec + exp(eta2+1.),2. );
                Q = sqrt( g/h );
                dQdeta = 0.5/Q*( dgdeta/h - g*dhdeta/(h*h) );
                d2Qdeta2 = -0.5*dQdeta/(Q*Q)*( dgdeta/h - g*dhdeta/(h*h) ) + 0.5/Q*( 2.*dhdeta/(h*h)*( g*dhdeta/h - dgdeta ) + d2gdeta2/h - g/(h*h)*d2hdeta2 );
                f_inf = 0.00352*pow(Z[i],2./3.)*b_ei1*sqrt( 1. + b_ei2/(x_r*x_r) );
                Ax = ( b_ei3 + 17.9*x_r*x_r )/( 1. + b_ei4*x_r*x_r );
                c_vei = n_i[i]*1e39*k_Bcgs*f_inf*s*Ax*pow( Q/Gamma, s-1. )*( eta*eta*d2Qdeta2 + (s-1.)*( eta*eta/Q*dQdeta*dQdeta - 2.*eta*dQdeta + Q ) );
            }
            else c_vei = 0.;

            k_Fn = pow(3.*pi*pi*n_nf[i],1./3.); //free neutron Fermi momentum in fm^{-1}
            E_Fn = mu_nf[i] - M_n;//sqrt( k_Fn*k_Fn*hbarc*hbarc + M_n*M_n ) - M_n; //mu_nf-M_n in MeV
            T_cl = 3.*E_Fn/(pi*pi); //in MeV
            pref = ( M_n + 10.*T_MeV )/( M_n + 5.*T_MeV );
            if( E_Fn < 1e-10 ) c_vn = 0.; //if no neutrons (E_Fn ~ 0) set c_vn = 0
            else c_vn = 1.5*pref*n_nf[i]*1e39*k_Bcgs*( 1. - exp( -T_MeV/(pref*T_cl) ) ); //in erg/K/cm^3
            if(tparams.SF == true){
                kF0 = k_Fn - k0;
                kF2 = k_Fn - k2;
                if(kF0 > 0. && kF2 < 0.){
                    Delta_n = Delta0*kF0*kF0/(kF0*kF0+k1)*kF2*kF2/(kF2*kF2+k3);
                    T_cn = 0.5669*Delta_n; //1S0 neutron superfluidity critical temperature in MeV
                    tau = T_MeV/T_cn; //ratio of T/neutron critical temperature
                    if(tau < 0.01) {
                        R_A = 0.;
                    }
                    else if(tau < 1.){
                        v_A = sqrt( 1. - tau )*( 1.456 - 0.157/sqrt(tau) + 1.764/tau );
                        R_A = pow( 0.4186 + sqrt( 1.014049 + 0.251001*v_A*v_A ),2.5 )*exp( 1.456 - sqrt( 2.119936 + v_A*v_A ) );
                    }
                    else R_A = 1.;
                }
                else{
                    R_A = 1.;
                }
            }

            c_v[i][j] += (c_vi_h + c_vi_ah + c_vei + R_A*c_vn)/s_0; //add ion and free neutron contributions to electron contribution to specific heat capacity.
        }
    }

    return;
}

/*
    Computes thermodynamic derivatives Omega_xy for simulations where T is held fixed.
    Fundamental variables of Omega are B, mu_e, T
    Inputs: Bmag: magnitude of magnetic field in reduced units
            mu: electron chemical potential in MeV. Taken to vary only in radial direction
            T: temperature in reduced units
            N_GC: number of ghost cells to exchange
            n_t: lapse function
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Outputs: M: magnetization M = -dOmega_e/dB = dP_e/dB in reduced units (i.e., units of B_0)
             chi: differential magnetic susceptibility chi = dM/dB = -d^2Omega_e/dB^2 = d^2P_e/dB^2 (dimensionless)
             M_mu: mixed susceptibility M_mu = -d^2Omega_e/dBdmu_e = d^2P_e/dBdmu_e = dM/dmu_e in MeV^-1 (divide by B_0 also or would have dimensions G/MeV)
*/
void Omega_xyFixedT(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & M, ScalarField & chi, ScalarField & M_mu)
{
    static double sqrtalpha_e = sqrt(alpha_e);
    static double M2 = M_e*M_e;
    static double eB_pref = unit_e*B_0*statCGtoMeV2; //constant prefactor for computing e*B in MeV^2
    static double pi2 = pi*pi;
    static double Tconv = k_B*T_0; //T_0 in MeV. Used to convert T in units of T_0 to MeV
    double eB, T_MeV, n_max, thermalBFactor, p_F; //fundamental charge times B in MeV^2, maximum occupied Landau level at T=0, dimensionless thermal damping factor for de Haas-van Alphen oscillations, Fermi momentum in MeV

    //T_B = 100*M*( sqrt(1+2*(n_max+1)*eB/eB_crit) - sqrt(1+2*n_max*eB/eB_crit) ); //critical temperature in MeV

    double dPdBTemp, d2PdB2Temp, d2PdBdmuTemp;

    for(size_t i=0; i<Bmag.shape()[0]-N_GC; i++){
        for(size_t j=0; j<Bmag.shape()[1]; j++){
            if( Bmag[i][j] < 1e-5 ){
                M[i][j] = 0.;
                chi[i][j] = 0.;
                M_mu[i][j] = 0.;
            }
            else {
                eB = Bmag[i][j]*eB_pref; //e*B in MeV^2
                p_F = sqrt( mu[i]*mu[i]-M2 ); //Fermi momentum in MeV
                // long double ratio = ((long double)mu[i]*(long double)mu[i]-(long double)M2)/(2.*(long double)eB);
                // n_max = (double)std::floor( ratio );
                // n_max = stable_floor( (mu[i]*mu[i]-M2)/(2.*eB) );
                n_max = std::floor( (mu[i]*mu[i]-M2)/(2.*eB) );
                T_MeV = Tconv*T[i][j]/n_t[i];
                thermalBFactor = 2.*pi2*mu[i]*T_MeV/eB;
                // if ( j <= 21 && static_cast<int>(std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][j]*eB_pref) )) != static_cast<int>(std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][Ntheta-j-1]*eB_pref) )) )
                // if ( j <= 21 && static_cast<int>(std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][j]*eB_pref) )) != static_cast<int>(std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][Ntheta-j-1]*eB_pref) )) )
                // {
                //     std:: cout << "i = " << i << ", j = " << j << ", Deltan_l = " << std::setprecision(15) << (mu[i]*mu[i]-M2)/(2.*eB)-n_max << ", Deltan_r = "
                //             << std::setprecision(15) << (mu[i]*mu[i]-M2)/(2.*Bmag[i][Ntheta-j-1]*eB_pref)-std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][Ntheta-j-1]*eB_pref) ) << std::endl;
                //     // std:: cout << "n_max_l = " << std::setprecision(15) << n_max << ", n_max_r = " << std::setprecision(15) << std::floor( (mu[i]*mu[i]-M2)/(2.*Bmag[i][Ntheta-j-1]*eB_pref) ) << std::endl;
                // }
                // std:: cout << "i = " << i << ", j = " << j << ", n_max = " << n_max << std::endl;

                dPdBTemp = 0.; //zero dPdBTemp, d2PdB2Temp and d2PdBdmuTemp
                d2PdB2Temp = 0.;
                d2PdBdmuTemp = 0.;

                if(thermalBFactor > 1.){
                    Omega_xyHighT_FixedT(eB, mu[i], T_MeV, p_F, OmB_Interpolators, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp);
                }
                else{
                    Omega_xyLowT_FixedT_EMSums(eB, mu[i], T_MeV, p_F, n_max, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp);
                    Omega_xyLowT_FixedT_MaxLL(eB, mu[i], T_MeV, n_max, OmB_Interpolators, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp);
                }

                M[i][j] = sqrtalpha_e*dPdBTemp*MeV2toGauss/B_0; //in MeV^2 -> G -> reduced units
                if( 1. - 4.*pi*alpha_e*d2PdB2Temp < 0.) chi[i][j] = 1./(4.*pi); //removes unstable regions (requires 1-4*pi*d2PdB2 = 1 + 4*pi*chi > 0)
                else chi[i][j] = alpha_e*d2PdB2Temp;
                M_mu[i][j] = sqrtalpha_e*( d2PdBdmuTemp )*MeV2toGauss/B_0; //in MeV -> G/MeV -> reduced units/MeV
            }
        }
    }

    return;

}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is held fixed.
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Outputs: dPdBTemp: dP_e/dB in MeV^2
             d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
             d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
*/
void Omega_xyLowT_FixedT_EMSums(double eB, double mu, double T, double p_F, double n_max, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp)
{
    //double T_B = M_e*( sqrt(1. + 2.*(n_max+1.)*eB/eB_crit) - sqrt(1. + 2.*n_max*eB/eB_crit) ); //critical temperature in MeV

    static double B_2 = 1./6., B_4 = -1./30.;
    static double pi2 = pi*pi;
    static double M2 = M_e*M_e;
    static double M4 = M2*M2;
    double npr = n_max - 1.;
    double p_F2 = p_F*p_F;
    double E_Fnpr2 = p_F2 - 2.*eB*npr;
    double E_Fnpr = sqrt( E_Fnpr2 );
    double InvE_Fnpr = 1./E_Fnpr;
    double MagMnpr_2 = M2 + 2.*eB*npr;
    double MagMnpr = sqrt( MagMnpr_2 );
    double InvMagMnpr = 1./MagMnpr;
    double E_F1 = sqrt( p_F2 - 2.*eB );
    double InvE_F1 = 1./E_F1;
    double MagM1 = sqrt( M2 + 2.*eB );
    double InvMagM1 = 1./MagM1;
    double mu2 = mu*mu;
    double mu4 = mu2*mu2;
    double eB2 = eB*eB;

    double MagMnpr_4, MagM1_2, MagM1_4, InvMagMnpr_2, InvMagMnpr_4, InvMagM1_2, InvMagM1_4;
    double InvE_Fnpr_3, InvE_Fnpr_5, InvE_Fnpr_7, InvE_F1_3, InvE_F1_5, InvE_F1_7, logFact_npr, logFact_1;

    if( n_max > 1.1 ){
        MagMnpr_4 = MagMnpr_2;
        MagM1_2 = M2 + 2.*eB;
        MagM1_4 = MagM1_2*MagM1_2;
        InvMagMnpr_2 = 1./MagMnpr_2;
        InvMagMnpr_4 = InvMagMnpr_2*InvMagMnpr_2;
        InvMagM1_2 = 1./MagM1_2;
        InvMagM1_4 = InvMagM1_2*InvMagM1_2;
        InvE_Fnpr_3 = InvE_Fnpr*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_5 = InvE_Fnpr_3*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_7 = InvE_Fnpr_5*InvE_Fnpr*InvE_Fnpr;
        InvE_F1_3 = InvE_F1*InvE_F1*InvE_F1;
        InvE_F1_5 = InvE_F1_3*InvE_F1*InvE_F1;
        InvE_F1_7 = InvE_F1_5*InvE_F1*InvE_F1;
        logFact_npr = log( (mu+E_Fnpr)*InvMagMnpr );
        logFact_1 = log( (mu+E_F1)*InvMagM1 );

        dPdBTemp += 1./(4.*pi2)*( mu*p_F-M2*log((mu+p_F)/M_e) )
                    + 0*1./(2.*pi2)*( npr*( mu*E_Fnpr - MagMnpr_2*logFact_npr ) - ( mu*E_F1 - MagM1_2*logFact_1 ) )
                    + 1./(4.*pi2)*( mu*E_Fnpr - (M2+4.*eB*npr)*logFact_npr + mu*E_F1 - (M2+4.*eB)*logFact_1 )
                    + B_2*eB/(2.*pi2)*( eB*npr*mu*InvMagMnpr_2*InvE_Fnpr - 2.*logFact_npr
                                          - eB*mu*InvMagM1_2*InvE_F1 + 2.*logFact_1 )
                    + B_4*eB2*eB*mu/(24.*pi2)*( ( 20.*mu2*(M2+eB*npr)*MagMnpr_2
                                                    - 3.*(4.*M2 + 3.*eB*npr)*MagMnpr_4
                                                    - 8.*mu4*(M2+eB*npr) )*InvMagMnpr_4*InvMagMnpr_2*InvE_Fnpr_5
                                                  - ( 20.*mu2*(M2+eB)*MagM1_2
                                                    - 3.*(4.*M2 + 3.*eB)*MagM1_4
                                                    - 8.*mu4*(M2+eB) )*InvMagM1_4*InvMagM1_2*InvE_F1_5 )
                        + mu*T*T/6.*( 0.5/p_F + ( npr*InvE_Fnpr - InvE_F1 ) + 0.5*( (p_F2 - eB*npr)*InvE_Fnpr_3 + (p_F2 - eB)*InvE_F1_3 )
                                             + B_2*eB/2.*( (2.*p_F2 - eB*npr)*InvE_Fnpr_5 - (2.*p_F2 - eB)*InvE_F1_5 ) );

        d2PdB2Temp += 1./pi2*( -npr*npr*logFact_npr + logFact_1)
                     + 1./(2.*pi2)*( npr*( eB*npr*mu*InvMagMnpr_2*InvE_Fnpr - 2.*logFact_npr )
                                                  + ( eB*mu*InvMagM1_2*InvE_F1 - 2.*logFact_1 ) )
                             + B_2/(2.*pi2)*( -eB*npr*mu*( 4.*M4 + M2*(13.*eB*npr - 4.*mu2) + 2.*eB*npr*(5.*eB*npr - 3.*mu2) )*InvMagMnpr_4*InvE_Fnpr_3
                                                      - 2.*logFact_npr
                                                + eB*mu*( 4.*M4 + M2*(13.*eB - 4.*mu2) + 2.*eB*(5.*eB - 3.*mu2) )*InvMagM1_4*InvE_F1_3
                                                      + 2.*logFact_1 )
                    + B_4*eB2*mu/(24.*pi2)*( -( 2.*mu2*MagMnpr_4*( 48.*M4 + 52.*M2*eB*npr + 17.*eB2*npr*npr )
                                                                 - 3.*MagMnpr_4*MagMnpr_2*( 12.*M4 + 8.*M2*eB*npr + 3.*eB2*npr*npr )
                                                                 + 8.*mu4*mu2*(3.*M4 + 4.*M2*eB*npr + 2*eB2*npr*npr)
                                                                 - 28.*mu4*MagMnpr_2*(3.*M4 + 4.*M2*eB*npr + 2.*eB2*npr*npr) )*InvMagMnpr_4*InvMagMnpr_4*InvE_Fnpr_7
                                                          + ( 2.*mu2*MagM1_4*( 48.*M4 + 52.*M2*eB + 17.*eB2 )
                                                                 - 3.*MagM1_4*MagM1_2*(12.*M4 + 8.*M2*eB + 3.*eB2)
                                                                 + 8.*mu4*mu2*(3.*M4 + 4.*M2*eB + 2.*eB2)
                                                                 - 28.*mu4*MagM1_2*(3.*M4 + 4.*M2*eB + 2.*eB2) )*InvMagM1_4*InvMagM1_4*InvE_F1_7 )
                                    + T*T*mu/6.*( (npr*npr*InvE_Fnpr_3 - InvE_F1_3) + 0.5*( (2.*p_F2-eB*npr)*InvE_Fnpr_5*npr + (2.*p_F2-eB)*InvE_F1_5 )
                                                                      + B_2/2.*( ( 2.*M4 - eB2*npr*npr + 4.*eB*npr*mu2 + 2.*mu4 - 4.*M2*(eB*npr + mu2) )*InvE_Fnpr_7
                                                                               - ( 2.*M4 - eB2 + 4.*eB*mu2 + 2.*mu4 - 4.*M2*(eB + mu2) )*InvE_F1_7 ) );;

        d2PdBdmuTemp += 1./(2.*pi2)*( ( p_F2-3.*eB*npr )*InvE_Fnpr + ( p_F2-3.*eB )*InvE_F1 )
                        + p_F/(2.*pi2) + 1./pi2*( npr*E_Fnpr - E_F1 )
                        + B_2*eB/(2.*pi2)*( -( 2.*p_F2 - 3.*eB*npr )*InvE_Fnpr_3 + ( 2.*p_F2 - 3.*eB )*InvE_F1_3 )
                        + B_4*eB2*eB/(8.*pi2)*( -( 4.*p_F2 - 3.*eB*npr )*InvE_Fnpr_7 + ( 4.*p_F2 - 3.*eB )*InvE_F1_7 )
                        - T*T/6.*( 0.5*M2*p_F2/pow(p_F2-M2,2.5) + npr*MagMnpr*InvE_Fnpr_3 - MagM1*InvE_F1_3
                                  + 0.5*( (mu2*(M2 + 4.*eB*npr) - MagMnpr_2*(M2 + eB*npr))*InvE_Fnpr_5 + (mu2*(M2 + 4.*eB) - MagM1_2*(M2 + eB))*InvE_F1_5 )
                                  + B_2*eB/2.*( ( 4.*mu4 - 2.*mu2*(M2 - 4.*eB*npr) - 2.*eB2*npr*npr - 5.*eB*npr*M2 - 2.*M4 )*InvE_Fnpr_7
                                              - ( 4.*mu4 - 2.*mu2*(M2 - 4.*eB) - 2.*eB2 - 5.*eB*M2 - 2.*M4 )*InvE_F1_7 ) );
    }
    else if( n_max > 0.5 ){
        //dPdB, d2PdB2, d2PdBdmu contribution from n=0
        dPdBTemp += 1./(2.*pi2)*( 0.5*mu*p_F - 0.5*M2*log( (mu+p_F)/M_e ) + pi2/6.*mu/p_F*T*T + 7.*pi2*pi2/120.*mu*M2/(p_F2*p_F2*p_F)*T*T*T*T );
        d2PdB2Temp += 0.;
        d2PdBdmuTemp += 1./(2.*pi2)*( p_F - pi2/6.*M2/(p_F2*p_F)*T*T - 7.*pi2*pi2/120.*M2*(4.*mu2 + M2)/(p_F2*p_F2*p_F2*p_F)*T*T*T*T );
    }
    else{
        //All contributions to dPdB, d2PdB2, d2PdBdmu when n_max = 0 come from "Omega_xyLowT_FixedT_MaxLL"
        dPdBTemp += 0.;
        d2PdB2Temp += 0.;
        d2PdBdmuTemp += 0.;
    }

    return;
}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is held fixed.
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
    Outputs: dPdBTemp: dP_e/dB in MeV^2
             d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
             d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
*/
void Omega_xyLowT_FixedT_MaxLL(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp)
{
    static double pi2 = pi*pi;
    static double pi4 = pi2*pi2;
    static double M2 = M_e*M_e;
    double M_n, M_np1, a, a1, b, pref;
    double G1_ab = 0, G2_ab = 0, H2_ab = 0, I1_ab = 0, I2_ab = 0;
    double G1_ab1 = 0, G2_ab1 = 0, H2_ab1 = 0, I1_ab1 = 0, I2_ab1 = 0;
    double LiN1_2, Li1_2, Li3_2, Li5_2, Li7_2; //temporary variables to hold interpolated polylogarithm function evaluations
    double KroneckerDelta = 0; //Kronecker delta function for n_max, 0; returns 1 if n_max = 0 and zero otherwise. Assume 0 and then correct otherwise.

    double a2, a3, z, diffz, zerr;

    M_n = sqrt(M2 + 2.*eB*n_max); //M_{n_max}
    M_np1 = sqrt(M2 + 2.*eB*(n_max+1)); //M_{n_max+1}

    // n = n_max contribution
    a = M_n/T;
    b = mu/T;

    if( b - a > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab = b*sqrt( b*b/(a*a)-1. )/(2.*a) + 0.5*log( sqrt( b*b/(a*a)-1. ) + b/a ) + pi2*b*(b*b-2.*a*a)/(6.*a*a*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(4.*a*a + b*b)/pow(b*b-a*a,3.5);
        G2_ab = log( sqrt( b*b/(a*a)-1. ) + b/a ) - pi2*b/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        H2_ab = -b/(a*sqrt( b*b/(a*a)-1. )) - pi2*b*a*a/(2.*pow(b*b-a*a,2.5)) - 7.*pi4/24.*b*a*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,4.5);
        I1_ab = sqrt( b*b/(a*a)-1. ) - pi2*a/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        I2_ab = -1./sqrt( b*b/(a*a)-1. ) - pi2*a*(2.*b*b + a*a)/(6.*pow(b*b-a*a,2.5)) - 7.*pi4/120.*a*(8.*b*b*b*b + 24.*a*a*b*b + 3.*a*a*a*a)/pow(b*b-a*a,4.5);
    }
    else if(a - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a))*exp(b-a);
        G1_ab = pref*( 1. + 7./(8.*a) + 57./(128.*a*a) - 195./(1024.*a*a*a) );
        G2_ab = pref*( 1. - 1./(8.*a) + 9./(128.*a*a) - 75./(1024.*a*a*a) );
        H2_ab = -pref*( a + 3./8. - 15./(128.*a) + 105./(1024.*a*a) );
        I1_ab = pref*( 1. + 3./(8.*a) - 15./(128.*a*a) + 105./(1024.*a*a*a) );
        I2_ab = -pref*( a - 1./8. + 9./(128.*a) - 75./(1024.*a*a) );
    }
    else if( sqrt(a*a+b*b) < 50. ){ //full interpolating functions
        G1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else{ //polylogarithms with argument -exp(b-a), using interpolating functions for speed

        diffz = mu - M_n;
        zerr = (mu-diffz) - M_n;
        z = diffz + zerr;
        // z = ( mu*mu - M_n*M_n )/( mu + M_n );
        // double x = (mu*mu-M2)/(2.*eB);
        // z = 2.*eB/(mu+M_n)*( 0.5 - 0.5/asin(tanh(5.))*asin( tanh(5.*sin(pi*x))*cos(pi*x) ) );
        // z = ( mu*mu*mu*mu - M_n*M_n*M_n*M_n )/( mu + M_n )/( mu*mu + M_n*M_n );
        // z = 0.01;//
        z = z/T;
        // if (b-a < 20.) {
            // LiN1_2term = sqrt(pi*a/2.)*gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, z, OmB_Interpolators.PolyLogN1_2_acc );
        // }
        // else {
        //     LiN1_2term = sqrt(pi/2.)*( - 1./sqrt(b/a-1.)/sqrt(pi) - sqrt(a)*pi*sqrt(pi)/8./pow(b-a,2.5) - sqrt(a)*49.*pi*pi*pi*sqrt(pi)/384./pow(b-a,4.5) );
        // }
        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, z, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, z, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, z, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, z, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, z, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a));
        a2 = a*a;
        a3 = a*a2;
        G1_ab = -pref*( Li1_2 + 7./(8.*a)*Li3_2 + 57./(128.*a2)*Li5_2 - 195./(1024.*a3)*Li7_2 );
        G2_ab = -pref*( Li1_2 - 1./(8.*a)*Li3_2 + 9./(128.*a2)*Li5_2 - 75./(1024.*a3)*Li7_2 );
        // H2_ab = LiN1_2term + pref*( 3./8.*Li1_2 - 15./(128.*a)*Li3_2 + 105./(1024.*a2)*Li5_2 );
        H2_ab = pref*( a*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a)*Li3_2 + 105./(1024.*a2)*Li5_2 );
        I1_ab = -pref*( Li1_2 + 3./(8.*a)*Li3_2 - 15./(128.*a2)*Li5_2 + 105./(1024.*a3)*Li7_2 );
        // I2_ab = LiN1_2term + pref*( - 1./8.*Li1_2 + 9./(128.*a)*Li3_2 - 75./(1024.*a2)*Li5_2 );
        I2_ab = pref*( a*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a)*Li3_2 - 75./(1024.*a2)*Li5_2 );
    }

    // n = n_max + 1 contribution

    a1 = M_np1/T;

    if( b - a1 > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab1 = b*sqrt( b*b/(a1*a1)-1. )/(2.*a1) + 0.5*log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) + pi2*b*(b*b-2.*a1*a1)/(6.*a1*a1*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(4.*a1*a1 + b*b)/pow(b*b-a1*a1,3.5);
        G2_ab1 = log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) - pi2*b/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        H2_ab1 = -b/(a1*sqrt( b*b/(a1*a1)-1. )) - pi2*b*a1*a1/(2.*pow(b*b-a1*a1,2.5)) - 7.*pi4/24.*b*a1*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,4.5);
        I1_ab1 = sqrt( b*b/(a1*a1)-1. ) - pi2*a1/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        I2_ab1 = -1./sqrt( b*b/(a1*a1)-1. ) - pi2*a1*(2.*b*b + a1*a1)/(6.*pow(b*b-a1*a1,2.5)) - 7.*pi4/120.*a1*(8.*b*b*b*b + 24.*a1*a1*b*b + 3.*a1*a1*a1*a1)/pow(b*b-a1*a1,4.5);
    }
    else if(a1 - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a1))*exp(b-a1);
        G1_ab1 = pref*( 1. + 7./(8.*a1) + 57./(128.*a1*a1) - 195./(1024.*a1*a1*a1) );
        G2_ab1 = pref*( 1. - 1./(8.*a1) + 9./(128.*a1*a1) - 75./(1024.*a1*a1*a1) );
        H2_ab1 = -pref*( a1 + 3./8. - 15./(128.*a1) + 105./(1024.*a1*a1) );
        I1_ab1 = pref*( 1. + 3./(8.*a1) - 15./(128.*a1*a1) + 105./(1024.*a1*a1*a1) );
        I2_ab1 = -pref*( a1 - 1./8. + 9./(128.*a1) - 75./(1024.*a1*a1) );
    }
    else if( sqrt(a1*a1+b*b) < 50. ){ //full interpolating functions
        G1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a1, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a1, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a1, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a1, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a1, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else {
        //polylogarithms with argument -exp(b-a), using interpolating functions for speed

        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, b-a1, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, b-a1, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, b-a1, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, b-a1, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, b-a1, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a1));
        G1_ab1 = -pref*( Li1_2 + 7./(8.*a1)*Li3_2 + 57./(128.*a1*a1)*Li5_2 - 195./(1024.*a1*a1*a1)*Li7_2 );
        G2_ab1 = -pref*( Li1_2 - 1./(8.*a1)*Li3_2 + 9./(128.*a1*a1)*Li5_2 - 75./(1024.*a1*a1*a1)*Li7_2 );
        H2_ab1 = pref*( a1*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a1)*Li3_2 + 105./(1024.*a1*a1)*Li5_2 );
        I1_ab1 = -pref*( Li1_2 + 3./(8.*a1)*Li3_2 - 15./(128.*a1*a1)*Li5_2 + 105./(1024.*a1*a1*a1)*Li7_2 );
        I2_ab1 = pref*( a1*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a1)*Li3_2 - 75./(1024.*a1*a1)*Li5_2 );
    }

    if(n_max < 0.5) KroneckerDelta = 1.;

    dPdBTemp += (2. - KroneckerDelta)/(2.*pi2)*( M_n*M_n*G1_ab - (M_n*M_n + eB*n_max)*G2_ab ) + 1./pi2*( M_np1*M_np1*G1_ab1 - (M_np1*M_np1 + eB*(n_max+1.))*G2_ab1 );
    d2PdB2Temp += -n_max/(2.*pi2)*(2. - KroneckerDelta)*( 2.*G2_ab + eB*n_max/(M_n*M_n)*H2_ab ) - (n_max+1.)/pi2*( 2.*G2_ab1 + eB*(n_max+1.)/(M_np1*M_np1)*H2_ab1 );
    d2PdBdmuTemp += (2. - KroneckerDelta)/(2.*pi2)*M_n*( I1_ab + eB*n_max/(M_n*M_n)*I2_ab ) + 1./pi2*M_np1*( I1_ab1 + eB*(n_max+1.)/(M_np1*M_np1)*I2_ab1 );

    return;

}

/*
    High-temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is held fixed.
    Consists of regular (non-oscillatory) and oscillatory terms with finite T thermal damping. Derived using the
    Poisson summation formula, and based on the Lifshitz-Kosevich formula.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Outputs: dPdBTemp: dP_e/dB in MeV^2
             d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
             d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
*/
void Omega_xyHighT_FixedT(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp)
{
    //T_B = 100*M*( sqrt(1+2*(n_max+1)*eB/eB_crit) - sqrt(1+2*n_max*eB/eB_crit) ); //critical temperature in MeV
    static double zeta3_2 = 2.61237534869; //zeta(3/2)
    static double sqrtpi = sqrt(pi);
    static double pi2 = pi*pi; //pi^2
    static double pi2_5 = pi*pi*sqrt(pi); //pi^2.5
    static double pi3 = pi2*pi; //pi^3
    static double sqrt2opi = sqrt(2./pi); //sqrt(2/pi)
    static double M2 = M_e*M_e;
    static double InvSq2 = 1./sqrt(2.); // 1/sqrt(2)
    static int p_max = 5; //maximum number of terms to include in high-T approximation sum

    double n, sqrtn, lambda, R_T, Sin, Cos, phase; //phase = phase of sines and cosines divided by pi
    double i1val, i2val, h1val, h2val;
    double SinP, CosP;

    double sqrteB = sqrt(eB);
    double eB_snap;
    double p_F2 = p_F*p_F;

    // First compute regular (non-oscillatory part) of the terms. Add them to zeroed variables.
    double x = p_F2/eB;
    if( x < 1e4 ){
        i1val = gsl_spline_eval( OmB_Interpolators.i_1, x, OmB_Interpolators.i_1_acc );
        i2val = gsl_spline_eval( OmB_Interpolators.i_2, x, OmB_Interpolators.i_2_acc );
    }
    else{
        i1val = sqrtpi/sqrt(x)/3.;
        i2val = sqrtpi/(6.*x*sqrt(x));
    }

    h1val = gsl_spline2d_eval( OmB_Interpolators.h_1, mu/M_e, eB/M2, OmB_Interpolators.h_1_acc_a, OmB_Interpolators.h_1_acc_b );
    h2val = gsl_spline2d_eval( OmB_Interpolators.h_2, mu/M_e, eB/M2, OmB_Interpolators.h_2_acc_a, OmB_Interpolators.h_2_acc_b );

    dPdBTemp += sqrteB/(2.*pi2_5)*h1val - M2/(4.*pi2_5*sqrteB)*h2val - 1./(8.*pi2_5)*sqrteB*( mu*i1val - sqrt2opi*M_e*zeta3_2 );
    d2PdB2Temp += 1./(2.*pi2_5*sqrteB)*h1val - M2/(2.*pi2_5*eB*sqrteB)*h2val - 5./(16.*pi2_5*sqrteB)*( mu*i1val - sqrt2opi*M_e*zeta3_2 );
    d2PdBdmuTemp += 3.*sqrteB/(8.*pi2_5)*i1val + p_F2/(4.*pi2_5*sqrteB)*i2val;

    // Next compute oscillatory part of the terms, involving a sum over p=0,...,p_max
    for(int p = 1; p <= p_max; p++) {
        n = static_cast<double>(p);
        sqrtn = sqrt(n);
        lambda = 2.*pi2*mu*T*n/eB; //lambda factor as defined in notes times pi
        if(lambda < 20.) R_T = lambda/sinh(lambda);
        else R_T = 2.*lambda*exp(-lambda); //use large lambda expansion for sinh for lambda > 20

        //snap eB to limit its precision- prevents amplification of floating point rounding error in trigonometric function evaluation
        eB_snap = std::ldexp(std::round(std::ldexp(eB, 27)), -27);
        phase = std::remainder( n*p_F2/eB_snap,2. ); //n*p_F2/eB_snap mod 2 (only relevant part of phase/pi)
        // std::cout << std::setprecision(15) << eB_snap << ", " << eB << std::endl;
        CosP = cos( phase*pi );
        SinP = sin( phase*pi );

        Sin = InvSq2*( CosP - SinP ); //sin( pi/4 - pi*n*p_F2/eB ) = 1/sqrt(2)*( cos(pi*n*p_F2/eB) - sin(pi*n*p_F2/eB) )
        Cos = InvSq2*( CosP + SinP ); //cos( pi/4 - pi*n*p_F2/eB ) = 1/sqrt(

        dPdBTemp += R_T*( sqrteB*p_F2/(4.*pi3*mu*n*sqrtn)*Sin - 5.*eB*sqrteB/(8.*pi3*pi*mu*n*n*sqrtn)*Cos );
        d2PdB2Temp += R_T*( p_F2*p_F2/(4.*pi2*mu*eB*sqrteB*sqrtn)*Cos + 3.*p_F2/(4.*pi3*mu*sqrteB*n*sqrtn)*Sin );
        d2PdBdmuTemp += -R_T*( 3.*sqrteB/(4*pi3*n*sqrtn)*Sin + p_F2/(2.*pi2*sqrteB*sqrtn)*Cos );

    }

    return;
}

/*
    Computes thermodynamic derivatives Omega_xy for simulations where T is varying.
    Fundamental variables of Omega are B, mu_e, T
    Inputs: Bmag: magnitude of magnetic field in reduced units
            mu: electron chemical potential in MeV. Taken to vary only in radial direction
            T: temperature in reduced units
            N_GC: number of ghost cells to exchange
            n_t: lapse function
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Outputs: M: magnetization M = -dOmega_e/dB = dP_e/dB in reduced units (i.e., units of B_0)
             chi: differential magnetic susceptibility chi = dM/dB = -d^2Omega_e/dB^2 = d^2P_e/dB^2 (dimensionless)
             M_mu: mixed susceptibility M_mu = -d^2Omega_e/dBdmu_e = d^2P_e/dBdmu_e = dM/dmu_e in MeV^-1 (divide by B_0 also or would have dimensions G/MeV)
             c_v_eB: specific heat capacity of magnetized electrons c_v = T*ds/dT = -T*d^2Omega_e/dT^2 = T*d^2P_e/dT^2 in reduced units.
             C_m: magnetocaloric coefficient C_m = ds/dB = -d^2Omega_e/dTdB = d^2P_e/dBdT in reduced units (i.e. times B_0/s_0)
*/
void Omega_xyVaryingT(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & M, ScalarField & chi, ScalarField & M_mu, ScalarField & c_v_eB, ScalarField & C_m)
{
    static double sqrtalpha_e = sqrt(alpha_e);
    static double M2 = M_e*M_e;
    static double eB_pref = unit_e*B_0*statCGtoMeV2; //constant prefactor for computing e*B in MeV^2
    double eB, n_max, thermalBFactor, p_F; //fundamental charge times B in MeV^2, maximum occupied Landau level at T=0, dimensionless thermal damping factor for de Haas-van Alphen oscillations, Fermi momentum in MeV
//    double T_B;
    double T_MeV; //local temperature in MeV

    double dPdBTemp, d2PdB2Temp, d2PdBdmuTemp, d2PdT2Temp, d2PdTdBTemp;


    for(size_t i=0; i<Bmag.shape()[0]-N_GC; i++){
        p_F = sqrt( mu[i]*mu[i]-M2 ); //Fermi momentum in MeV

        for(size_t j=0; j<Bmag.shape()[1]; j++){

            if(i < N_GC){
                T_MeV = T[N_GC+1-i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            else{
                T_MeV = T[i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }

//            c_v_eB[i][j] = ( T_MeV*mu[i]*p_F/3. )*k_B*MeVtoErg/(hbarc*hbarc*hbarc)*1e39/s_0;
            if( Bmag[i][j] < 1e-5 ){
                M[i][j] = 0.;
                chi[i][j] = 0.;
                M_mu[i][j] = 0.;
                c_v_eB[i][j] = ( T_MeV*mu[i]*p_F/3. )*k_B*MeVtoErg/(hbarc*hbarc*hbarc)*1e39/s_0;
                C_m[i][j] = 0.;
            }
            else{
                eB = Bmag[i][j]*eB_pref; //e*B in MeV^2
                n_max = std::floor( (mu[i]*mu[i]-M2)/(2.*eB) ); //highest occupied Landau level at T=0
                // double x = (mu[i]*mu[i]-M2)/(2.*eB);
                // int nmax = static_cast<int>(std::floor(x));

                thermalBFactor = 2.*pi*pi*mu[i]*T_MeV/eB;
                //T_B = M_E*( sqrt(1. + 2.*(n_max+1.)*eB/eB_crit) - sqrt(1. + 2.*n_max*eB/eB_crit) );

                dPdBTemp = 0.; //zero dPdBTemp, d2PdB2Temp and d2PdBdmuTemp
                d2PdB2Temp = 0.;
                d2PdBdmuTemp = 0.;
                d2PdT2Temp = 0.;
                d2PdTdBTemp = 0.;

                // if(thermalBFactor > 1.){
                    Omega_xyHighT_VaryingT(eB, mu[i], T_MeV, p_F, OmB_Interpolators, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp, d2PdT2Temp, d2PdTdBTemp);
                // }
                // else{
                //     Omega_xyLowT_VaryingT_EMSums(eB, mu[i], T_MeV, p_F, n_max, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp, d2PdT2Temp, d2PdTdBTemp);
                //     Omega_xyLowT_VaryingT_MaxLL(eB, mu[i], T_MeV, n_max, OmB_Interpolators, dPdBTemp, d2PdB2Temp, d2PdBdmuTemp, d2PdT2Temp, d2PdTdBTemp);
                // }

                M[i][j] = sqrtalpha_e*dPdBTemp*MeV2toGauss/B_0; //in MeV^2 -> G -> reduced units
                if( 1. - 4.*pi*alpha_e*d2PdB2Temp < 0. ) chi[i][j] = 1./(4.*pi); //removes unstable regions (requires 1-4*pi*d2PdB2 = 1 - 4*pi*chi > 0)
                else chi[i][j] = alpha_e*d2PdB2Temp;
                M_mu[i][j] = sqrtalpha_e*( d2PdBdmuTemp )*MeV2toGauss/B_0; //in MeV -> G/MeV -> reduced units/MeV
                c_v_eB[i][j] = ( T_MeV*d2PdT2Temp )*k_B*MeVtoErg/(hbarc*hbarc*hbarc)*1e39/s_0; //in MeV^3 -> MeV^4/K -> erg/(fm^3 K) -> erg/(cm^3 K) -> reduced units
//                if( T[i][j] < T_B ) c_v_eB[i][j] = ( k_B*T[i][j]*T_0*d2PdT2Temp )*k_B*MeVtoErg/pow(hbarc,3.)*1e39/s_0; //in MeV^3 -> MeV^4/K -> erg/(fm^3 K) -> erg/(cm^3 K) -> reduced units
//                else ( k_B*T[i][j]*T_0*mu[i]*p_F/3. )*k_B*MeVtoErg/pow(hbarc,3.)*1e39/s_0;
                C_m[i][j] = k_B*sqrtalpha_e*( d2PdTdBTemp )*MeV2toGauss*T_0/B_0; //in MeV -> MeV^2/K -> G/K -> reduced units
            }
        }
    }

    return;

}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying.
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Output: dPdBTemp: dP_e/dB in MeV^2
            d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
            d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
            d2PdT2Temp: d^2P/dT^2 in MeV^2
            d2PdTdBTemp: d^3P/dTdB in MeV
*/
void Omega_xyLowT_VaryingT_EMSums(double eB, double mu, double T, double p_F, double n_max, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp)
{
    //double T_B = M_e*( sqrt(1. + 2.*(n_max+1.)*eB/eB_crit) - sqrt(1. + 2.*n_max*eB/eB_crit) ); //critical temperature in MeV

    static double pi2 = pi*pi;
    static double pi4 = pi2*pi2;
    static double B_2 = 1./6., B_4 = -1./30.;
    static double M2 = M_e*M_e;
    static double M4 = M2*M2;
    double p_F2 = p_F*p_F;
    double npr = n_max - 1.;
    double E_Fnpr = sqrt( p_F2 - 2.*eB*npr );
    double InvE_Fnpr = 1./E_Fnpr;
    double MagMnpr = sqrt( M2 + 2.*eB*npr );
    double InvMagMnpr = 1./MagMnpr;
    double E_F1 = sqrt( p_F2 - 2.*eB );
    double InvE_F1 = 1./E_F1;
    double MagM1 = sqrt( M2 + 2.*eB );
    double InvMagM1 = 1./MagM1;
    double mu2 = mu*mu;
    double mu4 = mu*mu*mu*mu;
    double eB2 = eB*eB;

    double MagMnpr_2, MagMnpr_4, MagM1_2, MagM1_4, InvMagMnpr_2, InvMagMnpr_4, InvMagM1_2, InvMagM1_4;
    double InvE_Fnpr_3, InvE_Fnpr_5, InvE_Fnpr_7, InvE_Fnpr_9, InvE_F1_3, InvE_F1_5, InvE_F1_7, InvE_F1_9, logFact_npr, logFact_1;

    if( n_max > 1.1 ){
        MagMnpr_2 = M2 + 2.*eB*npr;
        MagMnpr_4 = MagMnpr_2*MagMnpr_2;
        MagM1_2 = M2 + 2.*eB;
        MagM1_4 = MagM1_2*MagM1_2;
        InvMagMnpr_2 = InvMagMnpr*InvMagMnpr;
        InvMagMnpr_4 = InvMagMnpr_2*InvMagMnpr_2;
        InvMagM1_2 = InvMagM1*InvMagM1;
        InvMagM1_4 = InvMagM1_2*InvMagM1_2;
        InvE_Fnpr_3 = InvE_Fnpr*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_5 = InvE_Fnpr_3*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_7 = InvE_Fnpr_5*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_9 = InvE_Fnpr_7*InvE_Fnpr*InvE_Fnpr;
        InvE_F1_3 = InvE_F1*InvE_F1*InvE_F1;
        InvE_F1_5 = InvE_F1_3*InvE_F1*InvE_F1;
        InvE_F1_7 = InvE_F1_5*InvE_F1*InvE_F1;
        InvE_F1_9 = InvE_F1_7*InvE_F1*InvE_F1;
        logFact_npr = log( (mu+E_Fnpr)*InvMagMnpr );
        logFact_1 = log( (mu+E_F1)*InvMagM1 );

        dPdBTemp += 1./(4.*pi2)*( mu*p_F-M2*log((mu+p_F)/M_e) )
                    + 1./(2.*pi2)*( npr*( mu*E_Fnpr - MagMnpr_2*logFact_npr )
                                          - ( mu*E_F1 - MagM1_2*logFact_1 ) )
                    + 1./(4.*pi2)*( mu*E_Fnpr - (M2+4.*eB*npr)*logFact_npr + mu*E_F1 - (M2+4.*eB)*logFact_1 )
                    + B_2*eB/(2.*pi2)*( eB*npr*mu*InvMagMnpr_2*InvE_Fnpr - 2.*logFact_npr
                                                  - eB*mu*InvMagM1_2*InvE_F1 + 2.*logFact_1 )
                    + B_4*eB2*eB*mu/(24.*pi2)*( ( 20.*mu2*(M2+eB*npr)*MagMnpr_2
                                                    - 3.*(4.*M2 + 3.*eB*npr)*MagMnpr_4
                                                    - 8.*mu4*(M2+eB*npr) )*InvMagMnpr_4*InvMagMnpr_2*InvE_Fnpr_5
                                                  - ( 20.*mu2*(M2+eB)*MagM1_2
                                                    - 3.*(4.*M2 + 3.*eB)*MagM1_4
                                                    - 8.*mu4*(M2+eB) )*InvMagM1_4*InvMagM1_2*InvE_F1_5 )
                    + mu*T*T/6.*( 0.5/p_F + ( npr*InvE_Fnpr - InvE_F1 ) + 0.5*( (p_F2 - eB*npr)*InvE_Fnpr_3 + (p_F2 - eB)*InvE_F1_3 )
                                                 + B_2*eB/2.*( (2.*p_F2 - eB*npr)*InvE_Fnpr_5 - (2.*p_F2 - eB)*InvE_F1_5 ) );
        d2PdB2Temp += 1./(pi2)*( -npr*npr*logFact_npr + logFact_1 )
                             + 1./(2.*pi2)*( npr*( eB*npr*mu*InvMagMnpr_2*InvE_Fnpr - 2.*logFact_npr )
                                                  + ( eB*mu*InvMagM1_2*InvE_F1 - 2.*logFact_1 ) )
                             + B_2/(2.*pi2)*( -eB*npr*mu*( 4.*M4 + M2*(13.*eB*npr - 4.*mu2) + 2.*eB*npr*(5.*eB*npr - 3.*mu2) )*InvMagMnpr_4*InvE_Fnpr_3
                                                      - 2.*logFact_npr
                                                + eB*mu*( 4.*M4 + M2*(13.*eB - 4.*mu2) + 2.*eB*(5.*eB - 3.*mu2) )*InvMagM1_4*InvE_F1_3
                                                      + 2.*logFact_1 )
                             + B_4*eB2*mu/(24.*pi2)*( -( 2.*mu2*MagMnpr_4*( 48.*M4 + 52.*M2*eB*npr + 17.*eB2*npr*npr )
                                                                             - 3.*MagMnpr_4*MagMnpr_2*( 12.*M4 + 8.*M2*eB*npr + 3.*eB2*npr*npr )
                                                                             + 8.*mu4*mu2*(3.*M4 + 4.*M2*eB*npr + 2*eB2*npr*npr)
                                                                             - 28.*mu4*MagMnpr_2*(3.*M4 + 4.*M2*eB*npr + 2.*eB2*npr*npr) )*InvMagMnpr_4*InvMagMnpr_4*InvE_Fnpr_7
                                                                      + ( 2.*mu2*MagM1_4*(48.*M4 + 52.*M2*eB + 17.*eB2)
                                                                             - 3.*MagM1_4*MagM1_2*(12.*M4 + 8.*M2*eB + 3.*eB2)
                                                                             + 8.*mu4*mu2*(3.*M4 + 4.*M2*eB + 2.*eB2)
                                                                             - 28.*mu4*MagM1_2*(3.*M4 + 4.*M2*eB + 2.*eB2) )*InvMagM1_4*InvMagM1_4*InvE_F1_7 )
                             + T*T*mu/6.*( (npr*npr*InvE_Fnpr_3 - InvE_F1_3) + 0.5*( (2.*p_F2-eB*npr)*InvE_Fnpr_5*npr + (2.*p_F2-eB)*InvE_F1_5 )
                                                          + B_2/2.*( ( 2.*M4 - eB2*npr*npr + 4.*eB*npr*mu2 + 2.*mu4 - 4.*M2*(eB*npr + mu2) )*InvE_Fnpr_7
                                                                   - ( 2.*M4 - eB2 + 4.*eB*mu2 + 2.*mu4 - 4.*M2*(eB + mu2) )*InvE_F1_7 ) );
        d2PdBdmuTemp += p_F/(2.*pi2) + 1./(pi2)*( npr*E_Fnpr - E_F1 )
                        + 1./(2.*pi2)*( ( p_F2-3.*eB*npr )*InvE_Fnpr + ( p_F2-3.*eB )*InvE_F1 )
                        + B_2*eB/(2.*pi2)*( -( 2.*p_F2 - 3.*eB*npr )*InvE_Fnpr_3 + ( 2.*p_F2 - 3.*eB )*InvE_F1_3 )
                        + B_4*eB2*eB/(8.*pi2)*( -( 4.*p_F2 - 3.*eB*npr )*InvE_Fnpr_7 + ( 4.*p_F2 - 3.*eB )*InvE_F1_7 )
                        - T*T/6.*( 0.5*M2*p_F2/pow(p_F2-M2,2.5) + npr*MagMnpr*InvE_Fnpr_3 - MagM1*InvE_F1_3
                                  + 0.5*( (mu2*(M2 + 4.*eB*npr) - MagMnpr_2*(M2 + eB*npr))*InvE_Fnpr_5 + (mu2*(M2 + 4.*eB) - MagM1_2*(M2 + eB))*InvE_F1_5 )
                                  + B_2*eB/2.*( ( 4.*mu4 - 2.*mu2*(M2 - 4.*eB*npr) - 2.*eB2*npr*npr - 5.*eB*npr*M2 - 2.*M4 )*InvE_Fnpr_7
                                              - ( 4.*mu4 - 2.*mu2*(M2 - 4.*eB) - 2.*eB2 - 5.*eB*M2 - 2.*M4 )*InvE_F1_7 ) );
        d2PdT2Temp += eB*mu/(6.*p_F) + eB*mu*( 1./3.*( -1./eB*( E_Fnpr - E_F1 ) + 0.5*( InvE_Fnpr + InvE_F1 )
                                                + 0.5*B_2*eB*( InvE_Fnpr_3 - InvE_F1_3 ) + 0.625*B_4*eB2*eB*( InvE_Fnpr_7 - InvE_F1_7 ) )
                                                + 7.*pi2*T*T/10.*( M2/(2.*p_F2*p_F2*p_F) -1./(3.*eB)*( (2.*mu2 - 3.*M2 - 6.*eB*npr)*InvE_Fnpr*InvE_Fnpr*InvE_Fnpr - (2.*mu2 - 3.*M2 - 6.*eB)*InvE_F1_3 )
                                                                        + 0.5*( (M2 + 2.*eB*npr)*InvE_Fnpr_5 + (M2 + 2.*eB)*InvE_F1_5 )
                                                                        + 0.5*B_2*eB*( (2.*mu2 + 3.*M2 + 6.*eB*npr)*InvE_Fnpr_7 - (2.*mu2 + 3.*M2 + 6.*eB)*InvE_F1_7 )
                                                                        + 4.375*B_4*eB2*eB*( (2.*mu2 + M2 + 2.*eB*npr)*InvE_Fnpr_9*InvE_Fnpr*InvE_Fnpr - (2.*mu2 + M2 + 2.*eB)*InvE_F1_9*InvE_F1*InvE_F1 ) ) );
        d2PdTdBTemp += T/6.*mu/p_F + T*mu/3.*( npr*InvE_Fnpr - InvE_F1
                                        + 0.5*( (E_Fnpr*E_Fnpr + eB*npr)*InvE_Fnpr*InvE_Fnpr*InvE_Fnpr + (E_F1*E_F1 + eB)*InvE_F1*InvE_F1*InvE_F1 )
                                        + 0.5*B_2*eB*( (2.*E_Fnpr*E_Fnpr + 3.*eB*npr)*InvE_Fnpr_5 - (2.*E_F1*E_F1 + 3.*eB)*InvE_F1_5 )
                                        + 0.625*B_4*eB2*eB*( (4.*E_Fnpr + 7.*eB*npr)*InvE_Fnpr_9 - (4.*E_F1 +7.*eB)*InvE_F1_9 ) );
    }
    else if( n_max > 0.5 ){
        //dPdB, d2PdB2, d2PdBdmu, d2PdT2, d2PdTdB contribution from n = 0
        dPdBTemp += 1./(2.*pi2)*( 0.5*mu*p_F - 0.5*M2*log( (mu+p_F)/M_e ) + pi2/6.*mu/p_F*T*T + 7.*pi4/120.*mu*M2/(p_F2*p_F2*p_F)*T*T*T*T );
        d2PdB2Temp += 0.;
        d2PdBdmuTemp += 1./(2.*pi2)*( p_F - pi2/6.*M2/(p_F2*p_F)*T*T - 7.*pi4/120.*M2*(4.*mu2 + M2)/(p_F2*p_F2*p_F2*p_F)*T*T*T*T );
        d2PdT2Temp += eB*mu/(6.*p_F);
        d2PdTdBTemp += T/6.*mu/p_F;
    }
    else{
        //All contributions to dPdB, d2PdB2, d2PdBdmu, d2PdT2, d2PdTdB when n_max = 0 come from "Omega_xyLowT_VaryingT_MaxLL"
        dPdBTemp += 0.;
        d2PdB2Temp += 0.;
        d2PdBdmuTemp += 0.;
        d2PdT2Temp += 0.;
        d2PdTdBTemp += 0.;
    }

    return;
}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying.
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
    Output: dPdBTemp: dP_e/dB in MeV^2
            d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
            d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
            d2PdT2Temp: d^2P/dT^2 in MeV^2
            d2PdTdBTemp: d^3P/dTdB in MeV
*/
void Omega_xyLowT_VaryingT_MaxLL(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp)
{
    static double pi2 = pi*pi;
    static double pi4 = pi2*pi2;
    static double M2 = M_e*M_e;
    double M_n, M_np1, a, a1, b, pref;
    double G1_ab, G2_ab, H2_ab, I1_ab, I2_ab;
    double G1_ab1, G2_ab1, H2_ab1, I1_ab1, I2_ab1;
    double LiN1_2, Li1_2, Li3_2, Li5_2, Li7_2; //temporary variables to hold interpolated polylogarithm function evaluations
    double KroneckerDelta = 0; //Kronecker delta function for n_max, 0; returns 1 if n_max = 0 and zero otherwise. Assume 0 and then correct otherwise.

    M_n = sqrt(M2 + 2.*eB*n_max); //M_{n_max}
    M_np1 = sqrt(M2 + 2.*eB*(n_max+1)); //M_{n_max+1}

    // n = n_max contribution

    a = M_n/T;
    b = mu/T;

    if( b - a > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab = b*sqrt( b*b/(a*a)-1. )/(2.*a) + 0.5*log( sqrt( b*b/(a*a)-1. ) + b/a ) + pi2*b*(b*b-2.*a*a)/(6.*a*a*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(4.*a*a + b*b)/pow(b*b-a*a,3.5);
        G2_ab = log( sqrt( b*b/(a*a)-1. ) + b/a ) - pi2*b/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        H2_ab = -b/(a*sqrt( b*b/(a*a)-1. )) - pi2*b*a*a/(2.*pow(b*b-a*a,2.5)) - 7.*pi4/24.*b*a*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,4.5);
        I1_ab = sqrt( b*b/(a*a)-1. ) - pi2*a/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        I2_ab = -1./sqrt( b*b/(a*a)-1. ) - pi2*a*(2.*b*b + a*a)/(6.*pow(b*b-a*a,2.5)) - 7.*pi4/120.*a*(8.*b*b*b*b + 24.*a*a*b*b + 3.*a*a*a*a)/pow(b*b-a*a,4.5);
    }
    else if(a - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a))*exp(b-a);
        G1_ab = pref*( 1. + 7./(8.*a) + 57./(128.*a*a) - 195./(1024.*a*a*a) );
        G2_ab = pref*( 1. - 1./(8.*a) + 9./(128.*a*a) - 75./(1024.*a*a*a) );
        H2_ab = -pref*( a + 3./8. - 15./(128.*a) + 105./(1024.*a*a) );
        I1_ab = pref*( 1. + 3./(8.*a) - 15./(128.*a*a) + 105./(1024.*a*a*a) );
        I2_ab = -pref*( a - 1./8. + 9./(128.*a) - 75./(1024.*a*a) );
    }
    else if( sqrt(a*a+b*b) < 50. ){ //full interpolating functions
        G1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else{ //polylogarithms with argument -exp(b-a), using interpolating functions for speed
        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, b-a, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, b-a, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, b-a, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, b-a, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, b-a, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a));
        G1_ab = -pref*( Li1_2 + 7./(8.*a)*Li3_2 + 57./(128.*a*a)*Li5_2 - 195./(1024.*a*a*a)*Li7_2 );
        G2_ab = -pref*( Li1_2 - 1./(8.*a)*Li3_2 + 9./(128.*a*a)*Li5_2 - 75./(1024.*a*a*a)*Li7_2 );
        H2_ab = pref*( a*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a)*Li3_2 + 105./(1024.*a*a)*Li5_2 );
        I1_ab = -pref*( Li1_2 + 3./(8.*a)*Li3_2 - 15./(128.*a*a)*Li5_2 + 105./(1024.*a*a*a)*Li7_2 );
        I2_ab = pref*( a*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a)*Li3_2 - 75./(1024.*a*a)*Li5_2 );

    }

    // n = n_max + 1 contribution

    a1 = M_np1/T;

    if( b - a1 > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab1 = b*sqrt( b*b/(a1*a1)-1. )/(2.*a1) + 0.5*log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) + pi2*b*(b*b-2.*a1*a1)/(6.*a1*a1*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(4.*a1*a1 + b*b)/pow(b*b-a1*a1,3.5);
        G2_ab1 = log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) - pi2*b/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        H2_ab1 = -b/(a1*sqrt( b*b/(a1*a1)-1. )) - pi2*b*a1*a1/(2.*pow(b*b-a1*a1,2.5)) - 7.*pi4/24.*b*a1*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,4.5);
        I1_ab1 = sqrt( b*b/(a1*a1)-1. ) - pi2*a1/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        I2_ab1 = -1./sqrt( b*b/(a1*a1)-1. ) - pi2*a1*(2.*b*b + a1*a1)/(6.*pow(b*b-a1*a1,2.5)) - 7.*pi4/120.*a1*(8.*b*b*b*b + 24.*a1*a1*b*b + 3.*a1*a1*a1*a1)/pow(b*b-a1*a1,4.5);
    }
    else if(a1 - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a1))*exp(b-a1);
        G1_ab1 = pref*( 1. + 7./(8.*a1) + 57./(128.*a1*a1) - 195./(1024.*a1*a1*a1) );
        G2_ab1 = pref*( 1. - 1./(8.*a1) + 9./(128.*a1*a1) - 75./(1024.*a1*a1*a1) );
        H2_ab1 = -pref*( a1 + 3./8. - 15./(128.*a1) + 105./(1024.*a1*a1) );
        I1_ab1 = pref*( 1. + 3./(8.*a1) - 15./(128.*a1*a1) + 105./(1024.*a1*a1*a1) );
        I2_ab1 = -pref*( a1 - 1./8. + 9./(128.*a1) - 75./(1024.*a1*a1) );
    }
    else if( sqrt(a1*a1+b*b) < 50. ){ //full interpolating functions
        G1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a1, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a1, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a1, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a1, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a1, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else{ //polylogarithms with argument -exp(b-a), using interpolating functions for speed

        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, b-a1, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, b-a1, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, b-a1, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, b-a1, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, b-a1, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a1));
        G1_ab1 = -pref*( Li1_2 + 7./(8.*a1)*Li3_2 + 57./(128.*a1*a1)*Li5_2 - 195./(1024.*a1*a1*a1)*Li7_2 );
        G2_ab1 = -pref*( Li1_2 - 1./(8.*a1)*Li3_2 + 9./(128.*a1*a1)*Li5_2 - 75./(1024.*a1*a1*a1)*Li7_2 );
        H2_ab1 = pref*( a1*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a1)*Li3_2 + 105./(1024.*a1*a1)*Li5_2 );
        I1_ab1 = -pref*( Li1_2 + 3./(8.*a1)*Li3_2 - 15./(128.*a1*a1)*Li5_2 + 105./(1024.*a1*a1*a1)*Li7_2 );
        I2_ab1 = pref*( a1*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a1)*Li3_2 - 75./(1024.*a1*a1)*Li5_2 );
    }

    if(n_max < 0.5) KroneckerDelta = 1.;

    dPdBTemp += (2. - KroneckerDelta)/(2.*pi2)*( M_n*M_n*G1_ab - (M_n*M_n + eB*n_max)*G2_ab ) + 1./pi2*( M_np1*M_np1*G1_ab1 - (M_np1*M_np1 + eB*(n_max+1.))*G2_ab1 );
    d2PdB2Temp += -n_max/(2.*pi2)*(2. - KroneckerDelta)*( 2.*G2_ab + eB*n_max/(M_n*M_n)*H2_ab ) - (n_max+1.)/pi2*( 2.*G2_ab1 + eB*(n_max+1.)/(M_np1*M_np1)*H2_ab1 );
    d2PdBdmuTemp += (2. - KroneckerDelta)/(2.*pi2)*M_n*( I1_ab + eB*n_max/(M_n*M_n)*I2_ab ) + 1./pi2*M_np1*( I1_ab1 + eB*(n_max+1.)/(M_np1*M_np1)*I2_ab1 );
    d2PdT2Temp += 1./(2.*pi2*T*T)*eB*( (2. - KroneckerDelta)*M_n*M_n*( 2.*G1_ab - G2_ab - ( b*b/(a*a) + 1. )*H2_ab + 2.*b/a*(I2_ab - I1_ab) )
                                        + 2.*M_np1*M_np1*( 2.*G1_ab1 - G2_ab1 - ( b*b/(a1*a1) + 1. )*H2_ab1 + 2.*b/a1*(I2_ab1 - I1_ab1) ) );
    d2PdTdBTemp += 1./(2.*pi2*T)*( (2. - KroneckerDelta)*( M_n*M_n*(2.*G1_ab - G2_ab - b/a*I1_ab) + eB*n_max*( H2_ab - b/a*I2_ab) )
                                    + 2.*( M_np1*M_np1*(2.*G1_ab1 - G2_ab1 - b/a1*I1_ab1) + eB*(n_max+1.)*( H2_ab1 - b/a1*I2_ab1 ) ) );

    return;

}

/*
    High-temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying.
    Consists of regular (non-oscillatory) and oscillatory terms with finite T thermal damping. Derived using the
    Poisson summation formula, and based on the Lifshitz-Kosevich formula.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Output: dPdBTemp: dP_e/dB in MeV^2
            d2PdB2Temp: d^2P_e/dB^2 (dimensionless)
            d2PdBdmuTemp: d^2P_e/dB/dmu_e in MeV
            d2PdT2Temp: d^2P/dT^2 in MeV^2
            d2PdTdBTemp: d^3P/dTdB in MeV
*/
void Omega_xyHighT_VaryingT(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp)
{

    //T_B = 100*M*( sqrt(1+2*(n_max+1)*eB/eB_crit) - sqrt(1+2*n_max*eB/eB_crit) ); //critical temperature in MeV

    static double zeta3_2 = 2.61237534869; //zeta(3/2)
    static double sqrtpi = sqrt(pi);
    static double pi2 = pi*pi; //pi^2
    static double pi2_5 = pi2*sqrt(pi); //pi^2.5
    static double pi3 = pi2*pi; //pi^3
    static double sqrt2opi = sqrt(2./pi); //sqrt(2/pi)
    static double M2 = M_e*M_e;
    static int p_max = 5; //maximum number of terms to include in high-T approximation sum

    double n, sqrtn, lambda, R_T, R_T2, R_T3, Sin, Cos, phase, sinhlambda;
    double i1val, i2val, i3val, h1val, h2val;

    double sqrteB = sqrt(eB);
    double p_F2 = p_F*p_F;

    // First compute regular (non-oscillatory part) of the terms. Add them to zeroed variables.
    // double x = p_F2/eB;
    // if( x < 1e4 ){
    //     i1val = gsl_spline_eval( OmB_Interpolators.i_1, x, OmB_Interpolators.i_1_acc );
    //     i2val = gsl_spline_eval( OmB_Interpolators.i_2, x, OmB_Interpolators.i_2_acc );
    //     i3val = gsl_spline_eval( OmB_Interpolators.i_3, x, OmB_Interpolators.i_3_acc );
    // }
    // else{
    //     i1val = sqrtpi/sqrt(x)/3.;
    //     i2val = sqrtpi/(6.*x*sqrt(x));
    //     i3val = sqrtpi/(4.*x*x*sqrt(x));
    // }
    //
    // h1val = gsl_spline2d_eval( OmB_Interpolators.h_1, mu/M_e, eB/(M2), OmB_Interpolators.h_1_acc_a, OmB_Interpolators.h_1_acc_b );
    // h2val = gsl_spline2d_eval( OmB_Interpolators.h_2, mu/M_e, eB/(M2), OmB_Interpolators.h_2_acc_a, OmB_Interpolators.h_2_acc_b );
    //
    // dPdBTemp += sqrteB/(2.*pi2_5)*h1val - M2/(4.*pi2_5*sqrteB)*h2val - 1./(8.*pi2_5)*sqrteB*( mu*i1val - sqrt2opi*M_e*zeta3_2 );
    // d2PdB2Temp += 1./(2.*pi2_5*sqrteB)*h1val - M2/(2.*pi2_5*eB*sqrteB)*h2val - 5./(16.*pi2_5*sqrteB)*( mu*i1val - sqrt2opi*M_e*zeta3_2 );
    // d2PdBdmuTemp += 3.*sqrteB/(8.*pi2_5)*i1val + p_F2/(4.*pi2_5*sqrteB)*i2val;
    // d2PdT2Temp += mu*p_F/3. - mu*sqrt(eB/pi)/6.*i2val;
    // d2PdTdBTemp += -mu*T/(12.*sqrt(pi*eB))*( i2val + 2.*p_F2/eB*i3val );

    // Next compute oscillatory part of the terms, involving a sum over p=0,...,p_max
    for(int p = 1; p <= p_max; p++){

        n = double(p);
        sqrtn = sqrt(n);
        lambda = 2.*pi*pi*mu*T*n/eB; //lambda factor as defined in notes times pi
        if(lambda < 20.){
            sinhlambda = sinh(lambda);
            R_T = lambda/sinhlambda;
            R_T2 = ( lambda/tanh(lambda) - 1. )/sinhlambda;
            R_T3 = 0.5*lambda*( 3.*lambda + lambda*cosh(2.*lambda) - 2.*sinh(2.*lambda) )/(sinhlambda*sinhlambda*sinhlambda);
        }
        else{ //use large lambda expansion for hyperbolic functions for lambda > 20
            R_T = 2.*lambda*exp(-lambda);
            R_T2 = 2.*( lambda - 1. )*exp(-lambda);
            R_T3 = 4.*lambda*( 3.*lambda*exp(-3.*lambda) + 0.5*lambda*exp(-lambda) - exp(-lambda) );
        }

        phase = std::remainder( (0.25 - n/eB*p_F2)*pi, 2.*pi );
        Sin = sin( phase );
        Cos = cos( phase );

        dPdBTemp += R_T*( sqrteB*p_F2/(4.*pi3*mu*n*sqrtn)*Sin - 5.*eB*sqrteB/(8.*pi3*pi*mu*n*n*sqrtn)*Cos );
        d2PdB2Temp += R_T*( p_F2*p_F2/(4.*pi2*mu*eB*sqrteB*sqrtn)*Cos + 3.*p_F2/(4.*pi3*mu*sqrteB*n*sqrtn)*Sin );
        d2PdBdmuTemp += -R_T*( 3.*sqrteB/(4*pi3*n*sqrtn)*Sin + p_F2/(2.*pi2*sqrteB*sqrtn)*Cos );
        d2PdT2Temp += -R_T3*eB/n*sqrteB/sqrtn/(2.*pi2*T)*Cos;
        d2PdTdBTemp += R_T2/(2.*pi*sqrteB*sqrtn)*( 3.*eB/(2.*pi*n)*Cos - mu*mu*Sin );
    }

    return;
}

/*
    Computes thermodynamic derivatives Omega_xy to compute only the specific heat capacity of Landau-quantized electrons where T is varying.
    Fundamental variables of Omega are B, mu_e, T.
    Inputs: Bmag: magnitude of magnetic field in reduced units
            mu: electron chemical potential in MeV. Taken to vary only in radial direction
            T: temperature in reduced units
            N_GC: number of ghost cells to exchange
            n_t: lapse function
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Outputs: c_v_eB: specific heat capacity of magnetized electrons c_v = T*ds/dT = -T*d^2Omega_e/dT^2 = T*d^2P_e/dT^2 in reduced units.
*/
void Omega_xyVaryingTC_veOnly(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & c_v_eB)
{
    static double M2 = M_e*M_e;
    static double eB_pref = unit_e*B_0*statCGtoMeV2; //constant prefactor for computing e*B in MeV^2
    double eB, n_max, thermalBFactor, p_F; //fundamental charge times B in MeV^2, maximum occupied Landau level at T=0, dimensionless thermal damping factor for de Haas-van Alphen oscillations, Fermi momentum in MeV

    double d2PdT2Temp;

    for(size_t i=0; i<Bmag.shape()[0]-N_GC; i++){
        p_F = sqrt( mu[i]*mu[i]-M2 ); //Fermi momentum in MeV
        for(size_t j=0; j<Bmag.shape()[1]; j++){
            if( Bmag[i][j] < 1e-5 ){
                c_v_eB[i][j] = ( k_B*T[i][j]/n_t[i]*T_0*mu[i]*p_F/3. )*k_B*MeVtoErg/(hbarc*hbarc*hbarc)*1e39/s_0;
            }
            else{
                eB = Bmag[i][j]*eB_pref; //e*B in MeV^2
                n_max = std::floor( (mu[i]*mu[i]-M2)/(2.*eB) ); //highest occupied Landau level at T=0
                thermalBFactor = 2.*pi*pi*mu[i]*k_B*T[i][j]/n_t[i]*T_0/eB;
                //T_B = M_E*( sqrt(1. + 2.*(n_max+1.)*eB/eB_crit) - sqrt(1. + 2.*n_max*eB/eB_crit) );
                d2PdT2Temp = 0.;

                if(thermalBFactor > 1.){
                    Omega_xyHighT_VaryingTC_veOnly(eB, mu[i], k_B*T[i][j]/n_t[i]*T_0, p_F, OmB_Interpolators, d2PdT2Temp);
                }
                else{
                    Omega_xyLowT_VaryingT_EMSumsC_veOnly(eB, mu[i], k_B*T[i][j]/n_t[i]*T_0, p_F, n_max, d2PdT2Temp);
                    //if( std::isnan(d2PdT2Temp) == true ) std::cout << "isnan in Low T approx. Euler-Maclaurin sums" << std::endl;
                    Omega_xyLowT_VaryingT_MaxLLC_veOnly(eB, mu[i], k_B*T[i][j]/n_t[i]*T_0, n_max, OmB_Interpolators, d2PdT2Temp);
                    //if( std::isnan(d2PdT2Temp) == true ) std::cout << "isnan in Low T approx. Max Landau level" << std::endl;
                }

                c_v_eB[i][j] = ( k_B*T[i][j]/n_t[i]*T_0*d2PdT2Temp )*k_B*MeVtoErg/(hbarc*hbarc*hbarc)*1e39/s_0; //in MeV^3 -> MeV^4/K -> erg/(fm^3 K) -> erg/(cm^3 K) -> reduced units
            }
        }
    }

    return;

}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying. Only computes d2PdT2 i.e., Landau-quantized electron specific heat capacity
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Output: d2PdT2Temp: d^2P/dT^2 in MeV^2
*/
void Omega_xyLowT_VaryingT_EMSumsC_veOnly(double eB, double mu, double T, double p_F, double n_max, double & d2PdT2Temp)
{
    //double T_B = M_e*( sqrt(1. + 2.*(n_max+1.)*eB/eB_crit) - sqrt(1. + 2.*n_max*eB/eB_crit) ); //critical temperature in MeV

    static double B_2 = 1./6., B_4 = -1./30.;
    double npr = n_max - 1.;
    double E_Fnpr = sqrt( p_F*p_F - 2.*eB*npr );
    double InvE_Fnpr = 1./E_Fnpr;
    double E_F1 = sqrt( p_F*p_F - 2.*eB );
    double InvE_F1 = 1./E_F1;
    static double M2 = M_e*M_e;
    static double pi2 = pi*pi;
    double mu2 = mu*mu;
    double eB2 = eB*eB;
    double pF2 = p_F*p_F;

    double InvE_Fnpr_3, InvE_Fnpr_5, InvE_Fnpr_7, InvE_Fnpr_9, InvE_F1_3, InvE_F1_5, InvE_F1_7, InvE_F1_9;

    if( n_max > 1.1 ){
        InvE_Fnpr_3 = InvE_Fnpr*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_5 = InvE_Fnpr_3*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_7 = InvE_Fnpr_5*InvE_Fnpr*InvE_Fnpr;
        InvE_Fnpr_9 = InvE_Fnpr_7*InvE_Fnpr*InvE_Fnpr;
        InvE_F1_3 = InvE_F1*InvE_F1*InvE_F1;
        InvE_F1_5 = InvE_F1_3*InvE_F1*InvE_F1;
        InvE_F1_7 = InvE_F1_5*InvE_F1*InvE_F1;
        InvE_F1_9 = InvE_F1_7*InvE_F1*InvE_F1;

        d2PdT2Temp += eB*mu/(6.*p_F) + eB*mu*( 1./3.*( -1./eB*( E_Fnpr - E_F1 ) + 0.5*( InvE_Fnpr + InvE_F1 )
                                                + 0.5*B_2*eB*( InvE_Fnpr_3 - InvE_F1_3 ) + 0.625*B_4*eB2*eB*( InvE_Fnpr_7 - InvE_F1_7 ) )
                                                + 7.*pi2*T*T/10.*( M2/(2.*pF2*pF2*p_F)  -1./(3.*eB)*( (2.*mu2 - 3.*M2 - 6.*eB*npr)*InvE_Fnpr*InvE_Fnpr*InvE_Fnpr - (2.*mu2 - 3.*M2 - 6.*eB)*InvE_F1_3 )
                                                                        + 0.5*( (M2 + 2.*eB*npr)*InvE_Fnpr_5 + (M2 + 2.*eB)*InvE_F1_5 )
                                                                        + 0.5*B_2*eB*( (2.*mu2 + 3.*M2 + 6.*eB*npr)*InvE_Fnpr_7 - (2.*mu2 + 3.*M2 + 6.*eB)*InvE_F1_7 )
                                                                        + 4.375*B_4*eB2*eB*( (2.*mu2 + M2 + 2.*eB*npr)*InvE_Fnpr_9*InvE_Fnpr*InvE_Fnpr - (2.*mu2 + M2 + 2.*eB)*InvE_F1_9*InvE_F1*InvE_F1 ) ) );
    }
    else if( n_max > 0.5 ){
        //d2PdT2 contribution from n = 0
        d2PdT2Temp += eB*mu/(6.*p_F);
    }
    else{
        //All contributions to d2PdT2 when n_max = 0 come from "Omega_xyLowT_VaryingT_MaxLL"
        d2PdT2Temp += 0.;
    }

    return;
}

/*
    Low temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying. Only computes d2PdT2 i.e., Landau-quantized electron specific heat capacity
    Euler-Maclaurin sum terms for n < n_max with Sommerfeld expansion finite T thermal corrections. The highest-occupied Landau level is treated separately from this term.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            n_max: maximum occupied Landau level at T=0
    Output: d2PdT2Temp: d^2P/dT^2 in MeV^2
*/
void Omega_xyLowT_VaryingT_MaxLLC_veOnly(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & d2PdT2Temp)
{
    static double pi2 = pi*pi;
    static double pi4 = pi2*pi2;
    static double M2 = M_e*M_e;
    double M_n, M_np1, a, a1, b, pref;
    double G1_ab, G2_ab, H2_ab, I1_ab, I2_ab;
    double G1_ab1, G2_ab1, H2_ab1, I1_ab1, I2_ab1;
    double LiN1_2, Li1_2, Li3_2, Li5_2, Li7_2; //temporary variables to hold interpolated polylogarithm function evaluations
    double KroneckerDelta = 0; //Kronecker delta function for n_max, 0; returns 1 if n_max = 0 and zero otherwise. Assume 0 and then correct otherwise.

    M_n = sqrt(M2 + 2.*eB*n_max); //M_{n_max}
    M_np1 = sqrt(M2 + 2.*eB*(n_max+1)); //M_{n_max+1}

    // n = n_max contribution

    a = M_n/T;
    b = mu/T;

    if( b - a > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab = b*sqrt( b*b/(a*a)-1. )/(2.*a) + 0.5*log( sqrt( b*b/(a*a)-1. ) + b/a ) + pi*pi*b*(b*b-2.*a*a)/(6.*a*a*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(4.*a*a + b*b)/pow(b*b-a*a,3.5);
        G2_ab = log( sqrt( b*b/(a*a)-1. ) + b/a ) - pi*pi*b/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        H2_ab = -b/(a*sqrt( b*b/(a*a)-1. )) - pi*pi*b*a*a/(2.*pow(b*b-a*a,2.5)) - 7.*pi4/24.*b*a*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,4.5);
        I1_ab = sqrt( b*b/(a*a)-1. ) - pi*pi*a/(6.*pow(b*b-a*a,1.5)) - 7.*pi4/120.*a*(4.*b*b + 3.*a*a)/pow(b*b-a*a,3.5);
        I2_ab = -1./sqrt( b*b/(a*a)-1. ) - pi*pi*a*(2.*b*b + a*a)/(6.*pow(b*b-a*a,2.5)) - 7.*pi4/120.*a*(8.*b*b*b*b + 24.*a*a*b*b + 3.*a*a*a*a)/pow(b*b-a*a,4.5);
    }
    else if(a - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a))*exp(b-a);
        G1_ab = pref*( 1. + 7./(8.*a) + 57./(128.*a*a) - 195./(1024.*a*a*a) );
        G2_ab = pref*( 1. - 1./(8.*a) + 9./(128.*a*a) - 75./(1024.*a*a*a) );
        H2_ab = -pref*( a + 3./8. - 15./(128.*a) + 105./(1024.*a*a) );
        I1_ab = pref*( 1. + 3./(8.*a) - 15./(128.*a*a) + 105./(1024.*a*a*a) );
        I2_ab = -pref*( a - 1./8. + 9./(128.*a) - 75./(1024.*a*a) );
    }
    else if( sqrt(a*a+b*b) < 50. ){ //full interpolating functions
        G1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else{ //polylogarithms with argument -exp(b-a), using interpolating functions for speed
        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, b-a, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, b-a, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, b-a, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, b-a, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, b-a, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a));
        G1_ab = -pref*( Li1_2 + 7./(8.*a)*Li3_2 + 57./(128.*a*a)*Li5_2 - 195./(1024.*a*a*a)*Li7_2 );
        G2_ab = -pref*( Li1_2 - 1./(8.*a)*Li3_2 + 9./(128.*a*a)*Li5_2 - 75./(1024.*a*a*a)*Li7_2 );
        H2_ab = pref*( a*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a)*Li3_2 + 105./(1024.*a*a)*Li5_2 );
        I1_ab = -pref*( Li1_2 + 3./(8.*a)*Li3_2 - 15./(128.*a*a)*Li5_2 + 105./(1024.*a*a*a)*Li7_2 );
        I2_ab = pref*( a*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a)*Li3_2 - 75./(1024.*a*a)*Li5_2 );

    }

    // n = n_max + 1 contribution

    a1 = M_np1/T;

    if( b - a1 > 50. ){ //T=0 plus Sommerfeld expansion correction
        G1_ab1 = b*sqrt( b*b/(a1*a1)-1. )/(2.*a1) + 0.5*log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) + pi*pi*b*(b*b-2.*a1*a1)/(6.*a1*a1*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(4.*a1*a1 + b*b)/pow(b*b-a1*a1,3.5);
        G2_ab1 = log( sqrt( b*b/(a1*a1)-1. ) + b/a1 ) - pi*pi*b/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*b*(2.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        H2_ab1 = -b/(a1*sqrt( b*b/(a1*a1)-1. )) - pi*pi*b*a1*a1/(2.*pow(b*b-a1*a1,2.5)) - 7.*pi4/24.*b*a1*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,4.5);
        I1_ab1 = sqrt( b*b/(a1*a1)-1. ) - pi*pi*a1/(6.*pow(b*b-a1*a1,1.5)) - 7.*pi4/120.*a1*(4.*b*b + 3.*a1*a1)/pow(b*b-a1*a1,3.5);
        I2_ab1 = -1./sqrt( b*b/(a1*a1)-1. ) - pi*pi*a1*(2.*b*b + a1*a1)/(6.*pow(b*b-a1*a1,2.5)) - 7.*pi4/120.*a1*(8.*b*b*b*b + 24.*a1*a1*b*b + 3.*a1*a1*a1*a1)/pow(b*b-a1*a1,4.5);
    }
    else if(a1 - b > 50.){ //Taylor series expansion of polylogarithms
        pref = sqrt(pi/(2.*a1))*exp(b-a1);
        G1_ab1 = pref*( 1. + 7./(8.*a1) + 57./(128.*a1*a1) - 195./(1024.*a1*a1*a1) );
        G2_ab1 = pref*( 1. - 1./(8.*a1) + 9./(128.*a1*a1) - 75./(1024.*a1*a1*a1) );
        H2_ab1 = -pref*( a1 + 3./8. - 15./(128.*a1) + 105./(1024.*a1*a1) );
        I1_ab1 = pref*( 1. + 3./(8.*a1) - 15./(128.*a1*a1) + 105./(1024.*a1*a1*a1) );
        I2_ab1 = -pref*( a1 - 1./8. + 9./(128.*a1) - 75./(1024.*a1*a1) );
    }
    else if( sqrt(a1*a1+b*b) < 50. ){ //full interpolating functions
        G1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_1, a1, b, OmB_Interpolators.LogG_1_acc_a, OmB_Interpolators.LogG_1_acc_b ) );
        G2_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogG_2, a1, b, OmB_Interpolators.LogG_2_acc_a, OmB_Interpolators.LogG_2_acc_b ) );
        H2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNH_2, a1, b, OmB_Interpolators.LogNH_2_acc_a, OmB_Interpolators.LogNH_2_acc_b ) );
        I1_ab1 = exp( gsl_spline2d_eval( OmB_Interpolators.LogI_1, a1, b, OmB_Interpolators.LogI_1_acc_a, OmB_Interpolators.LogI_1_acc_b ) );
        I2_ab1 = -exp( gsl_spline2d_eval( OmB_Interpolators.LogNI_2, a1, b, OmB_Interpolators.LogNI_2_acc_a, OmB_Interpolators.LogNI_2_acc_b ) );
    }
    else{ //polylogarithms with argument -exp(b-a), using interpolating functions for speed

        LiN1_2 = gsl_spline_eval( OmB_Interpolators.PolyLogN1_2, b-a1, OmB_Interpolators.PolyLogN1_2_acc );
        Li1_2 = gsl_spline_eval( OmB_Interpolators.PolyLog1_2, b-a1, OmB_Interpolators.PolyLog1_2_acc );
        Li3_2 = gsl_spline_eval( OmB_Interpolators.PolyLog3_2, b-a1, OmB_Interpolators.PolyLog3_2_acc );
        Li5_2 = gsl_spline_eval( OmB_Interpolators.PolyLog5_2, b-a1, OmB_Interpolators.PolyLog5_2_acc );
        Li7_2 = gsl_spline_eval( OmB_Interpolators.PolyLog7_2, b-a1, OmB_Interpolators.PolyLog7_2_acc );

        pref = sqrt(pi/(2.*a1));
        G1_ab1 = -pref*( Li1_2 + 7./(8.*a1)*Li3_2 + 57./(128.*a1*a1)*Li5_2 - 195./(1024.*a1*a1*a1)*Li7_2 );
        G2_ab1 = -pref*( Li1_2 - 1./(8.*a1)*Li3_2 + 9./(128.*a1*a1)*Li5_2 - 75./(1024.*a1*a1*a1)*Li7_2 );
        H2_ab1 = pref*( a1*LiN1_2 + 3./8.*Li1_2 - 15./(128.*a1)*Li3_2 + 105./(1024.*a1*a1)*Li5_2 );
        I1_ab1 = -pref*( Li1_2 + 3./(8.*a1)*Li3_2 - 15./(128.*a1*a1)*Li5_2 + 105./(1024.*a1*a1*a1)*Li7_2 );
        I2_ab1 = pref*( a1*LiN1_2 - 1./8.*Li1_2 + 9./(128.*a1)*Li3_2 - 75./(1024.*a1*a1)*Li5_2 );
    }

    if(n_max < 0.5) KroneckerDelta = 1.;

    d2PdT2Temp += 1./(2.*pi*pi*T*T)*eB*( (2. - KroneckerDelta)*M_n*M_n*( 2.*G1_ab - G2_ab - ( b*b/(a*a) + 1. )*H2_ab + 2.*b/a*(I2_ab - I1_ab) )
                                        + 2.*M_np1*M_np1*( 2.*G1_ab1 - G2_ab1 - ( b*b/(a1*a1) + 1. )*H2_ab1 + 2.*b/a1*(I2_ab1 - I1_ab1) ) );

    return;

}

/*
    High-temperature approximation for thermodynamic derivatives Omega_xy = -P_xy for simulations where T is varying. Only computes d2PdT2 i.e., Landau-quantized electron specific heat capacity
    Consists of regular (non-oscillatory) and oscillatory terms with finite T thermal damping. Derived using the
    Poisson summation formula, and based on the Lifshitz-Kosevich formula.
    Inputs: eB: fundamental charge times magnetic field in MeV^2
            mu: electron chemical potential in MeV
            T: temperature in MeV
            p_F: electron Fermi momentum in MeV
            OmB_Interpolators: instance of class OmegaB_deriv_Interpolators used to store interpolating functions
    Output: d2PdT2Temp: d^2P/dT^2 in MeV^2
*/
void Omega_xyHighT_VaryingTC_veOnly(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & d2PdT2Temp)
{

    //T_B = 100*M*( sqrt(1+2*(n_max+1)*eB/eB_crit) - sqrt(1+2*n_max*eB/eB_crit) ); //critical temperature in MeV
    static double sqrtpi = sqrt(pi);
    static int p_max = 5; //maximum number of terms to include in high-T approximation sum

    double n, sqrtn, lambda, R_T3, Cos, phase;
    double i2val;

    double sqrteB = sqrt(eB);
    double p_F2 = p_F*p_F;

    // First compute regular (non-oscillatory part) of the terms. Add them to zeroed variables.
    double x = p_F2/eB;
    if( x < 1e4 ){
        i2val = gsl_spline_eval( OmB_Interpolators.i_2, x, OmB_Interpolators.i_2_acc );
    }
    else{
        i2val = sqrtpi/(6.*x*sqrt(x));
    }

    d2PdT2Temp += mu*p_F/3. - mu*sqrt(eB)/sqrtpi/6.*i2val;

    // Next compute oscillatory part of the terms, involving a sum over p=0,...,p_max
    for(int p = 1; p <= p_max; p++){

        n = double(p);
        sqrtn = sqrt(n);
        lambda = 2.*pi*pi*mu*T*n/eB; //lambda factor as defined in notes times pi
        if(lambda < 20.){
            R_T3 = 0.5*lambda*( 3.*lambda + lambda*cosh(2.*lambda) - 2.*sinh(2.*lambda) )/pow(sinh(lambda),3.);
        }
        else{ //use large lambda expansion for hyperbolic functions for lambda > 20
            R_T3 = 4.*lambda*( 3.*lambda*exp(-3.*lambda) + 0.5*lambda*exp(-lambda) - exp(-lambda) );
        }

        phase = std::remainder( (0.25 - n/eB*p_F2)*pi, 2.*pi );
        Cos = cos( phase );

        d2PdT2Temp += -R_T3*eB/n*sqrteB/sqrtn/(2.*pi*pi*T)*Cos;
    }

    return;
}
