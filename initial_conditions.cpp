#include <vector>
#include "boost/multi_array.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "common.h"
#include "initial_conditions.h"

/*
    Initializes the cell face-averaged magnetic field
    Inputs: r, theta: vectors of the cell mid-points in each direction
    bparams: BParams class object containing properties of the initial magnetic field
    dm: Domain class object containing the size of the full domain as Lr and Ltheta
    N_GC: number of ghost cells added in each direction
    Deltar, Deltatheta: grid spacing in each direction
    Output: B: 3D multi_array whose first dimension corresponds to the three field components and the final two dimensions are the grid points
*/
void InitializeB(std::vector<double> & r, std::vector<double> & theta, const BParams & bparams, const Domain & dm, size_t N_GC, std::vector<double> & Deltar, double Deltatheta, VectorField & B)
{
    double Ri = dm.r_min;
    double Ro = dm.r_max;
    double B_pol_pole = bparams.B_pol_init; //value of poloidal field at pole of star (r=Ro, theta=0) in reduced units
    double B_tor_max = bparams.B_tor_init; //maximum value of toroidal field in reduced units

    //Solve for mu_tilde parameter in initial poloidal field as a function of Ri/Ro. This is required for the initial field to satisfy the
    //outer boundary condition
    double mu_tilde;

    int status;
    int iter = 0;
    int max_iter = 100;

    const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	gsl_function F;

	double n_r = 0.5*(dm.n_r[B.shape()[1]-N_GC-2]+dm.n_r[B.shape()[1]-N_GC-1]); //sqrt(g_{rr})=exp(lambda/2) at r=Ro

    struct muConstraint_params mu_params = {Ri,Ro,n_r};

	F.function = &muConstraint;
	F.params = &mu_params;

	//Fit based on computing roots for set of values of Ri/Ro
//	double mu_tilde_fit_guess = 17.19799773*exp( 0.24030211*(pow( Ri/Ro/0.9, 23.35487559 )-1.) );
    double mu_tilde_fit_guess;
//    if( n_r < 1.0001 ) mu_tilde_fit_guess = 17.19799772*exp( 0.24030211*(pow( Ri/Ro/0.9, 3.35487553 )-1.) );
//    else mu_tilde_fit_guess = 17.64841792*exp( 0.23199471*(pow( Ri/Ro/0.9, 23.57917689 )-1.) );
//    if( n_r < 1.0001 ) mu_tilde_fit_guess = 17.19799772*exp( 0.24030211*(pow( Ri/Ro/0.9, 3.35487553 )-1.) );
    mu_tilde_fit_guess = 17.62536936*exp( 0.23283267*(pow( Ri/Ro/0.9, 23.55347586 )-1.) );
	double n = ceil( mu_tilde_fit_guess/pi - 0.5 ); //search for mu_tilde between (n-0.5)*pi and (n+0.5)*pi for integer n
	double mutilde_lo = (n-0.5)*pi + 1e-8, mutilde_hi = (n+0.5)*pi - 1e-8;
	//if mu_tilde_fit_guess is poor, adjust mutilde_lo and mutilde_hi appropriately
	if( muConstraintf_r(mutilde_lo,Ri/Ro,n_r)*muConstraintf_r(mutilde_hi,Ri/Ro,n_r) > 0. ){
        mutilde_lo = (n+0.5)*pi + 1e-8;
        mutilde_hi = (n+1.5)*pi - 1e-8;
    }

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, mutilde_lo, mutilde_hi);

	double curr_r = 0.;

	//computes mu_tilde. If unsuccessful, code will exit.
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		curr_r = gsl_root_fsolver_root (s);
		mutilde_lo = gsl_root_fsolver_x_lower (s);
		mutilde_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (mutilde_lo, mutilde_hi, 1e-10,1e-9);

		 if (status == GSL_SUCCESS)
			mu_tilde = curr_r;
	}
	while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free (s);

    double B_pol = Ro*Ro*B_pol_pole/(2.*f_r(mu_tilde,Ro,n_r,Ro))*sqrt(4.*pi/3.); //set B_pol such that B_r = B_pol_pole at r=Ro and theta=0. Choice of sign here means theta=0 is the north magnetic pole

    // Integrate analytic initial magnetic field configuration to determine cell face averages of each component to initialize finite volume solver

    struct BConfig_params params = {B_pol, B_tor_max, Ri, Ro, mu_tilde,n_r};

    size_t subints = 2000; //maximum number of integration subintervals
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (subints);
    int key = 3; //integration rule key
    double epsabs = 1e-7;

    double result1, result2, error;

    gsl_function Fr_theta;
    Fr_theta.function = &InitialBr_theta;
    Fr_theta.params = &params;

    gsl_function Ftheta_r;
    Ftheta_r.function = &InitialBtheta_r;
    Ftheta_r.params = &params;

    gsl_function Fphi_r;
    Fphi_r.function = &InitialBphi_r;
    Fphi_r.params = &params;

    gsl_function Fphi_theta;
    Fphi_theta.function = &InitialBphi_theta;
    Fphi_theta.params = &params;

    //Compute cell face-averaged value of Br, Btheta, Bphi for every cell except for the top, given the functions InitialBr, InitialBtheta, InitialBphi, which computes the analytic initial profile for Br, Btheta, Bphi
    //
    // B[0][i][j] = Deltaphi*int_{theta_j-Deltatheta_j/2}^{theta_j+Deltatheta_j/2}dtheta*sin(theta)*(r_i - Deltar_i/2)^2*B_r(r_i-Deltar_i/2,theta)/C_ij_r
    // where C_ij_r = Deltaphi*(r_i - Deltar_i/2)^2*int_{theta_j-Deltatheta_j/2}^{theta_j+Deltatheta_j/2}dtheta*sin(theta) = Deltaphi*(r_i - Deltar_i/2)^2*(2*sin(theta_j)*sin(Deltatheta_j/2))
    //
    // B[1][i][j] = Deltaphi*sin(theta_j-Deltatheta_j/2)*int_{r_i-Deltar_i/2}^{r_i+Deltar_i/2}dr*exp(lambda(r)/2)*r*B_theta(r,theta_j-Deltatheta_j/2)/C_ij_theta
    // where C_ij_theta = Deltaphi*sin(theta_j-Deltatheta_j/2)*int_{r_i-Deltar_i/2}^{r_i+Deltar_i/2}dr*exp(lambda(r)/2)*r = Deltaphi*sin(theta_j-Deltatheta_j/2)*r_i*Deltar_i*exp(lambda(r_i)/2)
    //
    // B[2][i][j] = int_{theta_j-Deltatheta_j/2}^{theta_j+Deltatheta_j/2}dtheta*int_{r_i-Deltar_i/2}^{r_i+Deltar_i/2}dr*exp(lambda(r)/2)*r*B_phi(r,theta_j-Deltatheta_j/2)/C_ij_phi
    // where C_ij_phi = int_{theta_j-Deltatheta_j/2}^{theta_j+Deltatheta_j/2}dtheta*int_{r_i-Deltar_i/2}^{r_i+Deltar_i/2}dr*exp(lambda(r)/2)*r = r_i*Deltar_i*Deltatheta*exp(lambda(r_i)/2)
    for(size_t i=N_GC; i<r.size()-N_GC; i++){
        for(size_t j=N_GC; j<theta.size()-N_GC; j++){
//            gsl_integration_qag(&Fr_theta, theta[j]-Deltatheta/2., theta[j]+Deltatheta/2., 0, epsabs, subints, key, w, &result2, &error);
            if(i < r.size()-N_GC-1){
                result2 = InitialBr_theta_analytic(theta[j]-Deltatheta/2.,theta[j]+Deltatheta/2.);
                B[0][i][j] = result2*InitialBr_r( r[i]-Deltar[i]/2., &params )/( 2.*sin(theta[j])*sin(Deltatheta/2.) );
                gsl_integration_qag(&Ftheta_r, r[i]-Deltar[i]/2., r[i]+Deltar[i]/2., 0, epsabs, subints, key, w, &result1, &error);
                B[1][i][j] = result1*InitialBtheta_theta( theta[j]-Deltatheta/2., &params )/(dm.n_r[i]*Deltar[i]*dm.r[i]);
                result1 = InitialBphi_r_analytic(r[i]-Deltar[i]/2., r[i]+Deltar[i]/2., &params);
                result2 = InitialBphi_theta_analytic(theta[j]-Deltatheta/2.,theta[j]+Deltatheta/2.);
                // gsl_integration_qag(&Fphi_r, r[i]-Deltar[i]/2., r[i]+Deltar[i]/2., 0, epsabs, subints, key, w, &result1, &error);
                // gsl_integration_qag(&Fphi_theta, theta[j]-Deltatheta/2., theta[j]+Deltatheta/2., 0, epsabs, subints, key, w, &result2, &error);
                B[2][i][j] = result1*result2/(Deltar[i]*Deltatheta*r[i]);
            }
            else{
                result2 = InitialBr_theta_analytic(theta[j]-Deltatheta/2.,theta[j]+Deltatheta/2.);
                B[0][i][j] = result2*InitialBr_r( r[i-1]+Deltar[i]/2., &params )/( 2.*sin(theta[j])*sin(Deltatheta/2.) );
                gsl_integration_qag(&Ftheta_r, r[i], r[i]+Deltar[i], 0, epsabs, subints, key, w, &result1, &error);
                B[1][i][j] = 2.*result1*InitialBtheta_theta( theta[j]-Deltatheta/2., &params )/(0.5*(dm.n_r[i]+dm.n_r[i-1])*Deltar[i]*(dm.r[i]+0.5*Deltar[i])) - B[1][i-1][j]; //initialize Btheta continuous across the outer boundary. This gets changed by B_BoundaryConditions, but is needed to get the Btheta BC correct initially.
                B[2][i][j] = -B[2][i-1][j]; //initialize Bphi to be zero at the boundary. This gets changed by B_BoundaryConditions.
            }
        }
    }

    gsl_integration_workspace_free (w);

    return;
}

/*
        Function f(r) appearing in initial poloidal field (Aguilera D. N., Pons J. A., Miralles J. A., 2008, A&A, 486, 255)
        Inputs: mutilde: mu*Ro. Determined by imposing outer boundary condition
                Ro: outer radius in reduced units
                n_r: sqrt(g_rr) = exp(lambda/2) evaluated at r=Ro
                r: radial coordinate in reduced units
        Output: f: value of function f at r
*/
double f_r(double mutilde, double Ro, double n_r, double r)
{
    double f_1 = ( 4.*n_r*n_r/( 4. - 3.*( n_r - 1. ) ) + 2. )/3.; //for l = 1
    double fscr_1 = 3.*mutilde*(1.-f_1)/( 3.*( 1. - f_1  )-mutilde*mutilde );
    double b = ( tan(mutilde) - fscr_1 )/( fscr_1*tan(mutilde) + 1. );

    return mutilde*r/Ro*( sin(mutilde*r/Ro)/pow(mutilde*r/Ro,2.) - cos(mutilde*r/Ro)/(mutilde*r/Ro) + ( -cos(mutilde*r/Ro)/pow(mutilde*r/Ro,2.) - sin(mutilde*r/Ro)/(mutilde*r/Ro) )*b );
}

/*
        Function df(r)/dr appearing in initial poloidal field (Aguilera D. N., Pons J. A., Miralles J. A., 2008, A&A, 486, 255)
        Inputs: mutilde: mu*Ro. Determined by imposing outer boundary condition
                Ro: outer radius in reduced units
                n_r: sqrt(g_rr) = exp(lambda/2) evaluated at r=Ro
                r: radial coordinate in reduced units
        Output: df/dr: value of function df/dr at r
*/
double df_rdr(double mutilde, double Ro, double n_r, double r)
{
    double f_1 = ( 4.*n_r*n_r/( 4. - 3.*( n_r - 1. ) ) + 2. )/3.; //for l = 1
    double fscr_1 = 3.*mutilde*(1.-f_1)/( 3.*( 1. - f_1  )-mutilde*mutilde );
    double b = ( tan(mutilde) - fscr_1 )/( fscr_1*tan(mutilde) + 1. );

    return f_r(mutilde,Ro,n_r,r)/r + mutilde/Ro*( 2.*( cos(mutilde*r/Ro) + sin(mutilde*r/Ro)*b )/(mutilde*r/Ro)
                                                            + ( 1. - 2./pow(mutilde*r/Ro,2.) )*( sin(mutilde*r/Ro) - cos(mutilde*r/Ro)*b ) );
}

/*
        Constraint equation for poloidal field parameter mutilde
        Input: mutilde: dimensionless parameter in function f.
               params: muConstraint_params object containing Ri and Ro, the inner and outer radii in reduced units
        Output: value of constraint function
*/
double muConstraint(double mutilde, void *params)
{
	struct muConstraint_params *p = (struct muConstraint_params *) params;
	double Ri = p->Ri;
	double Ro = p->Ro;
	double n_r = p->n_r;

    return muConstraintf_r(mutilde,Ri/Ro,n_r);
}

/*
        Simplified version of the function f(r=Ri) used to compute mutilde given Ri.
        Inputs: mutilde: mu*Ro. Determined by imposing outer boundary condition
                R_ratio: inner radius divided by outer radius Ri/Ro < 1
                n_r: sqrt(g_rr) = exp(lambda/2) evaluated at r=Ro
        Output: funceval: value of function f at r=Ri. Should be zero.
*/
double muConstraintf_r(double mutilde, double R_ratio, double n_r)
{
    double funceval;

    if(n_r > 1.+1e-9){ //if considering general relativity
        double f_1 = ( 4.*n_r*n_r/( 4. - 3.*( n_r - 1. ) ) + 2. )/3.; //for l = 1
        double fscr_1 = 3.*mutilde*(1.-f_1)/( 3.*( 1. - f_1  )-mutilde*mutilde );
        funceval = ( fscr_1*tan(mutilde) + 1. )*( sin(mutilde*R_ratio)/(mutilde*R_ratio) - cos(mutilde*R_ratio) ) - ( tan(mutilde) - fscr_1 )*( cos(mutilde*R_ratio)/(mutilde*R_ratio) + sin(mutilde*R_ratio) );
    }
    else{ //if not considering general relativity
        funceval = ( sin(mutilde*R_ratio)/(mutilde*R_ratio) - cos(mutilde*R_ratio) ) - tan(mutilde)*( cos(mutilde*R_ratio)/(mutilde*R_ratio) + sin(mutilde*R_ratio) );
    }

    return funceval;
}

/*
        Initializes r-component of B field: r-dependent function
        Arguments: r: radial coordinate
                   params: parameters for function (needed by GSL even if empty)
        Output: function value B_r(r)
*/
double InitialBr_r(double r, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_pol = p -> B_pol;
    double mu = p -> mu_tilde;
    double Ro = p -> Ro;
    double n_r = p -> n_r;

    //Do not include extra factor of r arising from Jacobian of line integral over an arc with fixed r, since B_r is averaged over an arc of angular size Deltatheta at fixed r
    //so dividing by r*Deltatheta removes this extra factor
    return 2.*B_pol/(r*r)*f_r(mu,Ro,n_r,r);
}

/*
        Initializes r-component of B field times sin(theta): theta-dependent function
        Arguments: theta: polar angle coordinate
                   params: parameters for function (needed by GSL even if empty)
        Output: function value B_r(theta)
*/
double InitialBr_theta(double theta, void * params)
{
    return sqrt(3./(4.*pi))*cos(theta)*sin(theta);
}

/*
        Initializes r-component of B field times sin(theta): theta-dependent function
        Analytic integral if simple enough to compute
        Arguments: theta2, theta1: polar angle coordinates at upper and lower range of integral
        Output: function value B_r(theta)
*/
double InitialBr_theta_analytic(double theta1, double theta2)
{
    return sqrt(3./(4.*pi))*(sin(theta2)*sin(theta2)/2.-sin(theta1)*sin(theta1)/2.);
}

/*
        Initializes theta-component of B field: r-dependent function times Jacobian factor r*e^{\lambda/2}
        Arguments: r: radial coordinate
                   params: parameters for function (needed by GSL even if empty)
        Output: function value B_theta(r)
*/
double InitialBtheta_r(double r, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_pol = p -> B_pol;
    double mu = p -> mu_tilde;
    double Ro = p -> Ro;
    double n_r = p -> n_r;

    //Averaging over line integral along r at fixed theta, so no aditional r factor in Jacobian for line integral
    //The n_r*r = r*e^{\lambda/2} factor in the Jacobian cancels the factor e^{-\lambda/2}/r in the function Btheta
    return B_pol*df_rdr(mu,Ro,n_r,r);
}

/*
        Initializes theta-component of B field: theta-dependent function
        Arguments: theta: polar angle coordinate
                   params: parameters for function (needed by GSL even if empty)
        Output: function value B_theta(theta)
*/
double InitialBtheta_theta(double theta, void * params)
{
    return -sqrt(3./(4.*pi))*sin(theta);
}

/*
        Initializes phi-component of B field
        Arguments: r: radial coordinate
                  params: parameters for function (needed by GSL even if empty)
        Output: function value B_phi(r)
*/
double InitialBphi_r(double r, void * params)
{
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_tor_max = p -> B_tor_max;
    double Ri = p -> Ri;
    double Ro = p -> Ro;

    //Include extra factor of r arising from Jacobian of areal integral
    return r*32.*B_tor_max*sqrt(4.*pi/45.)*pow(r-Ri,2.)*pow(Ro-r,2.)/pow( Ro-Ri,4. );
}

/*
        Initializes phi-component of B field
        Arguments: theta: polar angle coordinate
                  params: parameters for function (needed by GSL even if empty)
        Output: function value B_phi(theta)
*/
double InitialBphi_theta(double theta, void * params)
{
    return -sqrt(45./(4.*pi))*sin(theta)*cos(theta);
}

/*
        Initializes phi-component of B field times sin(theta): r-dependent function
        Analytic integral if simple enough to compute
        Arguments: theta2, theta1: polar angle coordinates at upper and lower range of integral
        Output: function value B_r(theta)
*/
double InitialBphi_r_analytic(double r1, double r2, void * params) {
    struct BConfig_params *p = (struct BConfig_params *) params;
    double B_tor_max = p -> B_tor_max;
    double Ri = p -> Ri;
    double Ro = p -> Ro;

    return 32.*B_tor_max*sqrt(4.*pi/45.)/pow( Ro-Ri,4. )*( r2*r2/60.*( 10.*r2*r2*r2*r2 + 30.*Ri*Ri*Ro*Ro - 24.*r2*r2*r2*(Ri+Ro) -40.*r2*Ro*Ri*(Ri+Ro) + 15.*r2*r2*(Ri*Ri+4.*Ri*Ro+Ro*Ro) )
                                                         - r1*r1/60.*( 10.*r1*r1*r1*r1 + 30.*Ri*Ri*Ro*Ro - 24.*r1*r1*r1*(Ri+Ro) -40.*r1*Ro*Ri*(Ri+Ro) + 15.*r1*r1*(Ri*Ri+4.*Ri*Ro+Ro*Ro) ) );
}

/*
        Initializes phi-component of B field times sin(theta): theta-dependent function
        Analytic integral if simple enough to compute
        Arguments: theta2, theta1: polar angle coordinates at upper and lower range of integral
        Output: function value B_r(theta)
*/
double InitialBphi_theta_analytic(double theta1, double theta2) {

    return -sqrt(45./(4.*pi))*( -cos(2.*theta2) + cos(2.*theta1) )/4.;
}

/*
    Initializes the cell face-averaged temperature. Initializes ghost cells as well.
    Inputs: r, theta: vectors of the cell mid-points in each direction
    tparams: TParams class object containing properties of the initial temperature profile
    domain: Domain class object containing the size of the full domain as Lr and Ltheta
    N_GC: number of ghost cells added in each direction
    Deltar, Deltatheta: grid spacing in each direction
    Output: T: 1D multi_array whose first dimension corresponds to the three field components and the final two dimensions are the grid points
*/
void InitializeT(std::vector<double> & r, std::vector<double> & theta, const TParams& tparams, const Domain & domain, size_t N_GC, std::vector<double> & Deltar, double Deltatheta, ScalarField & T)
{

    double Ri = domain.r_min;
    double Ro = domain.r_max;

    double T_init = tparams.T_init;
    bool uniform_T = tparams.uniform_T;

    struct TConfig_params params = {T_init, Ri, Ro};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error;

    gsl_function T_r;
    T_r.function = &InitialT_r;
    T_r.params = &params;

    //Compute cell face-averaged value of Br, Btheta, Bphi for every cell except for the top, given the functions InitialBr, InitialBtheta, InitialBphi, which computes the analytic initial profile for Br, Btheta, Bphi
    for(size_t i=0; i<r.size()-1; i++){
        for(size_t j=0; j<theta.size(); j++){
            if(uniform_T == true){
                T[i][j] = T_init;
            }
            else{
                gsl_integration_qag(&T_r, r[i]-Deltar[i]/2., r[i]+Deltar[i]/2., 0, 1e-7, 1000, 3, w, &result, &error);
                T[i][j] = result/Deltar[i];
            }
        }
    }
    //Initializes outermost temperatures to zero (these are unused)
    for(size_t j=0; j<theta.size(); j++){
        T[r.size()-1][j] = 0.;
    }
    //Initialize the ghost cells at the r=rmin boundary equal to T_init. Needed for correct implementation of boundary conditions.
    for(size_t i=0; i<N_GC; i++){
        for(size_t j=N_GC; j<theta.size()-N_GC; j++){
            T[i][j] = T_init;
        }
    }

    gsl_integration_workspace_free (w);

    return;
}

/*
        Initializes temperature
        Arguments: r: radial coordinate
                params: parameters for function (needed by GSL even if empty)
        Output: function value B_phi(theta)
*/
double InitialT_r(double r, void * params)
{
    struct TConfig_params *p = (struct TConfig_params *) params;
    double T_init = p -> T_init;
    double Ri = p -> Ri;
    double Ro = p -> Ro;

    double rho_i = 1.22e14, rho_o = 1.41e10;
    double slope = 5.05;
    double constant = -10.;
    double power = 2.;

//    double rho_i = 1.22e14, rho_o = 1.41e10;
//    double slope = 9.;
//    double constant = 10.;
//    double power = 1.;


    return T_init*log( rho_i/rho_o*exp( -slope*( pow( ( Ri - r )/( Ri - Ro ),power ) ) ) + constant )/log( rho_i/rho_o + constant );
}
