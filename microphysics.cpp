#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_expint.h>

#include "common.h"
#include "microphysics.h"

/*
    Loads crust equation of state file and generates interpolating functions for crust microphysics properties
    Inputs: filename: name of EOS file, including file extension
            rho_cutoff: low density cutoff for EOS in g/cm^3. xmax will thus correspond to the radius at which rho=rho_cutoff
    Output: Interpolators: object of EOSInterpolation class to contain EOS interpolators and required auxiliary objects
            g_s14: surface gravitational acceleration in units of 1e14 cm/s^2
            r_s: outer radius of star (not simulation domain) in reduced units
            n_t_s: lapse function n_t = sqrt(-g_{tt}) = exp(nu/2) at outer radius of star
            n_b_nd: neutron drip baryon density in fm^-3
*/
void load_EOS(std::string CrustEOS_filename, double rho_cutoff, double & g_s14, double & r_s, double & n_t_s, double & n_b_nd, EOSInterpolation & Interpolators)
{
    double km_to_cm = 1e5; //converts radius in km to cm

	//Determines number of rows in data table
	std::ostringstream filename;
	filename << CrustEOS_filename;
	std::ifstream Datafile;
    Datafile.open(filename.str().c_str(),std::ifstream::in);
	std::string line, dummyLine;
    size_t N_dat = 0; //number of rows in data table. Initialize to zero
	getline(Datafile, dummyLine); //skips header (1 line)
	while (std::getline(Datafile, line))
		++N_dat;

	Datafile.close();

    //vectors to store values from data table
	std::vector<double> r_InterpVec; //radius in reduced units. Vector; will convert to array
	std::vector<double> gtt_InterpVec; //Metric factor g_{tt}=-e^{nu}. Vector; will convert to array
	std::vector<double> grr_InterpVec; //Metric factor g_{rr}=e^{lambda}. Vector; will convert to array
	std::vector<double> n_b_InterpVec; //Baryon number density in fm^{-3}. Vector; will convert to array
    std::vector<double> rho_InterpVec; //mass density in g/cm^3. Vector; will convert to array. Note: not equal to the mass-energy density called rho in the data table!
	std::vector<double> n_i_InterpVec; //ion number density in fm^{-3}. Vector; will convert to array
	std::vector<double> A_InterpVec; //Mass number of crustal nuclei (excluding dripped neutrons and protons). Vector; will convert to array
	std::vector<double> Z_InterpVec; //Atomic number of crustal nuclei, including dripped protons (called Z_eq in data table). Vector; will convert to array
	std::vector<double> n_e_InterpVec; //Electron number density in fm^{-3}. Vector; will convert to array
	std::vector<double> n_nf_InterpVec; //Number density of free (dripped) neutrons in fm^{-3}. Vector; will convert to array
    std::vector<double> mu_nf_InterpVec; //Free (dripped) neutron chemical potential in MeV. Vector; will convert to array

    //temporary variables to hold data table values
    double r_dt, M_dt, rho_dt, gtt_dt, grr_dt, n_b_dt, n_e_dt, n_nf_dt, mu_nf_dt;
    double Aeq_dt, A_dt, Zeq_dt, Z_dt; //Aeq = baryons per nucleus > A = mass number, Zeq = protons per nucleus > Z = atomic number. A and Z are per nucleus (outer crust) or "cluster" (inner crust)
    double dummy; //dummy variable to avoid saving all quantities from loaded data table

    int NDTrigger = 0; //trigger for determining neutron drip density

	Datafile.open(filename.str().c_str(),std::fstream::in);
	std::getline(Datafile, dummyLine); //skips header (1 line)
	//Gets background TOV solution values from data table
	for(size_t i=0; i<N_dat; i++)
	{
		Datafile >> r_dt >> dummy >> M_dt >> rho_dt >> gtt_dt >> grr_dt >> n_b_dt >> Aeq_dt >> A_dt >> Zeq_dt >> Z_dt >> n_e_dt >> n_nf_dt >> mu_nf_dt;

        //Determine neutron drip density
		if(n_nf_dt < 1e-30 && NDTrigger == 0){
            n_b_nd = ( n_b_dt + n_b_InterpVec.back() )/2.;
            NDTrigger = 1;
		}

		r_InterpVec.push_back(r_dt*km_to_cm/L_0);
		gtt_InterpVec.push_back(gtt_dt);
		grr_InterpVec.push_back(grr_dt);
		n_b_InterpVec.push_back(n_b_dt);
		rho_InterpVec.push_back(M_n*MeVtoErg*pow(1e13,3)/pow(c,2.)*n_b_dt);
		n_i_InterpVec.push_back(n_e_dt/Zeq_dt);
		//n_i_InterpVec.push_back(rho_dt/(Ap_dt*1.66053906892e-24)/1e39);
		A_InterpVec.push_back(A_dt);
		Z_InterpVec.push_back(Z_dt);
		n_e_InterpVec.push_back(n_e_dt);
        //n_e_InterpVec.push_back(rho_dt/(Ap_dt*1.66053906892e-24)/1e39*Z_dt);
		n_nf_InterpVec.push_back(n_nf_dt);
		mu_nf_InterpVec.push_back(mu_nf_dt);

	}
	Datafile.close();

	g_s14 = G*M_dt*M_solar/pow(r_dt*km_to_cm,2.)*sqrt(grr_dt)/1e14; //gravitational acceleration at the surface of the star in units of 1e14 cm/s^2
    r_s = r_dt*km_to_cm/L_0; //outer (surface) radius of star in reduced units
    n_t_s = sqrt(-gtt_dt); //lapse function sqrt(-g_{tt}) = exp(nu/2) at outer radius of star

	//converts vectors of data to arrays to use in interpolation
	double* r_Interp = &r_InterpVec[0]; //radius in reduced units
	double* gtt_Interp = &gtt_InterpVec[0]; //Metric factor g_{tt}=-e^{nu}
	double* grr_Interp = &grr_InterpVec[0]; //Metric factor g_{rr}=e^{lambda}
	double* n_b_Interp = &n_b_InterpVec[0]; //Baryon number density in fm^{-3}
    double* rho_Interp = &rho_InterpVec[0]; //mass density in g/cm^3
    double* n_i_Interp = &n_i_InterpVec[0]; //ion number density in fm^{-3}
	double*	A_Interp = &A_InterpVec[0]; //Mass number of crustal nuclei
	double* Z_Interp = &Z_InterpVec[0]; //Atomic number of crustal nuclei
	double* n_e_Interp = &n_e_InterpVec[0]; //Electron number density in fm^{-3}
	double* n_nf_Interp= &n_nf_InterpVec[0]; //Number density of free (dripped) neutrons in fm^{-3}
    double* mu_nf_Interp= &mu_nf_InterpVec[0]; //Chemical potential of free (dripped) neutrons in MeV

	//Generates spline and acc objects for each quantity to interpolate them and assigns them to EOSInterpolation object
	gsl_interp_accel *acc_gtt = gsl_interp_accel_alloc ();
	gsl_spline *spline_gtt = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_gtt, r_Interp, gtt_Interp, N_dat);
	Interpolators.gtt_spline = spline_gtt;
	Interpolators.gtt_acc = acc_gtt;

	gsl_interp_accel *acc_grr = gsl_interp_accel_alloc ();
	gsl_spline *spline_grr = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_grr, r_Interp, grr_Interp, N_dat);
	Interpolators.grr_spline = spline_grr;
	Interpolators.grr_acc = acc_grr;

	gsl_interp_accel *acc_n_b = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_b = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_b, r_Interp, n_b_Interp, N_dat);
	Interpolators.n_b_spline = spline_n_b;
	Interpolators.n_b_acc = acc_n_b;

    gsl_interp_accel *acc_rho = gsl_interp_accel_alloc ();
	gsl_spline *spline_rho = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_rho, r_Interp, rho_Interp, N_dat);
	Interpolators.rho_spline = spline_rho;
	Interpolators.rho_acc = acc_rho;

    gsl_interp_accel *acc_n_i = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_i = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_i, r_Interp, n_i_Interp, N_dat);
	Interpolators.n_i_spline = spline_n_i;
	Interpolators.n_i_acc = acc_n_i;

    gsl_interp_accel *acc_A = gsl_interp_accel_alloc ();
	gsl_spline *spline_A = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_A, r_Interp, A_Interp, N_dat);
	Interpolators.A_spline = spline_A;
	Interpolators.A_acc = acc_A;

    gsl_interp_accel *acc_Z = gsl_interp_accel_alloc ();
	gsl_spline *spline_Z = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_Z, r_Interp, Z_Interp, N_dat);
	Interpolators.Z_spline = spline_Z;
	Interpolators.Z_acc = acc_Z;

    gsl_interp_accel *acc_n_e = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_e = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_e, r_Interp, n_e_Interp, N_dat);
	Interpolators.n_e_spline = spline_n_e;
	Interpolators.n_e_acc = acc_n_e;

    gsl_interp_accel *acc_n_nf = gsl_interp_accel_alloc ();
	gsl_spline *spline_n_nf = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_n_nf, r_Interp, n_nf_Interp, N_dat);
	Interpolators.n_nf_spline = spline_n_nf;
	Interpolators.n_nf_acc = acc_n_nf;

    gsl_interp_accel *acc_mu_nf = gsl_interp_accel_alloc ();
	gsl_spline *spline_mu_nf = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_mu_nf, r_Interp, mu_nf_Interp, N_dat);
	Interpolators.mu_nf_spline = spline_mu_nf;
	Interpolators.mu_nf_acc = acc_mu_nf;

    Interpolators.R_cc = r_InterpVec[0]; //coordinate radius of the crust-core transition in reduced units
	Interpolators.RStar = r_InterpVec.back(); //coordinate radius of the star in reduced units

	std::reverse(rho_InterpVec.begin(), rho_InterpVec.end()); //flip order of energy density vector (need to do this since rho decreases with r)
	std::reverse(r_InterpVec.begin(), r_InterpVec.end()); //flip order of radial coordinate vector (need to do this since rho decreases with r)
	double* r_Interp_rev = &r_InterpVec[0]; //radius in cm
	double* rho_Interp_rev = &rho_InterpVec[0]; //energy density in g/cm^3

    gsl_interp_accel *acc_r = gsl_interp_accel_alloc ();
	gsl_spline *spline_r = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(spline_r, rho_Interp_rev, r_Interp_rev, N_dat);

	Interpolators.R_rhocutoff = gsl_spline_eval(spline_r, rho_cutoff, acc_r); //coordinate radius at cutoff density rho_cutoff

    gsl_spline_free (spline_r);
    gsl_interp_accel_free (acc_r);

	return;
}

/*
        Initializes general relativity factors using interpolating functions generated from
        Inputs: GR: true if general relativity is turned on and false otherwise
                EOS_Interps: object of EOSInterpolation class containing EOS interpolators and required auxiliary objects
                r: radial coordinates of cell centers in reduced units
                Nr: number of non-ghost cells in radial direction
                N_GC: number of ghost cells
                Output: n_t, n_r: lapse function sqrt(-g_{tt}) and sqrt(g_{rr})
*/
void GR_Factor_Initialize(bool GR, EOSInterpolation & EOS_Interps, std::vector<double> & r, size_t Nr, size_t N_GC, std::vector<double> & n_t, std::vector<double> & n_r)
{
    if( GR == true ){
        for(size_t i=0; i<Nr+2*N_GC; i++){
            if( i < N_GC ){
                //Mirrored values across inner boundary. Boundary is at mid-point between cells i=N_GC-1 and N_GC
                n_t.push_back( sqrt( -gsl_spline_eval(EOS_Interps.gtt_spline, r[2*N_GC-1-i], EOS_Interps.gtt_acc) ) );
                n_r.push_back( sqrt( gsl_spline_eval(EOS_Interps.grr_spline, r[2*N_GC-1-i], EOS_Interps.grr_acc) ) );
            }
            else if( i > Nr+N_GC-2 ){
                //Mirrored values across outer boundary. Boundary is at mid-point between cells i=Nr+N_GC-1 and Nr+N_GC-2
                n_t.push_back( sqrt( -gsl_spline_eval(EOS_Interps.gtt_spline, r[2*Nr-i+1], EOS_Interps.gtt_acc) ) );
                n_r.push_back( sqrt( gsl_spline_eval(EOS_Interps.grr_spline, r[2*Nr-i+1], EOS_Interps.grr_acc) ) );
            }
            else{
                n_t.push_back( sqrt( -gsl_spline_eval(EOS_Interps.gtt_spline, r[i], EOS_Interps.gtt_acc) ) );
                n_r.push_back( sqrt( gsl_spline_eval(EOS_Interps.grr_spline, r[i], EOS_Interps.grr_acc) ) );
            }
        }
        //Since i=Nr+N_GC-1 cell is outside of simulation domain, set its metric components equal to those of i=Nr+N_GC-2 cell
        n_t[Nr+N_GC-1] = n_t[Nr+N_GC-2];
        n_r[Nr+N_GC-1] = n_r[Nr+N_GC-2];
    }
    else{
        for(size_t i=0; i<Nr+2*N_GC; i++){
            n_t.push_back( 1. );
            n_r.push_back( 1. );
        }
    }

    return;
}

/*
    Sets superfluid gap parameters based on simulation parameters.
    Gap models from Table II. of Ho, Elshamouty, Heinke and Potekhin, PRC 91, 015806 (2015)
    Delta0 are in MeV, k0 and k2 are in fm^-1, k1 and k3 are in fm^-2
    Inputs: tparams: TParams object containing information about nucleon superfluid gaps.
*/
void SFgaps(TParams & tparams){

    // 1S0 neutron pairing in the crust
    if(tparams.SF_n_crust == "SFB"){
        tparams.Delta0_nCrust = 45.;
        tparams.k0_nCrust = 0.1;
        tparams.k1_nCrust = 4.5;
        tparams.k2_nCrust = 1.55;
        tparams.k3_nCrust = 2.5;
    }
    else if(tparams.SF_n_crust == "Ho2012"){
        tparams.Delta0_nCrust = 68.;
        tparams.k0_nCrust = 0.1;
        tparams.k1_nCrust = 4.;
        tparams.k2_nCrust = 1.7;
        tparams.k3_nCrust = 4.;
    }
    else
        std::cout << "Error: unimplemented crust neutron superfluid gap" << std::endl;

    // 1S0 proton pairing in the core
    if(tparams.SF_p_core == "CCDK"){
        tparams.Delta0_pCore = 102.;
        tparams.k0_pCore = 0.;
        tparams.k1_pCore = 9.;
        tparams.k2_pCore = 1.3;
        tparams.k3_pCore = 1.5;
    }
    else if(tparams.SF_p_core == "Ho2012"){
        tparams.Delta0_pCore = 120.;
        tparams.k0_pCore = 0.;
        tparams.k1_pCore = 9.;
        tparams.k2_pCore = 1.3;
        tparams.k3_pCore = 1.8;
    }
    else
        std::cout << "Error: unimplemented core proton superfluid gap" << std::endl;

    // 3P2, mJ=0 neutron ("type B") pairing in the core
    if(tparams.SF_n_core == "TToa"){
        tparams.Delta0_nCore = 2.1;
        tparams.k0_nCore = 1.1;
        tparams.k1_nCore = 0.6;
        tparams.k2_nCore = 3.2;
        tparams.k3_nCore = 2.4;
    }
    else if(tparams.SF_n_core == "SYHHP"){
        tparams.Delta0_nCore = 1.;
        tparams.k0_nCore = 2.08;
        tparams.k1_nCore = 0.04;
        tparams.k2_nCore = 2.7;
        tparams.k3_nCore = 0.013;
    }
    else if(tparams.SF_n_core == "EEHOr"){
        tparams.Delta0_nCore = 0.23;
        tparams.k0_nCore = 1.2;
        tparams.k1_nCore = 0.026;
        tparams.k2_nCore = 1.6;
        tparams.k3_nCore = 0.0080;
    }
    else
        std::cout << "Error: unimplemented core neutron superfluid gap" << std::endl;

    return;
}

/*
        Crust electrical and thermal conductivity in s^{-1}. Returns Ohmic diffusivity in reduced units and/or thermal conductivity in reduced units.
        Includes no Landau quantization effects

    For electron-ion scattering (thermal and electrical conductivity), uses results from Potekhin A. Y., Baiko D. A., Haensel P. and Yakovlev D. G., 1999, A&A 346, 345,
    Potekhin A. Y., 1999, A&A 351, 787, Gnedin, O. Y. et al, 2001, MNRAS 324, 725
    and reaction rates and transport in neutron stars by A. Schmitt and P. Shternin
    Use Pearson et al. MNRAS 481, 2994 (2018), Eq. (C21) and Table C10 for computing x_nuc=x_p for BSk24, fit to Douchin and Haensel, A&A 380, 151 (2001) data for SLy4.
    Note: Eq. (C21) has a typo in the denominator of the first term: 1 + p_4... should be 1 - p_4... (checked via their EOS code). Also replaced p_6 with 0.325 to improve fit
    For electron-electron scattering contribution (thermal conductivity) uses Shternin and Yakovlev PRD 74, 043004 (2006) with non-degenerate electron correction
    from Cassisi et al. ApJ, 661, 1094 (2007)
    For ion (phonon) contribution (thermal conductivity) uses Chugunov and Haensel MNRAS 381, 1143 (2007) with typo corrections based on Potekhin's conductivity code

    Inputs: T: temperature in reduced units
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
    Output: eta_O: ScalarField object of Ohmic diffusivity in reduced units
            kappa: ScalarField object of thermal conductivity in reduced units
            max_eta_T: maximum value of thermal diffusivity across simulation domain (reduced units)
*/
void ConductivityCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC, std::vector<double> & n_t,
                ScalarField & c_v, ScalarField & eta_O, ScalarField & kappa, double & max_eta_T)
{
    static double un1 = 2.78; //n = -1 frequency moment of bcc Coulombic lattice
    static double un2 = 12.973; //n = -2 frequency moment of bcc Coulombic lattice
    static double hbar = 6.582119569e-22; //reduced Planck constant in MeV*s
    static double n_drip = cparams.n_b_nd; //neutron drip density in fm^{-3} from equation of state
    static double kappa_pref = 1e39*pi*pi*k_Bcgs*c*c; //constant prefactor for computing thermal conductivities. 1e39 factor converts n_e to cm^{-3} and MeVtoErg/(c*c) converts mu_e to erg/c^2 = g

    double T_MeV, k_F, x_r, mu_e, a_i, q_BZ, Gamma, T_plasma, eta, eta_0, Z_cell;
    double v_F, beta, q_D, q_i2, k_TF2, w0, x_nuc, x_nuct, tp, q_s, s, w1, D;
    double G0, G2, G_s_sigma, G_s_kappa, w, sw, sw1, svF2, expi_sw, expi_sww, expi_sw1, expi_sw1w1, Lambda_1 = 0., Lambda_2 = 0., Lambda0_1 = 0., Lambda0_2 = 0.;
    double tau_ee, T_plasma_e, th, I_l, I_t, I_lt, A_fact, Bpow, C, C_1, C_2, T_TF, t, NonDegenerateFactor_ee; //for electron-electron scattering thermal conductivity
    double CoulombLogVal, tau_ei, tau_e;
    double kappa_i, kappa_ii, kappa_ie, LogFact, kappa_0, c_s, F_th, w_DW, w_form, y, exp_nw, exp_nwy2, Lambda_phe, Arho6Ap, L_ph; //for ionic thermal conductivity
    double sigma; //electrical conductivity in s^{-1}

    double Q = 0., s_imp, w1_imp, s_impvF2, Lambda_1_imp, Lambda_2_imp, Lambda0_1_imp, Lambda0_2_imp, CoulombLog_imp, tau_ei_imp;

    // parameters for specific heat capacity of a bcc lattice
    static double alpha[5] = {0.932446, 0.334547, 0.265764, 4.757014e-3, 4.7770935e-3};
    static double an[9] = {1., 0.1839, 0.593586, 5.4814e-3, 5.01813e-4, 0, 3.2947e-7, 0, 5.8356e-11};
    static double bn[8] = {261.66, 0, 7.07997, 0, 0.0409484, 3.97355e-4, 5.11148e-5, 2.19749e-6};
    double C_vi, th2, th3, th4, th5, th6, th7, th8, th9, th10, th11, Ath, Bth, Athp, Athpp, Bthp, Bthpp;

    for(size_t i=0; i<eta_O.shape()[0]-N_GC+1; i++){

        T_plasma = sqrt(756.40518016759142*Z[i]*Z[i]*n_i[i]/A[i]); //plasma temperature in MeV. 756.40518016759142 = 4*pi*alpha_e*(hbar*c)^3/(1 amu*c^2)
        eta_0 = 0.19/pow(Z[i],1./6.);
        a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm
        q_BZ = pow( 6.*pi*pi*n_i[i],1./3. ); //Brillouin zone boundary wavenumber in fm^{-1}

        //Set Q = Z_imp^2 in A. Potekhin's codes
        if(cparams.Z_impurity == "Carreau2020BSk24"){
            Q = Q_impCarreau2020BSk24(n_b[i],n_drip);
        }
        else if(cparams.Z_impurity == "Vigano2013Q100"){
            Q = Q_impVigano2013Q100(n_b[i]);
        }
        else{
            Q = 0.;
        }

        for(size_t j=0; j<eta_O.shape()[1]; j++){

            if(i < N_GC){
                T_MeV = T[N_GC+1-i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            else{
                T_MeV = T[i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            Gamma = alpha_e*hbarc*Z[i]*Z[i]/(T_MeV*a_i); //plasma coupling parameter (dimensionless)
            eta = T_MeV/T_plasma; //ratio of temperature to plasma temperature
            k_F = pow(3.*pi*pi*n_e[i][j],1./3.)*hbarc; //Fermi wavenumber = Fermi momentum in MeV
            x_r = k_F/M_e; //electron relativity parameter. Also used for p_0bar in some formulae
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
                Lambda_1_imp = 0.5*( log(1.+1./s_imp ) - 1./(1.+s_imp ) );
                Lambda_2_imp = v_F*v_F*( (2.*s_imp +1.)/(2.*s_imp +2.) - s_imp*log(1.+1./s_imp ) );
                Lambda0_1_imp = ( log(1.+1./s_imp) + 1./(1.+1./s_imp)*(1.-exp(-w1_imp)) - (1.+s_imp*w1_imp)*( expE1_rp(s_imp*w1_imp)-exp(-w1_imp)*expE1_rp(s_imp*w1_imp+w1_imp) ) )/2.;
                Lambda0_2_imp = ( v_F*v_F*(exp(-w1_imp)-1.+w1_imp)/w1_imp - s_impvF2/(1.+1./s_imp)*(1.-exp(-w1_imp)) - 2.*s_impvF2*log(1.+1./s_imp) + s_impvF2*(2.+s_imp*w1_imp)*( expE1_rp(s_imp*w1_imp)-exp(-w1_imp)*expE1_rp(s_imp*w1_imp+w1_imp) ) )/2.;
                CoulombLog_imp = (Lambda_1_imp - Lambda_2_imp) - (Lambda0_1_imp - Lambda0_2_imp); //Coulomb logarithm for impurity scattering
                tau_ei_imp = 5.699361922e-17/( Q/Z[i]*CoulombLog_imp*sqrt(1. + x_r*x_r) ); //Impurity scattering relaxation time in s
            }
            else tau_ei_imp = 1e100; //if Q = 0 (no impurities), set impurity scattering relaxation time to be extremely large so it does not affect result

            //Note: gsl_sf_expint_Ei = -int_{-x}^{inf}dt*exp(-t)/t
            if( w < 1e-15 ){
                Lambda_1 = 0.5*w*( 2. - 1./(1.+s) - 2.*s*log(1.+1./s) );
                Lambda_2 = v_F*v_F*0.5*w*( 1.5 - 3.*s - 1./(1.+s) + 3.*s*s*log(1.+1./s) );
            }
            else if( w > 100. ){
                Lambda_1 = 0.5*( log(1.+1./s) - 1./(1.+s) - 1./(sw*sw) );
                Lambda_2 = v_F*v_F*( (2.*s+1.)/(2.*s+2.) - s*log(1.+1./s) );
            }
            else if(sw < 10.){
                expi_sw = -gsl_sf_expint_Ei(-sw);
                expi_sww = -gsl_sf_expint_Ei(-sw-w);
                Lambda_1 = ( log(1.+1./s) + 1./(1.+1./s)*(1.-exp(-w)) - (1.+sw)*exp(sw)*( expi_sw - expi_sww ) )/2.;
                Lambda_2 = ( v_F*v_F*(exp(-w)-1.+w)/w - svF2/(1.+1./s)*(1.-exp(-w))
                            - 2.*svF2*log(1.+1./s) + svF2*(2.+sw)*exp(sw)*( expi_sw - expi_sww ) )/2.;
            }
            else{
                Lambda_1 = ( log(1.+1./s) + 1./(1.+1./s)*(1.-exp(-w)) - (1.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
                Lambda_2 = ( v_F*v_F*(exp(-w)-1.+w)/w - svF2/(1.+1./s)*(1.-exp(-w))
                            - 2.*svF2*log(1.+1./s) + svF2*(2.+sw)*( expE1_rp(sw)-exp(-w)*expE1_rp(sw+w) ) )/2.;
            }
            if( w1 < 1e-15 ){
                Lambda0_1 = 0.5*w1*( 2. - 1./(1.+s) - 2.*s*log(1.+1./s) );
                Lambda0_2 = v_F*v_F*0.5*w1*( 1.5 - 3.*s - 1./(1.+s) + 3.*s*s*log(1.+1./s) );
            }
            else if( w1 > 100. ){
                Lambda0_1 = 0.5*( log(1.+1./s) - 1./(1.+s) - 1./(sw1*sw1) );
                Lambda0_2 = v_F*v_F*( (2.*s+1.)/(2.*s+2.) - s*log(1.+1./s) );
            }
            else if(sw1 < 10.){
                expi_sw1 = -gsl_sf_expint_Ei(-sw1);
                expi_sw1w1 = -gsl_sf_expint_Ei(-sw1-w1);
                Lambda0_1 = ( log(1.+1./s) + 1./(1.+1./s)*(1.-exp(-w1)) - (1.+sw1)*exp(sw1)*( expi_sw1 - expi_sw1w1 ) )/2.;
                Lambda0_2 = ( v_F*v_F*(exp(-w1)-1.+w1)/w1 - svF2/(1.+1./s)*(1.-exp(-w1)) - 2.*svF2*log(1.+1./s) + svF2*(2.+sw1)*exp(sw1)*( expi_sw1 - expi_sw1w1 ) )/2.;
            }
            else{
                Lambda0_1 = ( log(1.+1./s) + 1./(1.+1./s)*(1.-exp(-w1)) - (1.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
                Lambda0_2 = ( v_F*v_F*(exp(-w1)-1.+w1)/w1 - svF2/(1.+1./s)*(1.-exp(-w1)) - 2.*svF2*log(1.+1./s) + svF2*(2.+sw1)*( expE1_rp(sw1)-exp(-w1)*expE1_rp(sw1+w1) ) )/2.;
            }

            // Thermal conductivity for electron-electron scattering. Note: u in Shternin and Yakovlev 2006 = v_F
            T_plasma_e = sqrt(4.*pi*alpha_e*hbarc*hbarc*hbarc*n_e[i][j]/mu_e); //electron plasma temperature in MeV
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
            Ath = an[0] + an[1]*th + an[2]*th2 + an[3]*th3 + an[4]*th4 + an[6]*th6+ an[8]*th8;
            Bth = bn[0] + bn[2]*th2 + bn[4]*th4 + bn[5]*th5 + bn[6]*th6 + bn[7]*th7 + alpha[3]*an[5]*th9 + alpha[4]*an[7]*th11;
            Athp = an[1] + 2.*an[2]*th + 3.*an[3]*th2 + 4.*an[4]*th3 + 6.*an[6]*th5 + 8.*an[8]*th7;
            Athpp = 2.*an[2] + 6.*an[3]*th + 12.*an[4]*th2 + 30.*an[6]*th4 + 56.*an[8]*th6;
            Bthp = 2.*bn[2]*th + 4.*bn[4]*th3 + 5.*bn[5]*th4 + 6.*bn[6]*th5 + 7.*bn[7]*th6 + 9.*alpha[3]*an[5]*th8 + 11.*alpha[4]*an[7]*th10;
            Bthpp = 2.*bn[2] + 12.*bn[4]*th2 + 20.*bn[5]*th3 + 30.*bn[6]*th4 + 42.*bn[7]*th5 + 72.*alpha[3]*an[5]*th7 + 110.*alpha[4]*an[7]*th9;
            C_vi = th2*( pow( alpha[0]/(exp(0.5*alpha[0]*th)-exp(-0.5*alpha[0]*th)),2. ) + pow( alpha[1]/(exp(0.5*alpha[1]*th)-exp(-0.5*alpha[1]*th)),2. )
                                     + pow( alpha[2]/(exp(0.5*alpha[2]*th)-exp(-0.5*alpha[2]*th)),2. ) + (Athpp*Bth*Bth-2.*Athp*Bthp*Bth+2.*Ath*Bthp*Bthp-Ath*Bth*Bthpp)/(Bth*Bth*Bth) ); //specific heat per ion divided by k_B (dimensionless)

            LogFact = log( 2. + 1./(sqrt(3.*Gamma)*Gamma) );
            kappa_0 = k_Bcgs*T_plasma*n_i[i]*a_i*a_i/hbarc*c*1e26; //normalization constant in units of erg/(cm*K*s). Numerical factor converts fm^-2 to cm^-2
            kappa_ii = kappa_0*sqrt( 16./(Gamma*Gamma*Gamma*Gamma*Gamma*LogFact*LogFact) + 0.16 + Gamma*Gamma/5929.*exp( 2./3.*T_plasma/T_MeV )  ); //ion-ion scattering part of ionic thermal conductivity in erg/(cm*K*s)

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
                CoulombLogVal = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_sigma*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                tau_ei = 5.699361922065672e-17/(Z[i]*CoulombLogVal*sqrt(1. + x_r*x_r)); //effective electron-ion relaxation time in s^{-1}
                tau_ei = tau_ei/( 1. + tau_ei/tau_ei_imp ); // corrects electron-ion relaxation times for impurity scattering
                sigma = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_ei; //returns electrical conductivity in s^{-1}. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                eta_O[i][j] = c*c*t_0/(4.*pi*sigma*L_0*L_0); //Ohmic diffusivity in reduced units
            }
            else if(cparams.EorT == "Thermal"){
                G_s_kappa = G0 + G2; //G_kappa. 1. + pow(x_nuc,2)*sqrt(2.*Z[i]) factor is from Gnedin et al. 2001.
                CoulombLogVal = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_kappa*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                tau_ei = 5.699361922065672e-17/(Z[i]*CoulombLogVal*sqrt(1. + x_r*x_r)); //effective electron-ion relaxation time in s^{-1}
                tau_ei = tau_ei/( 1. + tau_ei/tau_ei_imp ); // corrects electron-ion relaxation times for impurity scattering
                tau_e = tau_ei*tau_ee/( tau_ei + tau_ee ); //Includes electron-electron scattering in relaxation time for thermal conductivity
                kappa[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*tau_e )*t_0/(s_0*L_0*L_0); //electron-ion scattering thermal conductivity in reduced units. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                kappa[i][j] += kappa_i*t_0/(s_0*L_0*L_0); //adds ionic thermal conductivity to kappa.
                if( kappa[i][j]/c_v[i][j] > max_eta_T ) max_eta_T = kappa[i][j]/c_v[i][j]; //computes maximum value of thermal diffusivity
            }
            else if(cparams.EorT == "Both"){
                G_s_sigma = G0; //G_sigma
                CoulombLogVal = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_sigma*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                tau_ei = 5.699361922065672e-17/(Z[i]*CoulombLogVal*sqrt(1. + x_r*x_r)); //effective electron-ion relaxation time in s
                tau_ei = tau_ei/( 1. + tau_ei/tau_ei_imp ); // corrects electron-ion relaxation times for impurity scattering
                sigma = 1.2941734630093109e47*n_e[i][j]/mu_e*tau_ei; //returns electrical conductivity in s^{-1}. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                eta_O[i][j] = c*c*t_0/(4.*pi*sigma*L_0*L_0); //Ohmic diffusivity in reduced units

                G_s_kappa = G0 + G2; //G_kappa. Add correction factor to G_sigma. 1. + pow(x_nuc,2)*sqrt(2.*Z[i]) factor is from Gnedin et al. 2001.
                CoulombLogVal = ( (Lambda_1 - Lambda_2) - (Lambda0_1 - Lambda0_2) )*G_s_kappa*D; //Fitted formula for Coulomb logarithm from Gnedin et al 2001. Note we absorb the (v_F/c)^2 factor into the definition of Lambda_2
                tau_ei = 5.699361922065672e-17/(Z[i]*CoulombLogVal*sqrt(1. + x_r*x_r)); //effective electron-ion relaxation time in s
                tau_ei = tau_ei/( 1. + tau_ei/tau_ei_imp ); // corrects electron-ion relaxation times for impurity scattering
                tau_e = tau_ei*tau_ee/( tau_ei + tau_ee ); //Includes electron-electron scattering in relaxation time for thermal conductivity
                kappa[i][j] = ( n_e[i][j]/mu_e*T_MeV*kappa_pref/3.*tau_e )*t_0/(s_0*L_0*L_0); //electron-ion scattering thermal conductivity in reduced units. Numerical factors here are converting n_e to cm^{-3} and mu_e to erg/c^2 = g
                kappa[i][j] += kappa_i*t_0/(s_0*L_0*L_0); //adds ionic thermal conductivity to kappa.
                if( kappa[i][j]/c_v[i][j] > max_eta_T ) max_eta_T = kappa[i][j]/c_v[i][j]; //computes maximum value of thermal diffusivity
            }
            else{
                std::cout << "Invalid value for EorT in sigmaCalc" << std::endl;
            }
        }
    }

    //Mirrors kappa across lower boundary if it is being computed
//    if(cparams.EorT == "Thermal" || cparams.EorT == "Both"){
//        for(size_t i=0; i<N_GC; i++){
//            for(size_t j=0; j<eta_O.shape()[1]; j++){
//                kappa[N_GC-1-i][j] = kappa[N_GC+i][j];
//                eta_O[N_GC-1-i][j] = eta_O[N_GC+i][j];
//            }
//        }
//    }

    return;
}

/*
    Computes exponential integral exp(x)*E1(x)=exp(x)*\int_x^{inf}e^{-y}/y dy using rational-polynomial approximation
    Note that E1(x) = -Ei(-x)
    Eq. 5.1.55 and 5.1.56 from Abramowitz and Stegun, Handbook of Mathematical Functions (pg. 231)
*/
double expE1_rp(double x){

    double expE1;

    static double a_1 = 8.5733287401;
    static double a_2 = 18.0590169730;
    static double a_3 = 8.6347608925;
    static double a_4 = 0.2677737343;
    static double b_1 = 9.5733223454;
    static double b_2 = 25.6329561486;
    static double b_3 = 21.0996530827;
    static double b_4 = 3.9584969228;

    static double c_0 = -0.57721566;
    static double c_1 = 0.99999193;
    static double c_2 = -0.24991055;
    static double c_3 = 0.05519968;
    static double c_4 = -0.00976004;
    static double c_5 = 0.00107857;

    if( x > 1. ){
        expE1 = ( pow(x,4.) + a_1*pow(x,3.) + a_2*pow(x,2.) + a_3*x + a_4 )/( pow(x,4.) + b_1*pow(x,3.) + b_2*pow(x,2.) + b_3*x + b_4 )/x;
    }
    else{
        expE1 = exp(x)*( -log(x) + c_0 + c_1*x + c_2*pow(x,2.) + c_3*pow(x,3.) + c_4*pow(x,4.) + c_5*pow(x,5.) );
    }

    return expE1;
}

/*
    Computes impurity parameter Q=Z_imp^2 based on model Q100 of Vigano 2013 thesis (Fig. 5.10)
    Inputs: n_b: baryon number density in fm^-3
    Output: Q_imp (dimensionless) for a given n_b
 */
double Q_impVigano2013Q100(double n_b){

    static double n_bPasta1 = 0.0060; //lower baryon number density edge of pasta phase in fm^-3
    static double n_bPasta2 = 0.0472; //upper baryon number density edge of pasta phase in fm^-3

    double Q;

    if(n_b < n_bPasta1){
        Q = 0.1;
    }
    else if(n_b < n_bPasta2){
        Q = 0.1 + (100. - 0.1)*pow( ( n_b - n_bPasta1 )/( n_bPasta2 - n_bPasta1 ), 2. );
    }
    else{
        Q = 100.;
    }

    return Q;
}

/*
    Computes impurity parameter Q=Z_imp^2 based on calculation for BSk24 EoS from
    Carreau, Fantina and Gulminelli A&A 640, A77 (2020). Uses a fitting formula for both the outer and inner crust calculations
    (Fig. 5 and 6 respectively of Carreau+2020), smoothing over the calculated outer crust profile.

    Inputs: n_b: baryon number density in fm^-3
            n_drip: neutron drip density in fm^-3
    Output: Q_imp (dimensionless) for a given n_b
 */
double Q_impCarreau2020BSk24(double n_b, double n_drip){

    static double a[4] = {91.1311409,  0.518097003, -198.927735,  0.0483712375}; //fit parameters for n_b < n_drip
    static double b[8] = {0.163713777,  45.4854712,  644.814328, -3853.10412, 0.0684925928,  7.28712434, -31.7970860, -11.2961386}; //fit parameters for n_b > n_drip

    static double n_s = 0.16; //nuclear saturation density in fm^-3
    double u = n_b/n_s, Q;

    if(n_b < n_drip){
        Q = a[0]*pow(u,a[1]) + a[2]*u + a[3];
    }
    else if(n_b > n_drip){
        Q = ( b[0] + b[1]*u + b[2]*u*u + b[3]*u*u*u )/( b[4] + b[5]*u + b[6]*u*u + b[7]*u*u*u );
    }

    return Q;
}

/*
    Computes specific heat capacity at B=0 in reduced units. Includes electron, ion and free neutron contributions

    Thermal harmonic part of ionic heat capacity c_vi_h from Baiko, Potekhin and Yakovlev PRE 64, 057402 (2001), assuming bcc lattice
    Thermal anharmonic part of ionic heat capacity c_vi_ah from Baiko and Chugunov, MNRAS 510, 2628–2643 (2022)
    Electron polarization term c_vei from Potekhin and Chabrier Contrib. Plasma Phys. 50, 82 (2010) and Potekhin and Chabrier A&A, 550, A43 (2013) Appendix C.2
    Neutron heat capacity c_vn from Pastore, Chamel and Margueron, MNRAS 448, 1887 (2015)
    Neutron 1S0 pairing gap in crust parameterization taken from model of Kaminker, Haensel and Yakovlev A&A 373, L17 (2001)
        and Andersson, Comer and Glampedakis Nucl. Phys. A 763, 212 (2005)
    Specific pairing gap model taken from Schwenk, Friman and Brown Nucl. Phys. A 713, 191 (2003), SFB parametrization in Ho et al. PRC 91, 015806 (2015) Table II.
    Neutron 1S0 superfluidity control function R_A(T/T_cn) taken from Eq. (18), R_A of Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999)

    Inputs: T: redshifted temperature in reduced units
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
            c_v: ScalarField object containing specific heat capacity of magnetized electrons (B-dependent) in reduced units, computed by "Omega_xyVaryingTC_veOnly"
    Output: c_v: ScalarField object containing specific heat capacity in reduced units
*/
void c_vfuncCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, RadialScalarField & n_nf, RadialScalarField & mu_nf, const TParams & tparams, size_t N_GC, std::vector<double> & n_t, ScalarField & c_v)
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
    double k_Fn, E_Fn, T_cl, pref, x_r, mu_e, c_ve, c_vi_h, c_vi_ah, c_vei, c_vn;
    //parametrization of 1S0 neutron pairing gap. Delta0 in MeV, k0, k2 in fm^{-1}, k1, k3 in fm^{-2}
    double Delta0 = tparams.Delta0_nCrust, k0 = tparams.k0_nCrust, k1 = tparams.k1_nCrust, k2 = tparams.k2_nCrust, k3 = tparams.k3_nCrust;

    double kF0, kF2, Delta_n, T_cn, tau, v_A;
    double R_A = 1.; //superfluidity control function for neutron specific heat capacity. Taken to be 1 if no superfluidity.

    // std::fstream Dataout;
    // Dataout.open("Testing/ThermalPhysicsTesting/SpecificHeatCapacities/CrustCv_T5e8K.dat",std::fstream::out);
    // Dataout << "n_b (fm^-3)   c_ve (erg/K/cm^3)   c_vi (erg/K/cm^3)   c_vn (erg/K/cm^3)" << std::endl;
    for(size_t i=0; i<c_v.shape()[0]-N_GC+1; i++){

        a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm
        s = 1./( 1. + 0.01*pow( log(Z[i]),1.5 ) + 0.097/(Z[i]*Z[i]) );
        b_ei1 = 1. - 1.1866*pow( Z[i],-0.267 ) + 0.27/Z[i];
        b_ei2 = 1. + 2.25/pow( Z[i],1./3. )*( 1. + 0.684*pow( Z[i],5. ) + 0.222*pow( Z[i],6. ) )/( 1. + 0.222*pow( Z[i],6. ) );
        b_ei3 = 41.5/( 1. + log(Z[i]) );
        b_ei4 = 0.395*log(Z[i]) + 0.347/pow( Z[i],1.5 );

        for(size_t j=0; j<c_v.shape()[1]; j++) {
            if(i < N_GC){
                T_MeV = T[N_GC+1-i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            else{
                T_MeV = T[i][j]/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
            }
            mu_e = sqrt( pow(3.*pi*pi*n_e[i][j],2./3.)*hbarc*hbarc+M_e*M_e );
            x_r = pow(3.*pi*pi*n_e[i][j],1./3.)*hbarc/M_e;
            c_ve = mu_e*sqrt( mu_e*mu_e - M_e*M_e )*T_MeV/(3.*hbarc*hbarc*hbarc)*1e39*k_Bcgs; //in erg/K/cm^3

            Zeq = n_e[i][j]/n_i[i]; //proton number in unit cell. Only differs from Z near crust-core transition
            T_p = sqrt(T_pconst*Z[i]*Zeq*n_i[i]/A[i]); //plasma temperature in MeV
            // T_p = sqrt(T_pconst*Z[i]*Zeq*n_i[i]/( A[i] + 0.8*(n_b[i]/n_i[i]-A[i]) ) ); //plasma temperature in MeV
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
            Bthpp = 2.*bn[2] + 12.*bn[4]*th2 + 20.*bn[5]*th3 + 30.*bn[6]*th4 + 42.*bn[7]*th5 + 72.*alpha[3]*an[6]*th7 + 110.*alpha[4]*an[8]*th9;
            c_vi_h = n_i[i]*1e39*k_Bcgs*th2*( pow(alpha[0]/(exp(0.5*alpha[0]*th)-exp(-0.5*alpha[0]*th)),2.) + pow(alpha[1]/(exp(0.5*alpha[1]*th)-exp(-0.5*alpha[1]*th)),2.)
                         + pow(alpha[2]/(exp(0.5*alpha[2]*th)-exp(-0.5*alpha[2]*th)),2.) + (Athpp*Bth*Bth-2.*Athp*Bthp*Bth+2.*Ath*Bthp*Bthp-Ath*Bth*Bthpp)/(Bth*Bth*Bth) ); //in erg/K/cm^3

            A1 = -2.*A11/th*( 6.*A12*A12*th4 + 3.*A12*th2 + 1. )/pow( 1. + A12*th2,3. ) - 2.*A13/th*( 6.*A14*A14*th4 + 3.*A14*th2 + 1. )/pow( 1. + A14*th2,3. );
            A2 = 1.5*Acl[1]/th2*( 2. + 3.*A21*th4 )/pow( 1. + A21*th4,1.25 );
            A3 = 4.*Acl[2]/th3;
            c_vi_ah = n_i[i]*1e39*k_Bcgs*( A1/Gamma_q + A2/(Gamma_q*Gamma_q) + A3/(Gamma_q*Gamma_q*Gamma_q) ); //in erg/K/cm^3

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
            E_Fn = mu_nf[i]-M_n;//sqrt( k_Fn*k_Fn*hbarc*hbarc + M_n*M_n ) - M_n; //mu_nf-M_n in MeV
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
                        v_A = sqrt( 1. - tau )*( 1.456 - 0.157/sqrt(tau) + 1.764/tau ); //ratio of temperature-dependent gap to T, from Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999) Eq. (11)
                        R_A = pow( 0.4186 + sqrt( 1.014049 + 0.251001*v_A*v_A ),2.5 )*exp( 1.456 - sqrt( 2.119936 + v_A*v_A ) );
                    }
                    else R_A = 1.;
                }
                else{
                    R_A = 1.;
                }
            }
            // c_v[i][j] = (c_ve + c_vi_h + c_vi_ah + c_vei + c_vn)/s_0;
            c_v[i][j] = (c_ve + c_vi_h + c_vi_ah + c_vei + R_A*c_vn)/s_0;
        }
    }

    // Dataout.close();

    return;
}

/*
            Crust neutrino emissivity in reduced units

    Plasmon emissivity from Kantor and Gusakov, MNRAS 381, 1702 (2007)
    Electron-nucleus neutrino bremsstrahlung emissivity from Ofengeim et al, EPL 106, 31002 (2014)
    Neutrino synchrotron radiation emissivity from Bezchastnov et al, A&A 328, 409 (1997)
    Electron-positron pair annihilation emissivity (not currently implemented) from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (22,23,26,27), assuming degenerate electrons.
    Use Pearson et al, 2018, MNRAS 481, 2994, Eq. (C21) and Table C10 for computing R_p=sqrt(3/5)*x_p*a_i above neutron drip (x_p = x_nuc in ConductivityCalc) for BSk EOS.
    Note: Eq. (C21) has a typo in the denominator of the first term: 1 + p_4... should be 1 - p_4... (checked via their EOS code). Also replaced p_6 with 0.325 to improve fit
    Note: if you compare to Fig. 3 of Potekhin, Pons and Page, Space. Sci. Rev., 191, 239 (2015)
          there results were computed for T=2e9 K, not 1e9 K as labelled

    Inputs: T: redshifted temperature in reduced units
            n_e: electron number density in fm^{-3}
            Bmag: magnetic field magnitude in reduced units
            A: mass number of nuclei (excludes dripped neutrons and protons)
            Z: atomic number of nuclei (excludes dripped protons)
            n_i: ion number density in fm^{-3}
            n_b: baryon number density in fm^{-3}
            cparams: ConductParams object which conductivities to compute, which impurity model to use, whether to include quantization effects, neutron drip density
            tparams: TParams object containing information about core temperature.
            N_GC: number of ghost cells
            n_t: lapse function
            IMEX: true if using an IMEX timestepping method and require dq_nuCrustdT.
    Output: q_nuCrust: ScalarField object containing crust neutrino emissivity in reduced units
            dq_nuCrustdT: ScalarField object containing derivative of crust neutrino emissivity w.r.t. local temperature in reduced units
*/
void q_nuCrustCalc(ScalarField & T, ScalarField & n_e, ScalarField & Bmag, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC, std::vector<double> & n_t,
                    bool IMEX, ScalarField & q_nuCrust, ScalarField & dq_nuCrustdT)
{

    static double p[10] = {1.793, 0.0645, 0.433, 0.01139, 2.484e6, -0.6195, 9.632e-4, 0.4372, 1.614, 8.504e8};
    static double s[3] = {9.079, 1.399, -0.06592};
    static double r[3] = {0.3520, 1.195, -0.1060};
    static double q[5] = {0.7886, 0.2642, 1.024, 0.07839, 0.1784};
    static double zeta_3 = 1.20205690316; //Riemann zeta function zeta(3)
    static double zeta_5 = 1.03692775514; //Riemann zeta function zeta(5)
    static double b = 231., un1 = 2.798, un2 = 12.973, alpha = 0.23;
    static double Q_0 = 1.3858e21; //emissivity scale for plasmon decay in erg/s/cm^3
    static double Cm2 = 0.1748, Cp2 = 1.675; // C_-^2 and C_+^2, the difference/sum between the squares of the vector and axial vector constants squared
    static double amu = 931.49410242; //1 amu in MeV
//    static double Q_c = 1.023e23; //emissivity scale for pair annihilation in erg/s/cm^3
//    static double pk0[3] = {23.61, 32.11, 16.};
//    static double pk1[5] = {42.44, 140.8, 265.2, 287.9, 144.};
//    static double pk2[7] = {61.33, 321.9, 1153., 2624., 4468., 4600., 2304.};

    static double n_drip = cparams.n_b_nd; //neutron drip density in fm^{-3} from equation of state

    double T_MeV, a_i, rhobar, p_F, tbar, asy1, asy2, c1, c2;
    double C, d1, d2, D, f, v_F, v_star, beta;
    double beta6, a_1, a_2, b_1, b_2, f6, f7_5, asyt1, asyt2, asyl1, asyl2, W_t, W_l;
    double n_i34, R_p, xi, Gamma, tp, Lambda, tau_s, G, H, Gamma_s, Bconst, w_s, w_1s, R_NB, Lfit;
    double T_B, T_P, z, y_1, y_2, F_p_denom_factor, F_p, F_m_denom_factor, F_m, B13, T9, S_AB, S_BC;
    double q_plasmon = 0., q_bremsstrahlung = 0., q_synchrotron = 0.;
//    double q_pair = 0., y_r, U_m1, U_0, U_1, U_2, Phi_m1, Phi_0, Phi_1, Phi_2, tbar2, tbar3, tbar4, tbar5, tbar6, tbar7;

    double DeltaT_MeV = 1e-4; //fraction of T to use as denominator for finite-difference derivative of q_nuCrust
    //if using an IMEX timestepperm, will compite q_nuCrust twice with slightly different T and then compute dq_nuCrustdT using a finite difference
    size_t kmax = 0;
    if( IMEX == true ) kmax = 1;

    for(size_t i=0; i<q_nuCrust.shape()[0]-N_GC; i++){

        for(size_t j=0; j<q_nuCrust.shape()[1]; j++){

            for(size_t k=0; k<=kmax; k++){
                if(i < N_GC){
                    T_MeV = T[N_GC+1-i][j]*(1.+double(k)*DeltaT_MeV)/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
                }
                else{
                    T_MeV = T[i][j]*(1.+double(k)*DeltaT_MeV)/n_t[i]*T_0*k_B; //local temperature in MeV, converted from reduced units
                }
                T9 = T_MeV/k_B/1e9; //temperature in units of 1e9 K
                a_i = pow(3./(4.*pi*n_i[i]),1./3.); //ion sphere radius in fm

                rhobar = n_e[i][j]*1.e39*amu*MeVtoErg/(c*c); //effective mass density in g/cm^3 (Eq. (4) in Kantor and Gusakov 2007)
                p_F = hbarc/M_e*pow( 3.*pi*pi*n_e[i][j],1./3. ); //dimensionless Fermi momentum. Called x_r in some references

                tbar = T_MeV/M_e; //dimensionless ratio of thermal energy over electron mass

                asy1 = 4.*alpha_e/(3.*pi)*p_F*p_F*p_F/sqrt(1.+p_F*p_F);
                asy2 = 4.*pi*alpha_e/9.*p[1]*( tbar*tbar/p[1] + 1. + p[1]/(tbar*tbar) )*pow( 1. + p[2]/pow( tbar/sqrt(p[1]),p[0] ),-10. );
                c1 = p[3]*( 1. + p[4]*pow( rhobar,p[5] ) )/( 1. + p[6]*( 1 + p[4]*pow( rhobar,p[5] ) ) );
                c2 = p[7] + p[8]*rhobar/( p[9] + rhobar );

                C = 1. - c2*pow( c1*tbar,2. )/( 1. + pow(c1*tbar,2.) );
                d1 = 6./(pi*pi)*p_F*p_F*( 1. + p_F*p_F )/( 2.*p_F*p_F + 5. );
                d2 = pi*pi/6./( ( sqrt(1. + p_F*p_F) - 1. ) );
                D = tbar*tbar/( d1*sqrt( 1. + pow( d2*tbar,2. ) ) );
                f = pow( asy2*asy2 + pow( asy1*(1.-C*D),2. ),0.25 )/tbar; //plasma frequency times hbar divided by electron mass (dimensionless)
                v_F = p_F/sqrt( 1. + p_F*p_F );
                v_star = pow( ( v_F*v_F*v_F + s[0]*pow(tbar,s[1])*pow(rhobar,s[2]) )/( 1. + s[0]*pow(tbar,s[1])*pow(rhobar,s[2]) ),1./3. );
                beta = sqrt( 1.5/(v_star*v_star)*( 1. - (1.-v_star*v_star)/(2.*v_star)*log( (1.+v_star)/(1.-v_star) ) ) );

                beta6 = beta*beta*beta*beta*beta*beta + ( 3.375 - beta*beta*beta*beta*beta*beta )*pow( tbar,r[1] )*pow( rhobar,r[2] )/( r[0] + pow( tbar,r[1] )*pow( rhobar,r[2] ) ); //beta to the power of 6 using Eq. (31) from Kantor and Gusakov 2007
                a_1 = 4.*zeta_3*beta6;
                a_2 = 8./105. + ( 0.349 - 8./105. )*pow( v_star,10. );
                b_1 = sqrt(2.*pi)*pow( 1. + v_star*v_star/5.,-1.5 );
                b_2 = sqrt(0.5*pi)*pow( 3.*v_star*v_star/5.,-1.5 );
                f6 = f*f*f*f*f*f;
                asyt1 = a_1*f6;
                f7_5 = f6*f*sqrt(f);
                asyt2 = b_1*f7_5;
                asyl1 = a_2*f6*f*f;
                asyl2 = b_2*f7_5;
                W_t = asyt1 + asyt2*exp( q[2]/( pow(f,q[0]) + q[1] ) );
                W_l = ( asyl2*( asyl1 + q[3]*pow( 1. + q[4]*pow(v_star,2.5),3.5 )*pow(f,9.) ) )/( asyl2 + ( asyl1 + q[3]*pow( 1. + q[4]*pow(v_star,2.5),3.5 )*pow(f,9.) ) );
                q_plasmon = -Q_0*pow(tbar,9.)*( W_t + W_l )*exp(-f); //plasmon decay neutrino emissivity in erg/s/cm^3

                n_i34 = n_i[i]*1.e5; //ion number density in units of 1e34 cm^{-3}
                //Proton radius R_c = R_p = in fm given n_b in fm^{-3}. Equal to sqrt(3/5)*x_nuc*a_i.
                if(n_b[i] < n_drip) R_p = 1.15*pow(A[i],1./3.); //From Kaminker et al., Astron. Astrophys. 343, 1009 (1999). x_nuc = 1.15*pow(A[i],1./3.)/a_i
                else if(cparams.EOS == "BSk24") R_p = a_i*( 0.1035 + 1.944*pow( n_b[i],0.5717 ) )/( 1. - 608.*pow( n_b[i],3.143 ) ) + 0.0225*Z[i]*pow( n_b[i],1.26 ); //From Pearson et al. MNRAS 481, 2994 (2018) for BSk24 EOS. Eq. (C21) and Table C10. Use Z instead of Z_eq, only minor difference
                else if(cparams.EOS == "SLy4") R_p = a_i*( 0.08633 + 1.0173*pow( n_b[i],0.4857 ) )/( 1. - 237.6*pow( n_b[i],2.615 ) ) + 0.2732*Z[i]*pow( n_b[i],0.7222 ); //Pearson et al. 2018 formula fitted to Douchin and Haensel 2001 data
                else R_p = a_i*( 0.1035 + 1.944*pow(n_b[i],0.5717) )/( 1. - 608.*pow(n_b[i],3.143) ) + 0.0225*Z[i]*pow(n_b[i],1.26); //Uses BSk24 formula by default

                xi = 0.06665*R_p*pow( Z[i]*(n_i34),1./3. );
                Gamma = 0.5798*Z[i]*Z[i]/(T_MeV/k_B/1.e9)*pow( n_i34,1./3. );
                tp = 0.9914*(T_MeV/k_B/1.e9)/Z[i]*sqrt(A[i]/n_i34); //called tau in Ofengeim et al 2014
                Lambda = 1 + Z[i]*Z[i]/(35.5*35.5)/( 1. + Gamma/(223.*Z[i]) );
                tau_s = 0.095*pow( 2.*tp/(0.095*sqrt(tp*tp + 4.)),1./Lambda );
                G = 0.5 + 0.002*Z[i]*xi;
                H = ( 1.52 + 0.9*(1.-tau_s) )/( 1. + 0.5*xi );
                Gamma_s = Gamma*exp( -pow(tau_s,0.053)*pow( Z[i]/11.3,4./9. )/( 1. + pow(tau_s,G)*pow( Gamma/19.3/pow(Z[i],1.7),H ) ) );
                Bconst = 0.6*( pow(12.*pi*pi,1./3.)*un2*pow(Z[i],0.8)/( 1. + (xi/Z[i])*Gamma_s/( 200. + Gamma_s ) ) )/sqrt( Gamma_s*Gamma_s + 1037./pow( 1. + Gamma_s/204.,4. )*pow( Gamma_s,0.5+0.075*xi )/pow(Z[i],0.1) );
                w_s = Bconst*( 1. + un1/(2.*un2*tau_s)*exp( -9.1*(tau_s + 216./Gamma_s) ) );
                w_1s = Bconst*b*tau_s/sqrt( b*b*tau_s*tau_s + un2*un2*exp(-7.6*(tau_s + 18./Gamma_s)) );
                R_NB = ( 1. + 5.54e-3*Z[i] + 7.37e-5*Z[i]*Z[i] ); //non-Born correction

                Lfit = CoulombLogf( w_s + 4.*alpha*xi*xi ) - CoulombLogf( w_s - w_1s + 4.*alpha*xi*xi );
                q_bremsstrahlung = -5.362e15*Z[i]*Z[i]*n_i34*T9*T9*T9*T9*T9*T9*R_NB*Lfit; //neutrino bremsstrahlung emissivity in erg/s/cm^3

                B13 = Bmag[i][j]*B_0/1e13; //magnetic field magnitude in units of 1e13 G
                T_B = 1.34e9*B13/sqrt(1.+p_F*p_F)*k_B; //in MeV  //unit_e*Bmag[i][j]*B_0*statCGtoMeV2/(k_B*sqrt(1.+x_r*x_r)); //
                T_P = 1.5*T_B*p_F*p_F*p_F; //in MeV
                z = T_B/T_MeV;
                xi = T_P/T_MeV;
                y_1 = pow( pow( 1. + 3172.*pow( xi,2./3. ),2./3. ) - 1.,1.5 );
                y_2 = pow( pow( 1. + 172.2*pow( xi,2./3. ),2./3. ) - 1.,1.5 );
                F_p_denom_factor = 1. + 2.036e-4*y_1 + 7.405e-8*y_1*y_1;
                F_m_denom_factor = 1. + 3.356e-3*y_2 + 1.536e-5*y_2*y_2;
                F_p = 44.01*( 1. + 3.675e-4*y_1 )*( 1. + 3.675e-4*y_1 )/( F_p_denom_factor*F_p_denom_factor*F_p_denom_factor*F_p_denom_factor );
                F_m = 36.97*( 1. + 1.436e-2*y_2 + 1.024e-5*y_2*y_2 + 7.647e-8*y_2*y_2*y_2 )/( F_m_denom_factor*F_m_denom_factor*F_m_denom_factor*F_m_denom_factor*F_m_denom_factor );
                S_AB = 27.*xi*xi*xi*xi/(512.*pi*pi*zeta_5)*( F_p - (Cm2/Cp2)*F_m ); // 0.175/1.675 is the ratio of C_-^2/C_+^2
                S_BC = exp( -0.5*z )*( 1. + 0.4228*z + 0.1014*z*z + 6.240e-3*z*z*z )/( 1. + 0.4535*pow( z,2./3. ) + 0.03008*z - 0.05043*z*z + 0.004314*z*z*z );

                q_synchrotron = -9.04e14*B13*B13*T9*T9*T9*T9*T9*S_AB*S_BC; //neutrino synchrotron radiation emissivity in erg/cm^3/s

    //            Assume electrons are degenerate
    //            y_r = sqrt( 1. + p_F*p_F ); //dimensionless electron chemical potential (p_F is dimensionless electron Fermi momentum)
    //            U_m1 = ( y_r*p_F - log(p_F+y_r) )/(2.*pi*pi);
    //            U_0 = p_F*p_F*p_F/(3.*pi*pi);
    //            U_1 = ( y_r*p_F*(p_F*p_F+y_r*y_r) - log(p_F+y_r) )/(8.*pi*pi);
    //            U_2 = p_F*p_F*p_F*( 5. + 3.*p_F*p_F )/(15.*pi*pi);
    //            tbar2 = tbar*tbar;
    //            tbar3 = tbar2*tbar;
    //            tbar4 = tbar3*tbar;
    //            tbar5 = tbar4*tbar;
    //            tbar6 = tbar5*tbar;
    //            tbar7 = tbar6*tbar;
    //            Phi_m1 = tbar/(2.*pi*pi)*sqrt( tbar*(2.*pi+7.662*tbar+1.92*tbar2)/(1.+0.48*tbar) )*exp( -(1.+y_r)/tbar );
    //            Phi_0 = tbar/(2.*pi*pi)*sqrt( 2.*pi*tbar + pk0[0]*tbar2 + pk0[1]*tbar3 + pk0[2]*tbar4 )*exp( -(1.+y_r)/tbar );
    //            Phi_1 = tbar/(2.*pi*pi)*sqrt( 2.*pi*tbar + pk1[0]*tbar2 + pk1[1]*tbar3 + pk1[2]*tbar4 + pk1[3]*tbar5 + pk1[4]*tbar6 )*exp( -(1.+y_r)/tbar );
    //            Phi_2 = tbar/(2.*pi*pi)*sqrt( 2.*pi*tbar + pk2[0]*tbar2 + pk2[1]*tbar3 + pk2[2]*tbar4 + pk2[3]*tbar5 + pk2[4]*tbar6
    //                                            + pk2[5]*tbar7 + pk2[6]*tbar7*tbar )*exp( -(1.+y_r)/tbar );
    //
    //            q_pair = -Q_c/(36.*pi)*( Cp2*( 8.*(Phi_1*U_2 + Phi_2*U_1) - 2.*(Phi_m1*U_2 + Phi_2*U_m1) + 7.*(Phi_0*U_1 + Phi_1*U_0) + 5.*(Phi_0*U_m1 + Phi_m1*U_0) )
    //                                    + 9.*Cm2*( Phi_0*(U_1 + U_m1) + (Phi_m1 + Phi_1)*U_0 ) ); //neutrino synchrotron radiation emissivity in erg/cm^3/s

                //Define each of q_plasmon, q_bremsstrahlung and q_synchrotron to be negative so include in generalized heat equation as + (q_plasmon + q_bremsstrahlung + q_synchrotron)
                if(k == 0){
                    q_nuCrust[i][j] = (q_plasmon + q_bremsstrahlung + q_synchrotron)*t_0/(T_0*s_0); //converts to reduced units
                }
                else if(k == 1){
                    dq_nuCrustdT[i][j] = ( (q_plasmon + q_bremsstrahlung + q_synchrotron)*t_0/(T_0*s_0) - q_nuCrust[i][j] )/(DeltaT_MeV*T_MeV/(T_0*k_B));
                }
            }
        }
    }


    return;
}

/*
    Function appearing in fit formula for Coulomb logarithm used by Ofengeim et al 2014 for neutrino brehmsstrahlung emissivity
*/
double CoulombLogf(double x){
    return 0.5*( sqrt(pi/x) * erf( sqrt(x) ) + log(x) + exp(-x)*expE1_rp(x) );
}

/*
    Error function rational approximation
    Eq. 7.1.26 from Abramowitz and Stegun, Handbook of Mathematical Functions (pg. 299)
*/
double erf(double x){
    double sign, t, y;
    static double a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741, a4 = -1.453152027, a5 =  1.061405429, p = 0.3275911;

    // save the sign of x
    if( x >= 0. ) sign = 1.; // 1 if x >= 0 else -1
    else sign = -1.;
    x = abs(x);

    t = 1./(1. + p*x);
    y = 1. - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    return sign*y; // erf(-x) = -erf(x)
}

/*
        ThermPhysCore: Computes volume-integrated core specific heat capacity at fixed volume and
                        neutrino emissivity in reduced units evaluated at uniform redshifted temperature T. Also computes
                        thermal conductivity of electrons and neutrons at T at the core side
                        of the crust-core transition. Only need to invoke
                        this function once per simulation.

        Inputs: CoreEOS_filename: core EOS file name, including file extension
                GR: true if general relativity is turned on and false otherwise
                tparams: TParams object containing information about nucleon superfluid gaps.
                T: uniform redshifted temperature at which to compute core properties in reduced units
        Output: C_v: total heat capacity of core at T in reduced units
                Q_nu: total neutrino emissivity of core at T in reduced units
*/
void ThermPhysCore(std::string CoreEOS_filename, bool GR, const TParams & tparams, double T, double & C_v, double & Q_nu, double & kappa_eCC, double & kappa_nCC)
{

    double km_to_cm = 1e5; //converts radius in km to cm

	//Determines number of rows in data table
	std::ostringstream filename;
	filename << CoreEOS_filename;
	std::fstream Datafile;
    Datafile.open(filename.str().c_str(),std::fstream::in);
	std::string line, dummyLine;
    size_t N_dat = 0; //number of rows in data table. Initialize to zero
	getline(Datafile, dummyLine); //skips header (1 line)
	while (std::getline(Datafile, line))
		++N_dat;

	Datafile.close();

    //vectors to store values from data table
    std::vector<double> r_Vec; //radius in reduced units. Vector
    std::vector<double> n_t_Vec; //lapse function n_t = sqrt(-g_{tt}) = e^{nu/2}. Vector
    std::vector<double> n_r_Vec; //n_r = sqrt(g_{rr}) = e^{lambda/2}. Vector
    std::vector<double> n_b_Vec; //Baryon number density in fm^{-3}. Vector
    std::vector<double> n_e_Vec; //Electron number density in fm^{-3}. Vector
    std::vector<double> mu_n_Vec; //Free (dripped) neutron chemical potential in MeV. Vector
    std::vector<double> Y; //proton fraction of total baryons
    std::vector<double> Y_e; //electron fraction of total baryons
    std::vector<double> M_nEff; //neutron Landau effective mass in MeV
    std::vector<double> M_pEff; //proton Landau effective mass in MeV

    //temporary variables to hold data table values
    double r_dt, M_dt, gtt_dt, grr_dt, n_b_dt, M_nEff_dt, M_pEff_dt, n_e_dt, n_n_dt, mu_n_dt, Aeq_dt;
    double dummy; //dummy variable to avoid saving all quantities from loaded data table
//    double p_Fn, p_Fp, mu_e; //neutron Fermi momentum in MeV

	Datafile.open(filename.str().c_str(),std::fstream::in);
	std::getline(Datafile, dummyLine); //skips header (1 line)
	//Gets background TOV solution values from data table
	for(size_t i=0; i<N_dat; i++)
	{
		Datafile >> r_dt >> dummy >> M_dt >> dummy >> gtt_dt >> grr_dt >> n_b_dt >> Aeq_dt >> M_nEff_dt >> dummy >> M_pEff_dt >> n_e_dt >> n_n_dt >> mu_n_dt;

		r_Vec.push_back(r_dt*km_to_cm/L_0);
		if(GR == true){
            n_t_Vec.push_back( sqrt(-gtt_dt) );
            n_r_Vec.push_back( sqrt(grr_dt) );
        }
        else{
            n_t_Vec.push_back( 1. );
            n_r_Vec.push_back( 1. );
        }
		n_b_Vec.push_back(n_b_dt);
		n_e_Vec.push_back(n_e_dt);
		mu_n_Vec.push_back(mu_n_dt);
		Y.push_back( (n_b_dt - n_n_dt)/n_b_dt ); //n_b - n_n = n_p
		if( Aeq_dt < 1e-10 ) Y_e.push_back( n_e_dt/n_b_dt ); // Y_e = n_e/n_b
		else Y_e.push_back( (n_b_dt - n_n_dt)/n_b_dt );
		M_nEff.push_back( M_nEff_dt );
		M_pEff.push_back( M_pEff_dt );

	}
	Datafile.close();

    //Compute total core volume-integrated heat capacity and neutrino emissivity at T=T_0/e^{\nu/2} by volume integration. Use trapezoid rule.
    C_v = 0.;
    double Q_nuMUrca = 0., Q_nuDUrca = 0., Q_nu_NNbrems = 0., Q_PBF = 0.;
    double dr; //radial coordinate increment for integration

    for(size_t i=1; i<N_dat; i++){

        dr = (r_Vec[i]-r_Vec[i-1]);
        C_v += 4.*pi*dr*0.5*( pow(r_Vec[i-1],2.)*n_r_Vec[i-1]*c_vCore(T*T_0/n_t_Vec[i-1],n_b_Vec[i-1],Y[i-1],Y_e[i-1],M_nEff[i-1],M_pEff[i-1],tparams)
                              + pow(r_Vec[i],2.)*n_r_Vec[i]*c_vCore(T*T_0/n_t_Vec[i],n_b_Vec[i],Y[i],Y_e[i],M_nEff[i],M_pEff[i],tparams) );
        Q_nuMUrca += 4.*pi*dr*0.5*( pow(r_Vec[i-1],2.)*n_r_Vec[i-1]*pow(n_t_Vec[i-1],2.)*q_nuMUrca(T*T_0/n_t_Vec[i-1],n_b_Vec[i-1],Y[i-1],Y_e[i-1],M_nEff[i-1],M_pEff[i-1],tparams)
                                    + pow(r_Vec[i],2.)*n_r_Vec[i]*pow(n_t_Vec[i],2.)*q_nuMUrca(T*T_0/n_t_Vec[i],n_b_Vec[i],Y[i],Y_e[i],M_nEff[i],M_pEff[i],tparams) );
        Q_nuDUrca += 4.*pi*dr*0.5*( pow(r_Vec[i-1],2.)*n_r_Vec[i-1]*pow(n_t_Vec[i-1],2.)*q_nuDUrca(T*T_0/n_t_Vec[i-1],n_b_Vec[i-1],Y[i-1],Y_e[i-1],M_nEff[i-1],M_pEff[i-1],tparams)
                                    + pow(r_Vec[i],2.)*n_r_Vec[i]*pow(n_t_Vec[i],2.)*q_nuDUrca(T*T_0/n_t_Vec[i],n_b_Vec[i],Y[i],Y_e[i],M_nEff[i],M_pEff[i],tparams) );
        Q_nu_NNbrems += 4.*pi*dr*0.5*( pow(r_Vec[i-1],2.)*n_r_Vec[i-1]*pow(n_t_Vec[i-1],2.)*q_nu_NNbrems(T*T_0/n_t_Vec[i-1],n_b_Vec[i-1],Y[i-1],M_nEff[i-1],M_pEff[i-1],tparams)
                                    + pow(r_Vec[i],2.)*n_r_Vec[i]*pow(n_t_Vec[i],2.)*q_nu_NNbrems(T*T_0/n_t_Vec[i],n_b_Vec[i],Y[i],M_nEff[i],M_pEff[i],tparams) );
        if(tparams.SF == true){
            Q_PBF += 4.*pi*dr*0.5*( pow(r_Vec[i-1],2.)*n_r_Vec[i-1]*pow(n_t_Vec[i-1],2.)*q_nuPBF(T*T_0/n_t_Vec[i-1],n_b_Vec[i-1],Y[i-1],M_nEff[i-1],tparams)
                                    + pow(r_Vec[i],2.)*n_r_Vec[i]*pow(n_t_Vec[i],2.)*q_nuPBF(T*T_0/n_t_Vec[i],n_b_Vec[i],Y[i],M_nEff[i],tparams) );
        }
    }

    Q_nu = Q_nuMUrca + Q_nuDUrca + Q_nu_NNbrems + Q_PBF;

    kappa_eCC = kappa_eCore(n_b_Vec.back(),Y.back(), Y_e.back());
    kappa_nCC = kappa_nCore(n_b_Vec.back(), Y.back(), M_nEff.back(), M_pEff.back());

    return;
}

/*
            c_vCore: Core specific heat capacity at fixed volume in reduced units evaluated at T=T_0

    Specific heat capacity from Potekhin, Pons and Page, Space Science Rev. 191, 239 (2015) and Aguilera, Pons and Miralles A&A 486, 255–271 (2008)
    Neutron 3P2 and proton 1S0 pairing gaps in core parameterization taken from model of Kaminker, Haensel and Yakovlev A&A 373, L17 (2001)
        and Andersson, Comer and Glampedakis Nucl. Phys. A 763, 212 (2005)
    Specific pairing gap models:
    Neutron 3P2: Chen et al, Nucl. Phys A 555, 59 (1993) and Elgaroy et al. Nucl. Phys. A 604, 466 (1996), CCDK parametrization in Ho et al. PRC 91, 015806 (2015) Table II.
    Proton 1S0: Takatsuka and Tamagaki Prog. Theor. Phys. 112, 37 (2004), TToa parametrization in Ho et al. PRC 91, 015806 (2015) Table II.
    Proton 1S0 superfluidity control function R_A(T/T_cn) and neutron 3P2 superfluidity control function R_B(T/T_cn) taken from Eq. (18), R_A of Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999)
    Only consider mJ = 0 ("type B") 3P2 superfluidity for neutrons

        Inputs: T: local temperature in K
                n_b: baryon number density in fm^{-3}
                Y: proton fraction of total baryons
                Y_e: electron fraction of total baryons Y_e = n_e/n_b. Equals Y in the absence of muons
                M_nEff: neutron Landau effective mass in MeV
                M_pEff: proton Landau effective mass in MeV
                tparams: TParams object containing information about nucleon superfluid gaps.
        Output: c_ve + c_vn + c_vp: total specific heat capacity at fixed volume in reduced units evaluated at T_0
*/
double c_vCore(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams){

    double p_Fn, p_Fp, p_Fe, p_Fm = 0.;
    double c_ve, c_vm, c_vn, c_vp; //contributions to specific heat capacity from electrons, muons, neutrons, protons.

    double kFp0, kFp2, kFn0, kFn2, Delta_p, Delta_n, T_cp, T_cn, tau_p, tau_n, v_A = 0., v_B = 0.;
    double R_A = 1., R_B = 1.; //superfluid reduction factors for protons and neutrons. Default to 1, are only modified if SF == true

    //parametrization of 1S0 proton pairing gap. Delta0p in MeV, kp0, kp2 in fm^{-1}, kp1, kp3 in fm^{-2}
    double Delta0p = tparams.Delta0_pCore, kp0 = tparams.k0_pCore, kp1 = tparams.k1_pCore, kp2 = tparams.k2_pCore, kp3 = tparams.k3_pCore;
     //parametrization of 3P2 proton pairing gap. Delta0n in MeV, kn0, kn2 in fm^{-1}, kn1, kn3 in fm^{-2}
    double Delta0n = tparams.Delta0_nCore, kn0 = tparams.k0_nCore, kn1 = tparams.k1_nCore, kn2 = tparams.k2_nCore, kn3 = tparams.k3_nCore;

    p_Fn = pow( 3.*pi*pi*(1.-Y)*n_b, 1./3. )*hbarc; //neutron Fermi momentum in MeV
    p_Fp = pow( 3.*pi*pi*Y*n_b, 1./3. )*hbarc; //proton Fermi momentum in MeV

    p_Fe = pow( 3.*pi*pi*Y_e*n_b, 1./3. )*hbarc; //electron chemical potential in MeV
    if( Y_e < Y ) p_Fm = pow( 3.*pi*pi*(Y-Y_e)*n_b, 1./3. )*hbarc; //muon Fermi momentum in MeV if nonzero

    // Specific heat capacity in core evaluated at T=T_0. Scales linearly with T
    c_ve = k_Bcgs*p_Fe*sqrt( p_Fe*p_Fe + M_e*M_e )/(3.*hbarc*hbarc*hbarc)*1e39*k_B*T; //electron specific heat capacity in erg/cm^3/K evaluated at T
    c_vm = k_Bcgs*p_Fm*sqrt( p_Fm*p_Fm + M_m*M_m )/(3.*hbarc*hbarc*hbarc)*1e39*k_B*T; //muon specific heat capacity in erg/cm^3/K evaluated at T
    c_vn = k_Bcgs*p_Fn*sqrt( p_Fn*p_Fn + M_nEff*M_nEff )/(3.*hbarc*hbarc*hbarc)*1e39*k_B*T; //neutron specific heat capacity in erg/cm^3/K evaluated at T
    c_vp = k_Bcgs*p_Fp*sqrt( p_Fp*p_Fp + M_pEff*M_pEff )/(3.*hbarc*hbarc*hbarc)*1e39*k_B*T; //proton specific heat capacity in erg/cm^3/K evaluated at T

    if(tparams.SF == true){
        kFp0 = p_Fp/hbarc - kp0;
        kFp2 = p_Fp/hbarc - kp2;
        if(kFp0 > 0. && kFp2 < 0.){
            Delta_p = Delta0p*kFp0*kFp0/(kFp0*kFp0+kp1)*kFp2*kFp2/(kFp2*kFp2+kp3);
            T_cp = 0.5669*Delta_p; //1S0 proton superfluidity critical temperature in MeV
            tau_p = T*k_B/T_cp; //ratio of T/proton critical temperature
            if(tau_p < 1.){
                v_A = sqrt( 1. - tau_p )*( 1.456 - 0.157/sqrt(tau_p) + 1.764/tau_p ); //ratio of temperature-dependent gap to T, from Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999) Eq. (11)
                R_A = pow( 0.4186 + sqrt( 1.014049 + 0.251001*v_A*v_A ),2.5 )*exp( 1.456 - sqrt( 2.119936 + v_A*v_A ) );
            }
            else R_A = 1.;
        }

        kFn0 = p_Fn/hbarc - kn0;
        kFn2 = p_Fn/hbarc - kn2;
        if(kFn0 > 0. && kFn2 < 0.){
            Delta_n = Delta0n*kFn0*kFn0/(kFn0*kFn0+kn1)*kFn2*kFn2/(kFn2*kFn2+kn3);
            T_cn = 0.1187*Delta_n; //3P2, mJ=0 neutron superfluidity critical temperature in MeV
            tau_n = T*k_B/T_cn; //ratio of T/neutron critical temperature
            if(tau_n < 1.){
                v_B = sqrt( 1. - tau_n )*( 0.7893 + 1.188/tau_n ); //ratio of temperature-dependent gap to T, from Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999) Eq. (11)
                R_B = pow( 0.6893 + sqrt( 0.6241 + 0.07974976*v_B*v_B ),2. )*exp( 1.934 - sqrt( 3.740356 + v_B*v_B ) );
            }
            else R_B = 1.;
        }

    }

    return (c_ve + c_vm + R_A*c_vp + R_B*c_vn )/s_0;
}

/*
            q_nuMUrca: Core Modified Urca neutrino emissivity in reduced units evaluated at T

   Modified Urca neutrino emissivity from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (140-142)
   (also see Schmitt and Shternin arXiv:1711.06520v3).
   Same superfluid pairing gap models described in comment for c_vCore.
   Superfluid reduction factors for neutron MUrca (R^n_AB) and for proton MUrca (R^p_AB) for 1S0 ("type A") paired protons
   and 3P2, mJ=0 ("type B") paired neutrons given by Gusakov A&A 389, 702–715 (2002), Appendix Eq. (A7-A13)
   In-medium corrections from Shternin, Baldo and Haensel Phys. Lett B 786, 28 (2018), Eq. (8-9), (14)

        Inputs: T: local temperature in K
                n_b: baryon number density in fm^{-3}
                Y: proton fraction of total baryons
                Y_e: electron fraction of total baryons Y_e = n_e/n_b. Equals Y in the absence of muons
                M_nEff: neutron Landau effective mass in MeV
                M_pEff: proton Landau effective mass in MeV
                tparams: TParams object containing information about nucleon superfluid gaps.
        Output: q_Murca: modified Urca neutrino emissivity in reduced units evaluated at T
*/
double q_nuMUrca(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams){

    static double hbar = 6.582119569e-22; //hbar in MeV*s
    static double m_pi = 139.57039; //charged pion mass in MeV
    static double G2 = 1.35956e-22; //(reduced Fermi constant times cosine of Cabibbo angle)^2 in MeV^{-4}
    static double gA = 1.26; //axial-vector coupling constant
    static double fpiNN = 1.; //p-wave pi-N coupling constant
    double alpha_n = 1.16, beta_n = 0.68; //Correction and correlation factors (from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), below Eq. (140)). Default values without in-medium corrections
    double RdirA, RdirB; //factors in in-medium correction to alpha_n and beta_n

    double kFp0, kFp2, kFn0, kFn2, Delta_p, Delta_n, T_cp, T_cn, tau_p, tau_n, v_A = 0., v_B = 0.;
    double Rp_AB = 1., Rn_AB = 1.; //superfluid reduction factors for protons and neutrons. Default to 1, are only modified if SF == true

    //parametrization of 1S0 proton pairing gap. Delta0p in MeV, kp0, kp2 in fm^{-1}, kp1, kp3 in fm^{-2}
    double Delta0p = tparams.Delta0_pCore, kp0 = tparams.k0_pCore, kp1 = tparams.k1_pCore, kp2 = tparams.k2_pCore, kp3 = tparams.k3_pCore;
     //parametrization of 3P2 proton pairing gap. Delta0n in MeV, kn0, kn2 in fm^{-1}, kn1, kn3 in fm^{-2}
    double Delta0n = tparams.Delta0_nCore, kn0 = tparams.k0_nCore, kn1 = tparams.k1_nCore, kn2 = tparams.k2_nCore, kn3 = tparams.k3_nCore;

    double p_Fn, p_Fp, p_Fe;
    double mu_e, mu_m;
    double p_Fm = 0., q_MUrca = 0.;

    p_Fn = pow( 3.*pi*pi*(1.-Y)*n_b, 1./3. )*hbarc; //neutron Fermi momentum in MeV
    p_Fp = pow( 3.*pi*pi*Y*n_b, 1./3. )*hbarc; //proton Fermi momentum in MeV
    p_Fe = pow( 3.*pi*pi*Y_e*n_b, 1./3. )*hbarc; //electron Fermi momentum in MeV
    if( Y_e < Y ) p_Fm = pow( 3.*pi*pi*(Y-Y_e)*n_b, 1./3. )*hbarc; //muon Fermi momentum in MeV if nonzero

    if(tparams.SF == true){
        kFp0 = p_Fp/hbarc - kp0;
        kFp2 = p_Fp/hbarc - kp2;
        if(kFp0 > 0. && kFp2 < 0.){
            Delta_p = Delta0p*kFp0*kFp0/(kFp0*kFp0+kp1)*kFp2*kFp2/(kFp2*kFp2+kp3);
            T_cp = 0.5669*Delta_p; //1S0 proton superfluidity critical temperature in MeV
            tau_p = T*k_B/T_cp;
            if(tau_p < 1.) v_A = sqrt( 1. - tau_p )*( 1.456 - 0.157/sqrt(tau_p) + 1.764/tau_p );
        }
        else v_A = 0.;

        kFn0 = p_Fn/hbarc - kn0;
        kFn2 = p_Fn/hbarc - kn2;
        if(kFn0 > 0. && kFn2 < 0.){
            Delta_n = Delta0n*kFn0*kFn0/(kFn0*kFn0+kn1)*kFn2*kFn2/(kFn2*kFn2+kn3);
            T_cn = 0.1187*Delta_n; //3P2, mJ=0 neutron superfluidity critical temperature in MeV
            tau_n = T*k_B/T_cn;
            if(tau_n < 1.) v_B = sqrt( 1. - tau_n )*( 0.7893 + 1.188/tau_n );
        }
        else v_B = 0.;

        //If both v_A and v_B are too small (i.e., if the gap is zero) then set Rn_AB = Rp_AB = 1.
        //Otherwise, do calculation based on Gusakov A&A 389, 702–715 (2002), Appendix Eq. (A7-A13)
        if(v_A < 1e-15 && v_B < 1e-15){
            Rn_AB = 1.;
            Rp_AB = 1.;
        }
        else{
            if( sqrt(v_A*v_A+v_B*v_B) < 5. ){
                Rn_AB = Rn_AB_RegionIV(v_B,v_A);
                Rp_AB = Rp_AB_RegionIV(v_B,v_A);
            }
            else if( v_B > v_A ){
                Rn_AB = Rn_AB_RegionI(v_B,v_A);
                Rp_AB = Rp_AB_RegionI(v_B,v_A);
            }
            else if( v_B > v_A/3. ){
                Rn_AB = Rn_AB_RegionII(v_B,v_A);
                Rp_AB = Rp_AB_RegionII(v_B,v_A);
            }
            else{
                Rn_AB = Rn_AB_RegionIII(v_B,v_A);
                Rp_AB = Rp_AB_RegionIII(v_B,v_A);
            }
        }
    }

    //Modified Urca neutrino emissivity in erg/s/cm^3. Scales as T^8
    //Define to be negative so include as + q_MUrca in the generalized heat equation
    //Electronic MUrca process
    if( tparams.InMediumMUrca == true ){ //compute in-medium correction factors for electrons if enabled
        mu_e = sqrt( p_Fe*p_Fe + M_e*M_e );
//        RdirA = 2.*M_pEff*M_pEff*mu_e*mu_e/p_Fe*( p_Fe*( p_Fn*p_Fn + p_Fp*p_Fp - p_Fe*p_Fe )/( p_Fp*p_Fp*( p_Fn*p_Fn - ( p_Fp + p_Fe )*( p_Fp + p_Fe ) )*( p_Fn*p_Fn - ( p_Fp - p_Fe )*( p_Fp - p_Fe ) ) )
//                                                  + ( atanh( p_Fp/( p_Fn + p_Fe ) ) - atanh( p_Fp/( p_Fn - p_Fe ) ) )/( 2.*p_Fp*p_Fp*p_Fp ) );
        RdirA = 2.*M_pEff*M_pEff*mu_e*mu_e/p_Fe*( p_Fe*( p_Fn*p_Fn + p_Fp*p_Fp - p_Fe*p_Fe )/( p_Fp*p_Fp*( ( p_Fe - p_Fn - p_Fp )*( p_Fe + p_Fn - p_Fp )*( p_Fe - p_Fn + p_Fp )*( p_Fe + p_Fn + p_Fp ) ) )
                                                  + ( atanh( p_Fp/( p_Fn + p_Fe ) ) - atanh( p_Fp/( p_Fn - p_Fe ) ) )/( 2.*p_Fp*p_Fp*p_Fp ) );
        RdirB = 2.*M_nEff*M_nEff*mu_e*mu_e/(p_Fe*p_Fp)*( 2.*p_Fe*p_Fp/( p_Fe*p_Fe*p_Fe*p_Fe + ( p_Fn*p_Fn - p_Fp*p_Fp )*( p_Fn*p_Fn - p_Fp*p_Fp ) - 2.*p_Fe*p_Fe*( p_Fn*p_Fn + p_Fp*p_Fp ) ) );

        alpha_n = 0.6*RdirA + 0.2*RdirB;
        beta_n = 1.; //include all in-medium corrections in alpha_n
    }
    if( p_Fe + 3.*p_Fp - p_Fn > 0. ){ //if proton branch is allowed
        q_MUrca += -11513./(30240.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi,4. )*pow( k_B*T,8. )*alpha_n*beta_n/pow( hbarc,3. )*MeVtoErg*1e39/hbar*(
                    Rn_AB*M_nEff*M_nEff*M_nEff*M_pEff*p_Fp + Rp_AB*M_nEff*M_pEff*M_pEff*M_pEff*pow( p_Fe + 3.*p_Fp - p_Fn,2. )/(8.*p_Fe) );
    }
    else{ //if proton branch is forbidden
        q_MUrca += -11513./(30240.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi,4. )*pow( k_B*T,8. )*alpha_n*beta_n/pow( hbarc,3. )*MeVtoErg*1e39/hbar*Rn_AB*M_nEff*M_nEff*M_nEff*M_pEff*p_Fp;
    }

    //Muonic MUrca process
    if( Y_e < Y ){
        if( tparams.InMediumMUrca == true ){ //compute in-medium correction factors for muons if enabled
            mu_m = sqrt( p_Fm*p_Fm + M_m*M_m );
            RdirA = 2.*M_pEff*M_pEff*mu_m*mu_m/p_Fm*( p_Fm*( p_Fn*p_Fn + p_Fp*p_Fp - p_Fm*p_Fm )/( p_Fp*p_Fp*( ( p_Fm - p_Fn - p_Fp )*( p_Fm + p_Fn - p_Fp )*( p_Fm - p_Fn + p_Fp )*( p_Fm + p_Fn + p_Fp ) ) )
                                                      + ( atanh( p_Fp/( p_Fn + p_Fm ) ) - atanh( p_Fp/( p_Fn - p_Fm ) ) )/( 2.*p_Fp*p_Fp*p_Fp ) );
            RdirB = 2.*M_nEff*M_nEff*mu_m*mu_m/(p_Fm*p_Fp)*( 2.*p_Fm*p_Fp/( p_Fm*p_Fm*p_Fm*p_Fm + ( p_Fn*p_Fn - p_Fp*p_Fp )*( p_Fn*p_Fn - p_Fp*p_Fp ) - 2.*p_Fm*p_Fm*( p_Fn*p_Fn + p_Fp*p_Fp ) ) );

            alpha_n = 0.6*RdirA + 0.2*RdirB;
            beta_n = 1.; //include all in-medium corrections in alpha_n
        }
        if( p_Fm + 3.*p_Fp - p_Fn > 0. ){ //if proton branch is allowed
            q_MUrca += -11513./(30240.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi,4. )*pow( k_B*T,8. )*alpha_n*beta_n/pow( hbarc,3. )*MeVtoErg*1e39/hbar*(
                        Rn_AB*M_nEff*M_nEff*M_nEff*M_pEff*p_Fp + Rp_AB*M_nEff*M_pEff*M_pEff*M_pEff*pow( p_Fm + 3.*p_Fp - p_Fn,2. )/(8.*p_Fm) )*pow( Y/Y_e - 1.,1./3. );
        }
        else{ //if proton branch is forbidden
            q_MUrca += -11513./(30240.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi,4. )*pow( k_B*T,8. )*alpha_n*beta_n/pow( hbarc,3. )*MeVtoErg*1e39/hbar*Rn_AB*M_nEff*M_nEff*M_nEff*M_pEff*p_Fp*pow( Y/Y_e - 1.,1./3. );
        }
    }

    return q_MUrca*t_0/(T_0*s_0); //converts to reduced units.
}

/*
            q_nuMUrca: Core Direct Urca neutrino emissivity in reduced units evaluated at T

   Direct Urca neutrino emissivity from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001)
   (also see Schmitt and Shternin arXiv:1711.06520v3). Assume no muons in core
   Same superfluid pairing gap models described in comment for c_vCore.
   Superfluid reduction factor taken from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001) Eq. (199-206).
   We use approximation that R_D = min(R_AD, R_BD) (Eq. (204)) where R_AD is the reduction factor for 1S0 paired ("type A") protons
   and R_BD is the reduction factor for 3P2, mJ=0 ("type B") paired neutrons

        Inputs: T: local temperature in K
                n_b: baryon number density in fm^{-3}
                Y: proton fraction of total baryons
                Y_e: electron fraction of total baryons Y_e = n_e/n_b. Equals Y in the absence of muons
                M_nEff: neutron Landau effective mass in MeV
                M_pEff: proton Landau effective mass in MeV
                tparams: TParams object containing information about nucleon superfluid gaps.
        Output: q_Durca: direct Urca neutrino emissivity in reduced units evaluated at T
*/
double q_nuDUrca(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams){

    static double hbar = 6.582119569e-22; //hbar in MeV*s
    static double G2 = 1.35956e-22; //(reduced Fermi constant times cosine of Cabibbo angle)^2 in MeV^{-4}
    static double gA = 1.26; //axial-vector coupling constant

    double kFp0, kFp2, kFn0, kFn2, Delta_p, Delta_n, T_cp, T_cn, tau_p, tau_n, v_A = 0., v_B = 0.;
    double R_D = 1.; //superfluid reduction factor for direct Urca process. Defaults to 1, is only modified if SF == true
    double R_AD, R_BD; //superfluid reduction factors for type-A and type-B pairing, used for proton pairing alone and neutron pairing alone

    //parametrization of 1S0 proton pairing gap. Delta0p in MeV, kp0, kp2 in fm^{-1}, kp1, kp3 in fm^{-2}
    double Delta0p = tparams.Delta0_pCore, kp0 = tparams.k0_pCore, kp1 = tparams.k1_pCore, kp2 = tparams.k2_pCore, kp3 = tparams.k3_pCore;
     //parametrization of 3P2 proton pairing gap. Delta0n in MeV, kn0, kn2 in fm^{-1}, kn1, kn3 in fm^{-2}
    double Delta0n = tparams.Delta0_nCore, kn0 = tparams.k0_nCore, kn1 = tparams.k1_nCore, kn2 = tparams.k2_nCore, kn3 = tparams.k3_nCore;

    double p_Fn, p_Fp, mu_e, p_Fe, p_Fm = 0., q_DUrca;

    p_Fn = pow( 3.*pi*pi*(1.-Y)*n_b, 1./3. )*hbarc; //neutron Fermi momentum in MeV
    p_Fp = pow( 3.*pi*pi*Y*n_b, 1./3. )*hbarc; //proton Fermi momentum in MeV
    p_Fe = pow( 3.*pi*pi*Y_e*n_b, 1./3. )*hbarc; //electron Fermi momentum in MeV
    mu_e = sqrt( p_Fe*p_Fe + M_e*M_e );
    if( Y_e < Y ) p_Fm = pow( 3.*pi*pi*(Y-Y_e)*n_b, 1./3. )*hbarc; //muon Fermi momentum in MeV if nonzero

    //Direct Urca neutrino emissivity for both electronic and muonic processes in erg/s/cm^3. Zero if Y is below direct Urca threshold. Scales as T^6.
    //Define to be negative so include as + q_DUrca in the generalized heat equation
    if( p_Fp + p_Fe > p_Fn && p_Fp + p_Fm > p_Fn ){
        q_DUrca = 2.*( -457.*pi/10080.*G2*(1.+3.*gA*gA)*M_nEff*M_pEff*mu_e*pow( k_B*T,6. )/pow(hbarc,3.)*MeVtoErg*1e39/hbar );
    }
    else if( p_Fp + p_Fe > p_Fn ){ //if only electronic muonic processes allowed. Can never have muonic but not electronic
        q_DUrca = -457.*pi/10080.*G2*(1.+3.*gA*gA)*M_nEff*M_pEff*mu_e*pow( k_B*T,6. )/pow(hbarc,3.)*MeVtoErg*1e39/hbar;
    }
    else q_DUrca = 0.;

    if(tparams.SF == true){
        kFp0 = p_Fp/hbarc - kp0;
        kFp2 = p_Fp/hbarc - kp2;
        if(kFp0 > 0. && kFp2 < 0.){
            Delta_p = Delta0p*kFp0*kFp0/(kFp0*kFp0+kp1)*kFp2*kFp2/(kFp2*kFp2+kp3);
            T_cp = 0.5669*Delta_p; //1S0 proton superfluidity critical temperature in MeV
            tau_p = T*k_B/T_cp;
            if(tau_p < 1.){
                v_A = sqrt( 1. - tau_p )*( 1.456 - 0.157/sqrt(tau_p) + 1.764/tau_p );
                R_AD = pow( 0.2312 + sqrt( 0.59105344 + 0.02067844*v_A*v_A ),5.5 )*exp( 3.427 - sqrt( 11.744329 + v_A*v_A ) );
            }
            else R_AD = 1.;
        }
        else R_AD = 1.;

        kFn0 = p_Fn/hbarc - kn0;
        kFn2 = p_Fn/hbarc - kn2;
        if(kFn0 > 0. && kFn2 < 0.){
            Delta_n = Delta0n*kFn0*kFn0/(kFn0*kFn0+kn1)*kFn2*kFn2/(kFn2*kFn2+kn3);
            T_cn = 0.1187*Delta_n; //3P2, mJ=0 neutron superfluidity critical temperature in MeV
            tau_n = T*k_B/T_cn;
            if(tau_n < 1.){
                v_B = sqrt( 1. - tau_n )*( 0.7893 + 1.188/tau_n );
                R_BD = pow( 0.2546 + sqrt( 0.55562116 + 0.01648656*v_B*v_B ),5. )*exp( 2.701 - sqrt( 7.295401 + v_B*v_B ) );
            }
            else R_BD = 1.;
        }
        else R_BD = 1.;

        R_D = std::min(R_AD,R_BD);
    }

    return R_D*q_DUrca*t_0/(T_0*s_0); //converts to reduced units.
}

/*
            q_nu_nnbrems: Core nucleon-nucleon bremsstrahlung neutrino emissivity in reduced units evaluated at T

   Nucleon-nucleon bremsstrahlung neutrino emissivities from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (165-167).
   Same superfluid pairing gap models described in comment for c_vCore.
   Superfluid reduction factors for 3P2, mJ=0 ("type B") paired neutrons and 1S0 ("type A") paired protons given by:
   nn: R^nn_{nB} = R^{pp}_{pA}(v_n) from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (221) and (228) (orig. Yakovlev and Levenfish, A&A 29, 717 (1995) Eq. (59)).
   pp: R^pp = R^pp_{pA}(v_p) from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (221)
   np: R^{np}_{BA} = from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (199,205,229)

        Inputs: T: local temperature in K
                n_b: baryon number density in fm^{-3}
                Y: proton fraction of total baryons
                M_nEff: neutron Landau effective mass in MeV
                M_pEff: proton Landau effective mass in MeV
                tparams: TParams object containing information about nucleon superfluid gaps.
        Output: q_NNbrems: neutron-neutron bremsstrahlung neutrino emissivity in reduced units evaluated at T
*/
double q_nu_NNbrems(double T, double n_b, double Y, double M_nEff, double M_pEff, const TParams & tparams){

    static double hbar = 6.582119569e-22; //hbar in MeV*s
    static double m_pi0 = 134.9768; //neutral pion mass in MeV
    static double G2 = 1.35956e-22; //(reduced Fermi constant times cosine of Cabibbo angle)^2 in MeV^{-4}
    static double gA = 1.26; //axial-vector coupling constant
    static double fpiNN = 1.; //p-wave pi-N coupling constant
    //Correction factors (from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), below Eq. (167))
    static double alpha_nn = 0.59, beta_nn = 0.56;
    static double alpha_np = 1.06, beta_np = 0.66;
    static double alpha_pp = 0.11, beta_pp = 0.7;

    double kFp0, kFp2, kFn0, kFn2, Delta_p, Delta_n, T_cp, T_cn, tau_p, tau_n, v_A = 0., v_B = 0.;
    double Rnn_nB = 1., Rnp_BA = 1., Rpp_pA = 1.; //superfluid reduction factors for neutron-neutron, proton-neutron and proton-proton bremsstrahlung. Default to 1, are only modified if SF == true

    //parametrization of 1S0 proton pairing gap. Delta0p in MeV, kp0, kp2 in fm^{-1}, kp1, kp3 in fm^{-2}
    double Delta0p = tparams.Delta0_pCore, kp0 = tparams.k0_pCore, kp1 = tparams.k1_pCore, kp2 = tparams.k2_pCore, kp3 = tparams.k3_pCore;
     //parametrization of 3P2 proton pairing gap. Delta0n in MeV, kn0, kn2 in fm^{-1}, kn1, kn3 in fm^{-2}
    double Delta0n = tparams.Delta0_nCore, kn0 = tparams.k0_nCore, kn1 = tparams.k1_nCore, kn2 = tparams.k2_nCore, kn3 = tparams.k3_nCore;

    double q_nnbrems, q_npbrems, q_ppbrems; //nn, np, pp bremsstrahlung contributions to neutrino emissivity
    double N_v = 3.; //number of neutrino flavours (all three can be emitted in PBF)

    double p_Fn, p_Fp, M_nEff2, M_pEff2;
    p_Fn = pow( 3.*pi*pi*(1.-Y)*n_b, 1./3. )*hbarc; //neutron Fermi momentum in MeV
    p_Fp = pow( 3.*pi*pi*Y*n_b, 1./3. )*hbarc; //proton Fermi momentum in MeV

    M_nEff2 = M_nEff*M_nEff;
    M_pEff2 = M_pEff*M_pEff;

    if(tparams.SF == true){
        kFp0 = p_Fp/hbarc - kp0;
        kFp2 = p_Fp/hbarc - kp2;
        if(kFp0 > 0. && kFp2 < 0.){
            Delta_p = Delta0p*kFp0*kFp0/(kFp0*kFp0+kp1)*kFp2*kFp2/(kFp2*kFp2+kp3);
            T_cp = 0.5669*Delta_p; //1S0 proton superfluidity critical temperature in MeV
            tau_p = T*k_B/T_cp;
            if(tau_p < 1.) v_A = sqrt( 1. - tau_p )*( 1.456 - 0.157/sqrt(tau_p) + 1.764/tau_p ); //v_2 in Gusakov 2002's notation
        }
        else v_A = 0.;

        kFn0 = p_Fn/hbarc - kn0;
        kFn2 = p_Fn/hbarc - kn2;
        if(kFn0 > 0. && kFn2 < 0.){
            Delta_n = Delta0n*kFn0*kFn0/(kFn0*kFn0+kn1)*kFn2*kFn2/(kFn2*kFn2+kn3);
            T_cn = 0.1187*Delta_n; //3P2, mJ=0 neutron superfluidity critical temperature in MeV
            tau_n = T*k_B/T_cn;
            if(tau_n < 1.) v_B = sqrt( 1. - tau_n )*( 0.7893 + 1.188/tau_n ); //v_1 in Gusakov 2002's notation
        }
        else v_B = 0.;

        //If both v_A and v_B are too small (i.e., if the gap is zero) then set Rn_AB = Rp_AB = 1.
        //Otherwise, do calculation based on Gusakov A&A 389, 702–715 (2002), Appendix Eq. (A7-A13)
        if(v_A < 1e-15 && v_B > 1e-15){ //superfluid neutrons but not protons
            Rnn_nB = Rpp_pACalc(v_B);
            Rnp_BA = RD_BACalc(v_B,0.);
        }
        else if(v_A > 1e-15 && v_B < 1e-15){ //superconducting protons but not neutrons
            if (RD_ACalc(v_A) > 1e-15) {
                Rnp_BA = RD_BACalc(0,v_A)/RD_ACalc(v_A)*Rnp_pACalc(v_A);
            }
            else {
                Rnp_BA = 0.;
            }
            Rpp_pA = Rpp_pACalc(v_A);
        }
        else if(v_B > 1e-15 && v_A > 1e-15){ //both superfluid
            Rnn_nB = Rpp_pACalc(v_B);
            if (RD_ACalc(v_A) > 1e-15) {
                Rnp_BA = RD_BACalc(v_B,v_A)/RD_ACalc(v_A)*Rnp_pACalc(v_A);
            }
            else {
                Rnp_BA = 0.;
            }
            Rpp_pA = Rpp_pACalc(v_A);
        }
    }

    //Neutron-neutron bremsstrahlung neutrino emissivity in erg/s/cm^3. Scales as T^8
    //Define to be negative so include as + q_NNbrems in the generalized heat equation
    q_nnbrems = -41./(14175.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi0,4. )*M_nEff2*M_nEff2*p_Fn*pow( k_B*T,8. )*alpha_nn*beta_nn*N_v/pow(hbarc,3.)*MeVtoErg*1e39/hbar*Rnn_nB;
    q_npbrems = -82./(14175.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi0,4. )*M_nEff2*M_pEff2*p_Fp*pow( k_B*T,8. )*alpha_np*beta_np*N_v/pow(hbarc,3.)*MeVtoErg*1e39/hbar*Rnp_BA;
    q_ppbrems = -41./(14175.*2.*pi)*G2*gA*gA*pow( fpiNN/m_pi0,4. )*M_pEff2*M_pEff2*p_Fp*pow( k_B*T,8. )*alpha_pp*beta_pp*N_v/pow(hbarc,3.)*MeVtoErg*1e39/hbar*Rpp_pA;

    return ( q_nnbrems + q_npbrems + q_ppbrems )*t_0/(T_0*s_0); //converts to reduced units.
}

/*
            q_nuPBF: Core superfluid neutron pairing breaking and formation (PBF) neutrino emissivity in reduced units evaluated at T

   PBF Urca neutrino emissivity from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (236-241),
   with anomalous axial PBF contribution switched off (i.e., PBF emission is reduced) according to Leinson PRC 81, 025501 (2010).
   Same superfluid pairing gap models described in comment for c_vCore.

        Inputs: T: local temperature in K
                n_b: baryon number density in fm^{-3}
                Y: proton fraction of total baryons
                M_nEff: neutron Landau effective mass in MeV
                tparams: TParams object containing information about nucleon superfluid gaps.
        Output: q_PBF: modified Urca neutrino emissivity in reduced units evaluated at T
*/
double q_nuPBF(double T, double n_b, double Y, double M_nEff, const TParams & tparams){

    double kFn0, kFn2, Delta_n, T_cn, tau_n, v = 0., F_B;
    //parametrization of 3P2 proton pairing gap. Delta0n in MeV, kn0, kn2 in fm^{-1}, kn1, kn3 in fm^{-2}
    double Delta0n = tparams.Delta0_nCore, kn0 = tparams.k0_nCore, kn1 = tparams.k1_nCore, kn2 = tparams.k2_nCore, kn3 = tparams.k3_nCore;

    double p_Fn = pow( 3.*pi*pi*(1.-Y)*n_b, 1./3. )*hbarc; //neutron Fermi momentum in MeV
    double q_PBF;

    double a_nB = 4.17; //numerical factor for 3P2, mJ=0 ("type B") neutron superfluidity
    double N_v = 3.; //number of neutrino flavours (all three can be emitted in PBF)
    double CollectiveEffectSuppression = 0.19; //Suppression factor of PBF in "type B" neutron superfluidity from Eq. (97) of Leinson PRC 81, 025501 (2010).

    kFn0 = p_Fn/hbarc - kn0;
    kFn2 = p_Fn/hbarc - kn2;
    if(kFn0 > 0. && kFn2 < 0.){
        Delta_n = Delta0n*kFn0*kFn0/(kFn0*kFn0+kn1)*kFn2*kFn2/(kFn2*kFn2+kn3);
        T_cn = 0.1187*Delta_n; //3P2, mJ=0 neutron superfluidity critical temperature in MeV
        tau_n = T*k_B/T_cn;
        if(tau_n < 1.) v = sqrt( 1. - tau_n )*( 0.7893 + 1.188/tau_n ); //if tau_n > 1, v = 0 (set in variable declaration)
    }

    //Expression for F_B(v) from Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001), Eq. (241)
    F_B = ( 1.204*v*v + 3.733*v*v*v*v + 0.3191*v*v*v*v*v*v )/( 1. + 0.3511*v*v )*pow( 0.7591 + sqrt( 0.05803281 + 0.3145*v*v ),2. )*exp( 0.4616 - sqrt( 0.21307456 + 4.*v*v ) );

    //3P2 paired neutron PBF neutrino emissivity in erg/s/cm^3. Scales as T^7
    //Define to be negative so include as + q_PBF in the generalized heat equation
    q_PBF = -1.170e21*( M_nEff/M_n )*( p_Fn/M_n )*pow( T/1e9,7. )*a_nB*F_B*N_v*CollectiveEffectSuppression;

    return q_PBF*t_0/(T_0*s_0); //converts to reduced units.
}

/*
        Generate interpolating functions for core thermal parameters C_v (volume-integrated specific heat capacity) and Q_nu (volume-
        integrated neutrino emissivity). Both are in reduced units.
        Inputs: TVec: vector of uniform redshifted core temperature in reduced units
                C_vVec, Q_nuVec: vectors of volume-integrated specific heat capacity and neutrino emissivity of the core at the uniform redshifted temperatures in TVec. In reduced units.
                corethermalparams: instance of class CoreThermalParams used to store interpolating functions for core thermal parameters and heat flux prefactor
*/
void Compute_CoreThermalParams_Interpolators(std::vector<double> & TVec, std::vector<double> & C_vVec, std::vector<double> & Q_nuVec, CoreThermalParams & corethermalparams){

	//converts vectors of data to arrays to use in interpolation
	double* T_Interp = &TVec[0];
	double* C_v_Interp = &C_vVec[0];
	double* Q_nu_Interp = &Q_nuVec[0];

	size_t N_dat = TVec.size(); //number of entries in TVec, C_vVec and Q_nuVec

	//Generates spline and acc objects for each quantity to interpolate them and assigns them to EOSInterpolation object
	gsl_interp_accel *C_vCore_acc = gsl_interp_accel_alloc ();
	gsl_spline *C_vCore = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(C_vCore, T_Interp, C_v_Interp, N_dat);

	gsl_interp_accel *Q_nuCore_acc = gsl_interp_accel_alloc ();
	gsl_spline *Q_nuCore = gsl_spline_alloc (gsl_interp_steffen, N_dat);
	gsl_spline_init(Q_nuCore, T_Interp, Q_nu_Interp, N_dat);

    //Stores interpolation functions in instance of class "CoreThermalParams"
    //Also creates accelerator objects for each interpolator and stores them in this class instance.
    corethermalparams.C_vCore = C_vCore;
    corethermalparams.C_vCore_acc = C_vCore_acc;
    corethermalparams.Q_nuCore = Q_nuCore;
    corethermalparams.Q_nuCore_acc = Q_nuCore_acc;

    corethermalparams.HeatFlux_pref = 0.; //initialize heat flux prefactor to zero initially

    return;
}

/*
        Functions to compute modified Urca reduction functions using fits from Gusakov A&A 389, 702–715 (2002), Appendix Eq. (A7-A13).
        Rn_AB_Region1-IV functions are for neutron branch MUrca, Rp_AB_RegionI-IV are for proton branch MUrca
        Note that v_1 = v_n, v_2 = v_p in Gusakov's notation
        Inputs: v_n, v_p: ratios of temperature-dependent gap to temperature (dimensionless) for neutrons and protons. Computed using
                approximations in Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999), Eq. (11)
*/
double Rn_AB_RegionI(double v_n, double v_p){

    static double p1 = -0.719681, p2 = -0.024591, p3 = 0.297357, p4 = 1.260056, p5 = 0.100466;
    static double p6 = 0.148464, p7 = 0.253881, p8 = 140.3699, p9 = 0.132615, p10 = 0.280765;
    static double p11 = 0.375796, p12 = -0.096843, p13 = 3.100942, p14 = 0.275434, p15 = 0.330574;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double phi;
    if(v_p < 1e-15 ) phi = 0.5*pi;
    else if(v_n < 1e-15) phi = 0.;
    else phi = atan(v_n/v_p);
    double t = v_p*v_p/(v*v); //cos(phi)^2 where phi = atan(v_n/v_p)
    double y = sin( phi + p15 )*sin( phi + p15 );
    double z = cos( phi + p14 )*cos( phi + p14 );

    double A = p1 + p2*t*t + p4/( 1. + p3*t ) - p5*t;
    double B = p6 + p7*t*t + p9*t/( 1. + p8*t*t ) - p10*t;
    double C = p11 - p12/( y*pow( 1. + p13*z*z,3. ) );

    return exp( -A*v*v/pow( 1. + B*v*v,C ) );
}

double Rn_AB_RegionII(double v_n, double v_p){

    static double p1 = -6.475443, p2 = -1.186294, p3 = 0.591347, p4 = 6.953996, p5 = 3.366945;
    static double p6 = -9.172994, p7 = -2.675793, p8 = 1.053679, p9 = 10.38526, p10 = 7.138369;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double z = v_p*v_p/(v*v); //cos(phi)^2 where phi = atan(v_n/v_p)

    double A = p1 + p2*z*z + p4/( 1. + p3*z ) + p5*z;
    double B = 0.035;
    double C = p6 + p7*z*z + p9/( 1. + p8*z ) + p10*z;

    return exp( -A*v*v/pow( 1. + B*v*v,C ) );
}

double Rn_AB_RegionIII(double v_n, double v_p){

    static double p1 = 0.316041, p2 = -289.2964, p3 = 2480.961, p4 = -268.8219, p5 = 1984.115;
    static double p6 = 3503.094, p7 = 0.331551, p8 = -0.265977, p9 = 1098.324, p10 = 65528.01;
    static double p11 = 0.024500, p12 = 0.120536, p13 = 89.79866, p14 = 5719.134, p15 = 285.8473;
    static double p16 = 0.402111, p17 = 16657.19;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double phi;
    if(v_p < 1e-15 ) phi = 0.5*pi;
    else if(v_n < 1e-15) phi = 0.;
    else phi = atan(v_n/v_p);
    double phi2 = phi*phi;
    double phi3 = phi2*phi;
    double phi4 = phi3*phi;

    double A = p12*( 1. + p15*phi2 + p16*phi3 + p17*phi4 )/( 1. + p13*phi2 + p14*phi3 );
    double B = p11 + p7/( 1. + p8*phi2 + p9*phi3 + p10*phi4 );
    double C = p1*( 1. + p4*phi2 + p5*phi3 + p6*phi4 )/( 1. + p2*phi2 + p3*phi3 );

    return exp( -A*v*v/pow( 1. + B*v*v,C ) );
}

double Rn_AB_RegionIV(double v_n, double v_p){

    static double p1 = 0.565001, p2 = 0.087929, p3 = 0.006756, p4 = 1.667194e-4, p5 = 3.782805e-6;
    static double p6 = 0.173165, p7 = 1.769413e-5, p8 = 7.710124e-8, p9 = 0.001695;

    double v_n2 = v_n*v_n;
    double v_p2 = v_p*v_p;
    double v_n4 = v_n2*v_n2;
    double v_n6 = v_n4*v_n2;
    double v_n8 = v_n4*v_n4;

    double A = p1*v_n2 + p2*v_p2 + p3*v_n2*v_p2+ p4*v_n6;
    double B = sqrt( 1. + p6*v_n2 + p7*v_p2 + p5*v_n8 + p8*v_n6 );
    double C = 1. + p9*v_p2*v_p2;

    return C*exp( -A/B );
}

double Rp_AB_RegionI(double v_n, double v_p){

    static double p1 = 0.288203, p2 = -0.124974, p3 = 17.39273, p4 = 0.083392, p5 = 0.059046;
    static double p6 = 0.028084, p7 = -0.019990, p8 = 28.37210, p9 = 0.244471, p10 = -0.610470;
    static double p11 = 0.023288, p12 = 0.475196, p13 = -0.180420, p14 = 25.51325, p15 = 0.281721;
    static double p16 = -0.080480, p17 = -0.191637;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double phi;
    if(v_p < 1e-15 ) phi = 0.5*pi;
    else if(v_n < 1e-15) phi = 0.;
    else phi = atan(v_n/v_p);
    double t = v_p*v_p/(v*v); //cos(phi)^2 where phi = atan(v_n/v_p)
    double y = sin( phi + p15 )*sin( phi + p15 );

    double A = p1 + p2*phi + p4*phi/pow( 1. + p3*t*phi,2. ) + p5*t*phi*phi;
    double B = p6 + p7*phi + p9*phi/pow( 1. + p8*t*phi + p10*y*t,2. ) + p11*phi*phi;
    double C = p12 - p13*t + p16/pow( 1. + p14*t*t,2. ) + p17*t*phi;

    return exp( -A*v*v/pow( 1. + B*v*v,C ) );
}

double Rp_AB_RegionII(double v_n, double v_p){

    static double p1 = 0.398261, p2 = -0.054952, p3 = -0.084964, p4 = -0.036240, p5 = -0.168712;
    static double p6 = -0.704750, p7 = -0.066981, p8 = 1.223731, p9 = 0.363094, p10 = -0.357641;
    static double p11 = 0.869196, p12 = -0.364248, p13 = 2.668230, p14 = -0.765093, p15 = -4.198753;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double phi;
    if(v_p < 1e-15 ) phi = 0.5*pi;
    else if(v_n < 1e-15) phi = 0.;
    else phi = atan(v_n/v_p);
    double z = v_p*v_p/(v*v); //cos(phi)^2 where phi = atan(v_n/v_p)

    double A = p1 + p2*z*z + p4/( 1. + p3*z ) + p5*phi;
    double B = p6 + p7*z*z + p8/( 1. + p9*z ) + p10*phi;
    double C = p11 - p12*phi*phi + p13*phi*phi/( 1. + p14*z ) + p15*phi;

    //This function has a problem with B becoming negative, probably due to a lack of digits in Table 3
    //of Gusakov A&A 389, 702–715 (2002). So we set the output equal to zero in this case since Lim x-> 0 exp( -A*v*v/x )->0
    double output;
    if( 1. + B*v*v > 0.){
        output = exp( -A*v*v/pow( 1. + B*v*v,C ) );
    }
    else output = 0.;

    return output;
}

double Rp_AB_RegionIII(double v_n, double v_p){

    static double p1 = 0.387542, p2 = -195.5462, p3 = 3032.985, p4 = -189.0452, p5 = 3052.617;
    static double p6 = 442.6031, p7 = 0.041901, p8 = -0.022201, p9 = 5608.168, p10 = -10761.76;
    static double p11 = 0.064643, p12 = 0.296253, p13 = 106.3387, p14 = -75.36126, p15 = 84.65801;
    static double p16 = 0.530223, p17 = -86.76801;

    double v = sqrt( v_n*v_n + v_p*v_p );
    double phi;
    if(v_p < 1e-15) phi = 0.5*pi;
    else if(v_n < 1e-15) phi = 0.;
    else phi = atan(v_n/v_p);
    double phi2 = phi*phi;
    double phi3 = phi2*phi;
    double phi4 = phi3*phi;

    double A = p12*( 1. + p15*phi2 + p16*phi3 + p17*phi4 )/( 1. + p13*phi2 + p14*phi3 );
    double B = p11 + p7/( 1. + p8*phi2 + p9*phi3 + p10*phi4 );
    double C = p1*( 1. + p4*phi2 + p5*phi3 + p6*phi4 )/( 1. + p2*phi2 + p3*phi3 );

    //This function has a problem with B becoming negative, probably due to a lack of digits in Table 3
    //of Gusakov A&A 389, 702–715 (2002). So we set the output equal to zero in this case since Lim x-> 0 exp( -A*v*v/x )->0
    double output;
    if( 1. + B*v*v > 0.){
        output = exp( -A*v*v/pow( 1. + B*v*v,C ) );
    }
    else output = 0.;

    return output;
}

double Rp_AB_RegionIV(double v_n, double v_p){

    static double p1 = 0.272730, p2 = 0.165858, p3 = 0.005903, p4 = 2.555386e-5, p5 = 2.593057e-7;
    static double p6 = 0.023930, p7 = 0.006180, p8 = 1.289532e-5, p9 = 0.005368;

    double v_p2 = v_p*v_p;
    double v_p4 = v_p2*v_p2;
    double v_p6 = v_p4*v_p2;
    double v_n2 = v_n*v_n;

    double A = p1*v_p2 + p2*v_n2 + p3*v_p2*v_n2 + p4*v_p6 + p5*v_n2*v_n2*v_n2;
    double B = 1. + p6*v_p2 + p7*v_n2 + p8*v_p4;
    double C = 1. + p9*v_n2*v_n2;

    return C*exp( -A/B );
}

/*
        Functions to compute proton-proton, proton-neutron and neutron-neutron bremsstrahlung neutrino emissivity superfluid reduction factors.
        From Yakovlev, Kaminker, Gnedin and Haensel, Physics Reports 354, 1 (2001):
        RD_BA: Eq. (205) (direct Urca reduction factor)
        RD_A: Eq. (199) (direct Urca reduction factor)
        Rnp_pA: Eq. (220)
        Rpp_pA: Eq. (221)
        Inputs: v: ratio of temperature-dependent gap to temperature (dimensionless). Takes either v_n = v_1 or v_p = v_2. Computed using
                approximations in Yakovlev, Levenfish and Shibanov, Phys. Usp. 42, 737 (1999), Eq. (11)
                or both v_n, v_p
*/
double RD_ACalc(double v_p){

    return pow( 0.2312 + sqrt( 0.59105344 + 0.02067844*v_p*v_p ),5.5 )*exp( 3.427 - sqrt( 11.744329 + v_p*v_p ) );
}

double RD_BACalc(double v_n, double v_p){

    double v_n2 = v_n*v_n;
    double v_p2 = v_p*v_p;

    double numer = 1e4 - 2.839*v_p2*v_p2 - 5.022*v_n2*v_n2;
    double denom = 1e4 + 757.0*v_p2 + 1494.*v_n2 + 211.1*v_n2*v_p2 + 0.4832*v_n2*v_n2*v_p2*v_p2;

    return numer/denom;
}

double Rnp_pACalc(double v){

    double a = 0.9982 + sqrt( 3.24e-6 + 0.14554225*v*v );
    double b = 0.3949 + sqrt( 0.36614601 + 0.07107556*v*v );

    return 1./2.732*( a*exp( 1.306 - sqrt( 1.705636 + v*v ) ) + 1.732*b*b*b*b*b*b*b*exp( 3.303 - sqrt( 10.909809 + 4.*v*v ) ) );
}

double Rpp_pACalc(double v){

    double c = 0.1747 + sqrt( 0.68112009 + 0.006293249*v*v );
    double d = 0.7333 + sqrt( 0.07112889 + 0.02815684*v*v );

    return 0.5*( c*c*exp( 4.228 - sqrt( 17.875984 + 4.*v*v ) ) + pow( d,7.5 )*exp( 7.762 - sqrt( 60.248644 + 9.*v*v ) ) );
}

/*
            kappa_eCore: electron thermal conductivity in reduced units

   Electron thermal conductivity from Shternin and Yakovlev PRD 75, 103004 (2007)

           Inputs: n_b: baryon number density in fm^{-3}
                   Y: proton fraction of total baryons
                   Y_e: electron fraction of total baryons Y_e = n_e/n_b. Equals Y in the absence of muons
           Output: kappa_e: electron contribution to core thermal conductivity in reduced units. Temperature-independent.
*/
double kappa_eCore(double n_b, double Y, double Y_e){

    double zeta_3 = 1.20205690316; //Riemann zeta function zeta(3)
    double p_Fe = pow( 3.*pi*pi*Y_e*n_b, 1./3. )*hbarc; //electron Fermi momentum in MeV
    double p_Fm = 0.;
    if( Y_e < Y ) p_Fm = pow( 3.*pi*pi*(Y-Y_e)*n_b, 1./3. )*hbarc; //muon Fermi momentum in MeV if nonzero

    //Thermal conductivity of electrons in the core. Assumes no superfluid/superconducting component and no muons
    double kappa_e = pi*pi/(54.*zeta_3)*k_Bcgs*c*p_Fe*p_Fe/(hbarc*hbarc*alpha_e/1e26); // in erg/cm/K/s. 1e26 factor converts fm^{-2} to cm^{-2}. Independent of T
    double kappa_m = pi*pi/(54.*zeta_3)*k_Bcgs*c*p_Fm*p_Fm/(hbarc*hbarc*alpha_e/1e26); // in erg/cm/K/s. 1e26 factor converts fm^{-2} to cm^{-2}. Independent of T

    return (kappa_e + kappa_m)*t_0/(s_0*L_0*L_0); //converts to reduced units
}

/*
            kappa_nCore: neutron thermal conductivity in reduced units evaluated at T=T_0

   Neutron thermal conductivity from Baiko, Haensel and Yakovlev, A&A 374, 151–163 (2001)

           Inputs: n_b: baryon number density in fm^{-3}
                   Y: proton fraction
                   M_nEff: neutron Dirac effective mass in MeV
                   M_pEff: proton Dirac effective mass in MeV
           Output: kappa_n: neutron contribution to core thermal conductivity in reduced units evaluated at T_0
*/
double kappa_nCore(double n_b, double Y, double M_nEff, double M_pEff){

    static double hbar = 6.582119569e-22; //hbar in MeV*s
    double k_Fn, k_Fp, u, Kn2, Kp1, Kp2, Sn2, Sp1, Sp2;
    double nu_nn, nu_np, tau_n, kappa_n;

    //Thermal conductivity of neutrons in the core. Assumes no superfluid/superconducting component and no muons

    k_Fn = pow( 3.*pi*pi*n_b*(1.-Y),1./3. ); // neutron Fermi wave number in fm^{-1}
    k_Fp = pow( 3.*pi*pi*n_b*Y,1./3. ); // proton Fermi wave number in fm^{-1}

    u = k_Fn - 1.556;
    Kn2 = pow( M_n/M_nEff,2. )*( 0.4891 + 1.111*u*u - 0.2283*pow(u,3.) + 0.01589*k_Fp - 0.02099*k_Fp*k_Fp + 0.2773*u*k_Fp );
    u = k_Fn - 2.126;
    Kp1 = pow( M_p/M_pEff,2. )*( 0.04377 + 1.100*u*u + 0.1180*pow(u,3.) + 0.1626*k_Fp + 0.3871*u*k_Fp - 0.2990*pow(u,4.) );
    u = k_Fn - 2.116;
    Kp2 = pow( M_p/M_pEff,2. )*( 0.0001313 + 1.248*u*u + 0.2403*pow(u,3.) + 0.3257*k_Fp + 0.5536*u*k_Fp - 0.3237*pow(u,4.) + 0.09786*u*u*k_Fp );
    Sn2 = 7.880/pow( k_Fn,1.5 )*( 1 - 0.0788*k_Fn + 0.0883*pow(k_Fn,3.) )/(1 - 0.1114*k_Fn ); //in mb
    Sp1 = 0.8007*k_Fp/pow( k_Fn,2. )*( 1 + 31.28*k_Fp - 0.0004285*k_Fp*k_Fp + 26.85*k_Fn + 0.08012*k_Fn*k_Fn )/( 1 - 0.5898*k_Fn + 0.2368*k_Fn*k_Fn + 0.5838*k_Fp*k_Fp + 0.884*k_Fn*k_Fp ); //in mb
    Sp2 = 0.3820*pow( k_Fp,4. )/pow( k_Fn,5.5 )*( 1 + 102.0*k_Fp + 53.91*k_Fn )/( 1 - 0.7087*k_Fn + 0.2537*k_Fn*k_Fn + 9.404*k_Fp*k_Fp - 1.589*k_Fn*k_Fp ); //in mb

    nu_nn = 64./5.*M_nEff*pow( M_nEff/M_n,2. )*pow( k_B*T_0,2. )*Kn2*Sn2/10./pow(hbarc,2.)/hbar; //n-n collision frequency in s^{-1}
    nu_np = 64./5.*M_pEff*(M_nEff/M_n)*(M_pEff/M_p)*pow( k_B*T_0,2. )*(Kp1*Sp1 + Kp2*Sp2)/10./pow(hbarc,2.)/hbar; //n-p collision frequency in s^{-1}
    tau_n = 1./(nu_nn + nu_np);  //neutron relaxation time in s
    kappa_n = pi*pi*k_B*k_B*T_0*n_b*(1.-Y)/( 3.*M_nEff )*tau_n*MeVtoErg*1e39*c*c;// in erg/cm/K/s. Because tau_n ~ 1/T^2, scales as 1/T

    return kappa_n*t_0/(s_0*L_0*L_0); //converts to reduced units
}
