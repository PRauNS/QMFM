#ifndef MICROPHYSICS_BDEP_H_INCLUDED
#define MICROPHYSICS_BDEP_H_INCLUDED

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class OmegaB_deriv_Interpolation {
	public:
        gsl_spline* PolyLogN1_2; //Polylogarithm_{-1/2}
        gsl_spline* PolyLog1_2; //Polylogarithm_{1/2}
        gsl_spline* PolyLog3_2; //Polylogarithm_{3/2}
        gsl_spline* PolyLog5_2; //Polylogarithm_{5/2}
        gsl_spline* PolyLog7_2; //Polylogarithm_{7/2}
        gsl_spline2d* LogG_1; //ln(G_1(a,b))
        gsl_spline2d* LogG_2; //ln(G_2(a,b))
        gsl_spline2d* LogNH_2; //ln(-H_2(a,b))
        gsl_spline2d* LogI_1; //ln(I_1(a,b))
        gsl_spline2d* LogNI_2; //ln(-I_1(a,b))
        gsl_spline2d* h_1; //h_1(a,b)
        gsl_spline2d* h_2; //h_2(a,b)
        gsl_spline* i_1; //i_1(x)
        gsl_spline* i_2; //i_2(x)
        gsl_spline* i_3; //i_3(x)
        gsl_interp_accel* PolyLogN1_2_acc;
        gsl_interp_accel* PolyLog1_2_acc;
        gsl_interp_accel* PolyLog3_2_acc;
        gsl_interp_accel* PolyLog5_2_acc;
        gsl_interp_accel* PolyLog7_2_acc;
        gsl_interp_accel* LogG_1_acc_a;
        gsl_interp_accel* LogG_1_acc_b;
        gsl_interp_accel* LogG_2_acc_a;
        gsl_interp_accel* LogG_2_acc_b;
        gsl_interp_accel* LogNH_2_acc_a;
        gsl_interp_accel* LogNH_2_acc_b;
        gsl_interp_accel* LogI_1_acc_a;
        gsl_interp_accel* LogI_1_acc_b;
        gsl_interp_accel* LogNI_2_acc_a;
        gsl_interp_accel* LogNI_2_acc_b;
        gsl_interp_accel* h_1_acc_a;
        gsl_interp_accel* h_1_acc_b;
        gsl_interp_accel* h_2_acc_a;
        gsl_interp_accel* h_2_acc_b;
        gsl_interp_accel* i_1_acc;
        gsl_interp_accel* i_2_acc;
        gsl_interp_accel* i_3_acc;

    ~OmegaB_deriv_Interpolation() {
        gsl_spline_free (PolyLogN1_2);
        gsl_interp_accel_free (PolyLogN1_2_acc);
        gsl_spline_free (PolyLog1_2);
        gsl_interp_accel_free (PolyLog1_2_acc);
        gsl_spline_free (PolyLog3_2);
        gsl_interp_accel_free (PolyLog3_2_acc);
        gsl_spline_free (PolyLog5_2);
        gsl_interp_accel_free (PolyLog5_2_acc);
        gsl_spline_free (PolyLog7_2);
        gsl_interp_accel_free (PolyLog7_2_acc);
        gsl_spline2d_free (LogG_1);
        gsl_interp_accel_free (LogG_1_acc_a);
        gsl_interp_accel_free (LogG_1_acc_b);
        gsl_spline2d_free (LogG_2);
        gsl_interp_accel_free (LogG_2_acc_a);
        gsl_interp_accel_free (LogG_2_acc_b);
        gsl_spline2d_free (LogNH_2);
        gsl_interp_accel_free (LogNH_2_acc_a);
        gsl_interp_accel_free (LogNH_2_acc_b);
        gsl_spline2d_free (LogI_1);
        gsl_interp_accel_free (LogI_1_acc_a);
        gsl_interp_accel_free (LogI_1_acc_b);
        gsl_spline2d_free (LogNI_2);
        gsl_interp_accel_free (LogNI_2_acc_a);
        gsl_interp_accel_free (LogNI_2_acc_b);
        gsl_spline2d_free (h_1);
        gsl_interp_accel_free (h_1_acc_a);
        gsl_interp_accel_free (h_1_acc_b);
        gsl_spline2d_free (h_2);
        gsl_interp_accel_free (h_2_acc_a);
        gsl_interp_accel_free (h_2_acc_b);
        gsl_spline_free (i_1);
        gsl_interp_accel_free (i_1_acc);
        gsl_spline_free (i_2);
        gsl_interp_accel_free (i_2_acc);
        gsl_spline_free (i_3);
        gsl_interp_accel_free (i_3_acc);
    }

};

struct HighTInterp_params{
    double x; // x = (mu_e^2-M_e^2)/eB
};

void Compute_OmegaB_Interpolators(OmegaB_deriv_Interpolation & OmB_Interpolators);
double i_1Integrand(double y, void * params);
double i_2Integrand(double y, void * params);
double i_3Integrand(double y, void * params);

void Compute_mu_e(ScalarField & n_e, RadialScalarField & mu_e);
void ConductivityCalcB(ScalarField & T, ScalarField & Bmag, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC,
    std::vector<double> & n_t, ScalarField & c_v, ScalarField & eta_O_par, ScalarField & eta_O_perp, ScalarField & kappa_par, ScalarField & kappa_perp, ScalarField & kappa_H, double & max_eta_T, const Process & process, const Domain & dm);

void c_vfuncCalcB(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, RadialScalarField & n_nf, RadialScalarField & mu_nf, const TParams & tparams, size_t N_GC, std::vector<double> & n_t, ScalarField & c_v);

void Omega_xyFixedT(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & M, ScalarField & chi, ScalarField & M_mu);
void Omega_xyLowT_FixedT_EMSums(double eB, double mu, double T, double p_F, double n_max, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp);
void Omega_xyLowT_FixedT_MaxLL(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp);
void Omega_xyHighT_FixedT(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp);

void Omega_xyVaryingT(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & M, ScalarField & chi, ScalarField & M_mu, ScalarField & c_v_eB, ScalarField & C_m);
void Omega_xyLowT_VaryingT_EMSums(double eB, double mu, double T, double p_F, double n_max, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp);
void Omega_xyLowT_VaryingT_MaxLL(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp);
void Omega_xyHighT_VaryingT(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & dPdBTemp, double & d2PdB2Temp, double & d2PdBdmuTemp, double & d2PdT2Temp, double & d2PdTdBTemp);

void Omega_xyVaryingTC_veOnly(ScalarField & Bmag, RadialScalarField & mu, ScalarField & T, size_t N_GC, std::vector<double> & n_t, OmegaB_deriv_Interpolation & OmB_Interpolators, ScalarField & c_v_eB);
void Omega_xyLowT_VaryingT_EMSumsC_veOnly(double eB, double mu, double T, double p_F, double n_max, double & d2PdT2Temp);
void Omega_xyLowT_VaryingT_MaxLLC_veOnly(double eB, double mu, double T, double n_max, OmegaB_deriv_Interpolation & OmB_Interpolators, double & d2PdT2Temp);
void Omega_xyHighT_VaryingTC_veOnly(double eB, double mu, double T, double p_F, OmegaB_deriv_Interpolation & OmB_Interpolators, double & d2PdT2Temp);

#endif // MICROPHYSICS_BDEP_H_INCLUDED
