#ifndef MICROPHYSICS_H_INCLUDED
#define MICROPHYSICS_H_INCLUDED

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

class EOSInterpolation{
	public:
        gsl_spline* n_e_spline;
        gsl_spline* A_spline;
        gsl_spline* Z_spline;
        gsl_spline* n_i_spline;
        gsl_spline* n_b_spline;
        gsl_spline* rho_spline;
        gsl_spline* n_nf_spline;
        gsl_spline* mu_nf_spline;
        gsl_spline* gtt_spline;
        gsl_spline* grr_spline;
        gsl_interp_accel* n_e_acc;
        gsl_interp_accel* A_acc;
        gsl_interp_accel* Z_acc;
        gsl_interp_accel* n_i_acc;
        gsl_interp_accel* n_b_acc;
        gsl_interp_accel* rho_acc;
        gsl_interp_accel* n_nf_acc;
        gsl_interp_accel* mu_nf_acc;
        gsl_interp_accel* gtt_acc;
        gsl_interp_accel* grr_acc;
        double R_cc; //radius of crust-core transition in cm
        double RStar; //radius of star in cm
        double R_rhocutoff; //radius of star at cutoff density (minimum density considered in crust) in cm

    ~EOSInterpolation() {
        gsl_spline_free (n_e_spline);
        gsl_interp_accel_free (n_e_acc);
        gsl_spline_free (A_spline);
        gsl_interp_accel_free (A_acc);
        gsl_spline_free (Z_spline);
        gsl_interp_accel_free (Z_acc);
        gsl_spline_free (n_i_spline);
        gsl_interp_accel_free (n_i_acc);
        gsl_spline_free (n_b_spline);
        gsl_interp_accel_free (n_b_acc);
        gsl_spline_free (rho_spline);
        gsl_interp_accel_free (rho_acc);
        gsl_spline_free (n_nf_spline);
        gsl_interp_accel_free (n_nf_acc);
        gsl_spline_free (mu_nf_spline);
        gsl_interp_accel_free (mu_nf_acc);
        gsl_spline_free (gtt_spline);
        gsl_interp_accel_free (gtt_acc);
        gsl_spline_free (grr_spline);
        gsl_interp_accel_free (grr_acc);
    }

};

void load_EOS(std::string CrustEOS_filename, double rho_cutoff, double & g_s14, double & r_s, double & n_t_s, double & n_b_nd, EOSInterpolation & Interpolators);
void GR_Factor_Initialize(bool GR, EOSInterpolation & EOS_Interps, std::vector<double> & r, size_t Nr, size_t N_GC, std::vector<double> & n_t, std::vector<double> & n_r);
void SFgaps(TParams & tparams);
void ConductivityCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC, std::vector<double> & n_t,
                ScalarField & c_v, ScalarField & eta_O, ScalarField & kappa, double & max_eta_T);
double expE1_rp(double x);
double Q_impVigano2013Q100(double n_b);
double Q_impCarreau2020BSk24(double n_b, double n_drip);

void c_vfuncCalc(ScalarField & T, ScalarField & n_e, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, RadialScalarField & n_nf, RadialScalarField & mu_nf, const TParams & tparams, size_t N_GC, std::vector<double> & n_t, ScalarField & c_v);
void q_nuCrustCalc(ScalarField & T, ScalarField & n_e, ScalarField & Bmag, RadialScalarField & A, RadialScalarField & Z, RadialScalarField & n_i, RadialScalarField & n_b, ConductParams & cparams, const TParams & tparams, size_t N_GC, std::vector<double> & n_t,
                    bool IMEX, ScalarField & q_nuCrust, ScalarField & dq_nuCrustdT);
double CoulombLogf(double x);
double erf(double x);

void ThermPhysCore(std::string CoreEOS_filename, bool GR, const TParams & tparams, double T, double & C_v, double & Q_nu, double & kappa_eCC, double & kappa_nCC);
double c_vCore(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams);
double q_nuMUrca(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams);
double q_nuDUrca(double T, double n_b, double Y, double Y_e, double M_nEff, double M_pEff, const TParams & tparams);
double q_nu_NNbrems(double T, double n_b, double Y, double M_nEff, double M_pEff, const TParams & tparams);
double q_nuPBF(double T, double n_b, double Y, double M_nEff, const TParams & tparams);
void Compute_CoreThermalParams_Interpolators(std::vector<double> & TVec, std::vector<double> & C_vVec, std::vector<double> & Q_nuVec, CoreThermalParams & corethermalparams);
double Rn_AB_RegionI(double v_n, double v_p);
double Rn_AB_RegionII(double v_n, double v_p);
double Rn_AB_RegionIII(double v_n, double v_p);
double Rn_AB_RegionIV(double v_n, double v_p);
double Rp_AB_RegionI(double v_n, double v_p);
double Rp_AB_RegionII(double v_n, double v_p);
double Rp_AB_RegionIII(double v_n, double v_p);
double Rp_AB_RegionIV(double v_n, double v_p);
double RD_ACalc(double v_p);
double RD_BACalc(double v_n, double v_p);
double Rnp_pACalc(double v);
double Rpp_pACalc(double v);

double kappa_eCore(double n_b, double Y, double Y_e);
double kappa_nCore(double n_b, double Y, double M_nEff, double M_pEff);

#endif // MICROPHYSICS_H_INCLUDED
