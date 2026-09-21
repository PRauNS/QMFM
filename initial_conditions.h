#ifndef INITIAL_CONDITIONS_H_INCLUDED
#define INITIAL_CONDITIONS_H_INCLUDED

struct muConstraint_params{
    double Ri, Ro, n_r;
};

struct BConfig_params{
    double B_pol, B_tor_max, Ri, Ro, mu_tilde, n_r;
};

struct TConfig_params{
    double T_init, Ri, Ro;
};

void InitializeB(std::vector<double> & r, std::vector<double> & theta, const BParams & bparams, const Domain & domain, size_t N_GC, std::vector<double> & Deltar, double Deltatheta, VectorField & B);
double f_r(double mutilde, double Ro, double n_r, double r);
double df_rdr(double mutilde, double Ro, double n_r, double r);
double muConstraint(double mutilde, void *params);
double muConstraintf_r(double mutilde, double R_ratio, double n_r);
double InitialBr_r(double r, void * params);
double InitialBr_theta(double theta, void * params);
double InitialBr_theta_analytic(double theta1, double theta2);
double InitialBtheta_r(double r, void * params);
double InitialBtheta_theta(double theta, void * params);
double InitialBphi_r(double r, void * params);
double InitialBphi_theta(double theta, void * params);
double InitialBphi_r_analytic(double r1, double r2, void * params);
double InitialBphi_theta_analytic(double theta1, double theta2);

void InitializeT(std::vector<double> & r, std::vector<double> & theta, const TParams& tparams, const Domain & domain, size_t N_GC, std::vector<double> & Deltar, double Deltatheta, ScalarField & T);
double InitialT_r(double r, void * params);

#endif // INITIAL_CONDITIONS_H_INCLUDED
