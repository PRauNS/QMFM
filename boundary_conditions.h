#ifndef BOUNDARY_CONDITIONS_H_INCLUDED
#define BOUNDARY_CONDITIONS_H_INCLUDED

void B_Symmetrize(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm);
void B_BoundaryConditions(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm);
void B_BoundaryConditionsIntermediate(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm);
void Btheta_BC_Calc(std::vector<double> & Br_BC, std::vector<double> & Btheta_BC, size_t Ntheta, std::vector<int> & Ntheta_locs, std::vector<int> & starts, double n_r, int world_rank, MPI_Comm comm1D);
void Bphi_Outer_Sym(std::vector<double> & Bphi_outer, std::vector<double> & Bphi_outer_sym, size_t Ntheta, std::vector<int> & Ntheta_locs, std::vector<int> & starts, int world_rank, MPI_Comm comm1D);
void Bphi_BoundaryConditions(VectorField & B, const BParams & bparams, size_t Ntheta, size_t N_GC, const Process & process, std::vector<int> & Ntheta_locs, std::vector<int> & starts, const Domain & dm);
void E_BoundaryConditions(VectorField & phE, VectorField & cE, const BParams & bparams, size_t N_GC, const Process & process, const Domain & dm);
void T_BoundaryConditions(ScalarField & T, VectorField & q, double Tcore, TransCoeffs & tC, MagCoeffs & mC, TParams & tparams, const Domain & domain, size_t N_GC, const Process & process);
void T_sGudmundsson1983(double Tb8, double g_s14, double & Ts, double & Ts_star, double & dTsdTb);
void T_sPotekhin2001(double Tb8, double g_s14, double & Ts, double & Ts_star, double & dTsdTb);
void T_sPotekhinMag2001(double Tb9, double B12, double g_s14, double theta_B, double & Ts, double & Ts_star, double & dTsdTb);
void T_sPotekhin2015(double Tb9, double B12, double g_s14, double theta_B, double & Ts, double & Ts_star, double & dTsdTb);

void T_InnerBC_r(double & b_0_term, double & d_0_term, double Tcore, const Domain & dm, const CoreThermalParams & corethermalparams, double a_ii);
//void T_InnerBC_r(double & Tcore_update, double Tcore, const Domain & dm, const CoreThermalParams & corethermalparams, double a_ii);
void T_OuterBC_r(size_t j, double & b_n_term, double & d_n_term, ScalarField & T, TransCoeffs & tC, MagCoeffs & mC, const TParams & tparams, const Domain & dm, size_t N_GC);
void T_InteriorBC_theta(size_t j_in, ScalarField & T, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const Domain & dm, size_t N_GC, double a_ii, ScalarField & T_f);
// double T_InteriorBC_theta(size_t i_in, size_t j_in, ScalarField & T, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const Domain & dm, size_t N_GC, double a_ii);

#endif // BOUNDARY_CONDITIONS_H_INCLUDED
