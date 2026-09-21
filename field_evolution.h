#ifndef FIELD_EVOLUTION_H_INCLUDED
#define FIELD_EVOLUTION_H_INCLUDED

void Compute_J(VectorField & B, VectorField & J, ScalarField & T, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process);
void Compute_J_Poloidal(VectorField & B, VectorField & J, ScalarField & T, MagCoeffs & mC, size_t N_GC, const Domain & domain, const Process & process);
void B_torEvolve(ScalarField & Qphi, ScalarField & q_SH, VectorField & B, VectorField & J, ScalarField & T, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process);
double Godunov(double u_L, double u_R, double beta_L, double B_L, double beta_R, double B_R);
double Reconstruct_Btheta(VectorField & B, size_t i, size_t j, double vr_av_phi, double Deltar_1, double Deltar_2, double Deltar_3, const Domain & dm);
double Reconstruct_Br(VectorField & B, size_t i, size_t j, double vtheta_av_phi, double Deltatheta_1, double Deltatheta_2, double Deltatheta_3, const Domain & dm);
double Reconstruct_BrEth(VectorField & B, size_t i, size_t j, double vtheta_av_phi, double Deltatheta_1, double Deltatheta_2, double Deltatheta_3, const Domain & dm);
void Compute_E(VectorField & B, VectorField & Bn, VectorField & phE, VectorField & cE, VectorField & J, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process);
void Compute_EMF(ScalarField & Qr, ScalarField & Qtheta, VectorField & phE, size_t N_GC, const Domain & dm, const Process & process);
void Compute_q(ScalarField & T, VectorField & q, TransCoeffs & trC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process);
void TEvolve(ScalarField & QT, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH, TransCoeffs & trC, ThermCoeffs & thC, TParams & tparams, size_t N_GC, const Domain & dm);
double JouleHeating(size_t i, size_t j, VectorField & J, VectorField & phE, VectorField & cE, TransCoeffs & tC, const Domain & dm);
double QuasiJouleHeating(size_t i, size_t j, VectorField & J, VectorField & phE, VectorField & cE, TransCoeffs & tC, MagCoeffs & mC, const Domain & dm);
double TCoreEvolve(ScalarField & T, double Tcore, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, const CoreThermalParams & corethermalparams);
double TCore_starEvolve(ScalarField & T, double Tcore, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, const CoreThermalParams & corethermalparams, double a_ijIm);
double TCoreEvolveImplicit(ScalarField & T, double Tcore, double Tcore_init, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, const Process & process, CoreThermalParams & corethermalparams, double a_ii);

void Compute_Tstar_IMEX(ScalarField & T_star, ScalarField & T_init, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH,
                        ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, TParams & tparams, size_t N_GC, const Domain & dm, const Process & process, double a_ijIm);
void Compute_KT_IMEX(ScalarField & KT, ScalarField & T, VectorField & q, VectorField & J, VectorField & phE, VectorField & cE, ScalarField & q_SH, TransCoeffs & trC, ThermCoeffs & thC, MagCoeffs & mC, TParams & tparams, size_t N_GC, const Domain & dm, const Process & process);
void Compute_Tint_IMEX(ScalarField & T_int, ScalarField & T, ScalarField & T_star, ScalarField & KT, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const TParams & tparams, const CoreThermalParams & corethermalparams, size_t N_GC,
                        const Domain & dm, const Process & process, double Tcore, double a_ii);
void Compute_T_IMEX(ScalarField & T_f, ScalarField & T_int, ScalarField & T, ThermCoeffs & thC, TransCoeffs & trC, MagCoeffs & mC, const Process & process, size_t N_GC, const Domain & dm, double a_ii);

#endif // FIELD_EVOLUTION_H_INCLUDED
