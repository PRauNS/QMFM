#ifndef CONSERVATION_CHECKS_H_INCLUDED
#define CONSERVATION_CHECKS_H_INCLUDED

double divB_Calculator(VectorField & B, size_t N_GC, const Domain & dm);
void divB_Monitor(VectorField & B, size_t N_GC, const Domain & dm, MPI_Comm comm1D, int world_rank, std::vector<int> & Ntheta_locs);
void EnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & q_SH, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & JouleH, double & PoyntingF, double & SHSum);
void QuasiEnergyConservation(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & q_SH, TransCoeffs & tC, MagCoeffs & mC, size_t N_GC, const Domain & dm, double & U_B, double & QuasiJouleH, double & PoyntingF,
    double & JouleH, double & SH);
void T_avCalc(ScalarField & T, size_t N_GC, double & T_av, double & T_min, double & T_max);
void Qnu_CrustIntegratedCalc(ScalarField & q_nuCrust, size_t N_GC, const Domain & dm, double & Q_nuCrustIntegrated);

void EnergyConservationLocal_Initialize(VectorField & B, ScalarField & UB_Init, ScalarField & JHPF_Cumulative, size_t N_GC, const Domain & dm, double t, int world_rank);
void EnergyConservationLocal(VectorField & B, VectorField & phE, VectorField & cE, VectorField & J, ScalarField & DeltaE, ScalarField & UB_Init, ScalarField & JHPF_Cumulative, TransCoeffs & tC, MagCoeffs & mC,
                            size_t N_GC, const Domain & dm, double t, const Process & process, int world_rank, size_t iter, size_t ECons_cadence, const BParams & bparams);

#endif // CONSERVATION_CHECKS_H_INCLUDED
