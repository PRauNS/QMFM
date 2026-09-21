#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <mpi.h>
#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "boost/multi_array.hpp"

extern const double pi;
extern const double c; //speed of light in cm/s
extern const double hbarc; //hbar times c in MeV*fm
extern const double unit_e; //elementary charge in units of statcoulomb
extern const double yr; //1 year in seconds
extern const double G; //Gravitational extern constant in dyn*cm^2/g^2
extern const double k_B; //Boltzmann extern constant in MeV/K
extern const double k_Bcgs; //Boltzmann extern constant in cgs units (erg/K)
extern const double M_e; //electron mass in MeV
extern const double M_m; //Muon mass in MeV (105.658375 MeV)
extern const double M_N; //nucleon mass in MeV
extern const double M_n; //neutron mass in MeV
extern const double M_p; //proton mass in MeV
extern const double M_solar; //solar mass in g
extern const double alpha_e; //electromagnetic fine structure extern constant (dimensionless)
extern const double B_crit; //quantum critical magnetic field (Schwinger field) in G
extern const double eB_crit; //critical magnetic field times elementary charge in MeV^2
extern const double sigma_SB; //Stefan-Boltzmann extern constant in erg/cm^2/s/K^4

extern const double MeVtoErg; //conversion factor from MeV to erg
extern const double statCGtoMeV2; //conversion factor from statCoulomb*Gauss to MeV^2: 1 statCoulomb*Gauss to __ MeV/fm*hbarc = ___ MeV^2
extern const double MeV2toGauss; //Conversion factor from MeV^2 to G: 1 MeV^2 = 4.002719868e16/(197.3269804)^(3/2) G = 1.444027592e13 G (hbarc=c=1 all energy units)

extern const double B_0; //characteristic magnetic field (G)
extern const double n_e0; //characteristic electron density (fm^{-3})
extern const double L_0; //characteristic length scale (cm)
extern const double t_0; //characteristic timescale (s). Taken as the Hall time-scale with length scale 10 m and field 10^{15} G
extern const double T_0; //characteristic temperature (K)
extern const double s_0; //characteristic entropy density (erg/K/cm^3)
extern const double E_0; //characteristic electric field (statV/cm)

typedef boost::multi_array<double, 3> VectorField; //type definition for three-component vector fields defined over two spatial dimensions. First index is components (0=x, 1=y, 2=z), second and third indices are x and y coordinates
typedef boost::multi_array_types::index_range range; //range of indices used in array slicing
typedef boost::multi_array<double, 2> ScalarField; //type definition for scalar fields defined over two spatial dimensions. Indices are x and y coordinates
typedef boost::multi_array<double, 1> RadialScalarField; //type definition for scalar fields defined over one ("radial" or x) spatial dimension. Index is x coordinate

struct Process {

    std::string Timestep_method;
    int world_rank, num_procs, nbrleft, nbrright;
    size_t MyS, MyE;
    MPI_Comm comm1D;

};

struct Domain {
    size_t Nr, Ntheta, N_GC;
    double Lr, Ltheta, Deltatheta;
	double Deltatheta_dfactor;
    std::vector<double> Deltar;
    std::vector<double> r, theta;
    std::vector<double> n_t, n_r; //lapse function n_t = sqrt(-g_{tt}) and square root of radial part of metric tensor n_r = sqrt(g_{rr})
    double n_router; //n_r = sqrt(g_{rr}) at r=Ro
    double r_min, r_max; //inner and outer radius of simulation domain in reduced units.
    double r_star; //radius of star in reduced units.
    double Deltat;
    std::vector<int> Ntheta_locs, starts;
    double t;
    bool InitiallyEquatoriallySymmetric; //whether initial magnetic field profile is equatorially symmetric (true) or not (false)
    MPI_Datatype stridetype_Vec, stridetype_Vec_q, stridetype_Sca; //stridetypes for exchanging data between parallelized partial domains for VectorField and ScalarField objects
};

struct SimParams {

    std::string CrustEOS; //File name for crust EOS data table
    std::string CoreEOS; //File name for core EOS data table
    std::string Timestep_method; //timestepping method
    bool IMEX; //true if timestepping method is an IMEX method, false otherwise
    bool varying_mesh = false; //true for varying cell-center spacing in "radial"- (x-)direction or false otherwise. Defaults to false if unspecified.
    size_t Nr = 0, Ntheta = 0; //resolution in x-and y-directions
    double r_min = -1., r_max = -1., theta_min = 1., theta_max = -1.; //limits of simulation domain in reduced units
    double t_max = 0.; //maximum time in reduced units
    double k_CB = 0.; //Courant number for magnetic evolution
    double k_CT = 0.; //Courant number for thermal evolution
    double rho_cutoff = 0.; //cutoff (energy) density in g/cm^3 (lowest density to include in simulation domain)
    size_t saves_number = 0; //maximum number of snapshots to save to H5 file
    size_t ECons_cadence = 0; //cadence to print energy conservation information to terminal
    bool GR = true; //true for general relativity turned on or false otherwise. Defaults to true if unspecified.
    std::string OutputFile; //File name for output H5 file (excluding file extension)
    bool divBCheck = false; //Explicitly show div(B) check every ECons_cadence. Defaults to false if unspecified.
    double C_hyp = 0.1; //toroidal B hyperdiffusivity dimensionless prefactor <1 (usually <~0.1, default value).
    bool InitiallyEquatoriallySymmetric; //whether initial magnetic field profile is equatorially symmetric (true) or not (false)
};

struct BParams {

    double B_pol_init = 1., B_tor_init = 0.5; //initial magnitude of poloidal and toroidal fields in reduced units.
    std::string B_perp_lower, B_parallel_lower, B_pol_upper, B_tor_upper; //boundary conditions for the magnetic field
    std::string E_perp_lower, E_parallel_lower, E_perp_upper, E_parallel_upper; //lower boundary conditions for electric field
    bool quantization = false; //Whether to include Landau quantization effects (true) or not (false). Defaults to false if unspecified.
    std::string Z_impurity = "zero"; //form of impurities to include. Currently implement three models: "zero" (Q=0, no impurities, default), "Vigano2013Q100" (Q modeled after impurity parameter used in Vigano 2013, PhD thesis, University of Alacant, Spain) or "Carreau2020BSk24" (Q for BSk24 equation of state from Carreau et al. A&A, A77 (2020))
    bool PolarAxisQuantization = true; //Whether to include Landau quantization effects in current at polar axis (true) or not (false). Defaults to true if unspecified.

};

struct TParams {

    double T_init = 3e8/T_0; //initial temperature value. Defaults to 3e8 K (converted to reduced units)
    double g_s14 = 1.; //gravitational acceleration at surface of the star in units of 1e14 cm/s^2
    double r_s = 12.; //outer (surface) radius of star in reduced units
    double n_t_s = 1.; //lapse function exp(nu/2) at outer radius of star
    bool uniform_T = true; //true for uniform initial temperature profile, false for non-uniform initial temperature profile.
    bool fixed_T = false; //true for fixed (non-evolving) temperature profile, false for evolving temperature with time.
    double JouleHeating = 1.; //1 for Joule heating turned on and 0 for it turned off. Defaults to 1 if unspecified.
    bool SF = false; //Whether to include superfluidity (true) or not (false). Defaults to false if unspecified.
    bool InMediumMUrca = false; //Whether to include in-medium corrections to the modified Urca neutrino emissivity in the core (true) or not (false). Defaults to false if unspecified.
    std::string T_lower, T_upper; //boundary conditions for temperature
    bool conductivity_anisotropy = false; //Whether to include anisotropic conductivity (true) or not (false). Defaults to false if unspecified.
    double Ts4_av = 1e7/T_0; //average (local surface temperature)^4 in reduced units across part of surface within current process. Updated by T_BoundaryConditionsInit and T_BoundaryConditions functions.
    double Tcore; //redshifted temperature of the isothermal core in reduced units. Updated within RKStep function in main.cpp. Only contains values for full timesteps (not intermediate values).
    double Qnu_Crust; //partial domain-integrated neutrino emissivity in
    double B_pol; //magnitude of magnetic field at poles in reduced units. Updated within RKStep function in main.cpp
    //Models for superfluid gap to use. Defaults to "SFB" for neutrons in crust, "CCDK" for protons in core, "TToa" for neutrons in core. Models defined in Ho et al. PRC 91, 015806 (2015) or Ho, Glampedakis and Andersson MNRAS 422, 2632 (2012).
    std::string SF_n_crust = "SFB"; //crust neutrons (1S0 pairing)
    std::string SF_p_core = "CCDK"; //core protons (1S0 pairing)
    std::string SF_n_core = "TToa"; //core neutrons (3P2, mJ=0 pairing)
    double Delta0_nCrust = 45., k0_nCrust = 0.1, k1_nCrust = 4.5, k2_nCrust = 1.55, k3_nCrust = 2.5; //parametrization of neutron 1S0 superfluid gap in crust. Delta0 in MeV, k0 and k2 in fm^-1, k1 and k2 in fm^-2. Defaults to "SFB" values.
    double Delta0_pCore = 102., k0_pCore = 0., k1_pCore = 9., k2_pCore = 1.3, k3_pCore = 1.5; //parametrization of proton 1S0 superfluid gap in core. Delta0 in MeV, k0 and k2 in fm^-1, k1 and k2 in fm^-2. Defaults to "CCDK" values.
    double Delta0_nCore = 2.1, k0_nCore = 1.1, k1_nCore = 0.6, k2_nCore = 3.2, k3_nCore = 2.4; //parametrization of neutron 3P2 superfluid gap in core. Delta0 in MeV, k0 and k2 in fm^-1, k1 and k2 in fm^-2. Defaults to "TToa" values.

};

struct ConductParams {

    std::string EOS; //which EOS is being used. Filled in "load_EOS" function
    std::string EorT = "both"; //whether to conduct electrical conductivity parameters ("Electrical"), thermal conductivity parameters ("Thermal"), or both ("Both", default)
    bool quantization = false; //Whether to include Landau quantization effects (true) or not (false). Defaults to false if unspecified.
    std::string Z_impurity = "zero"; //form of impurities to include. Currently implement two models: "zero" (Q=0, no impurities, default) or "Vigano2013" (Q modeled after impurity parameter used in Vigano 2013, PhD thesis, University of Alacant, Spain).
    double n_b_nd = 1e-3; //neutron drip baryon number density in fm^-3. Default to 1e-3 fm^-3 (too high, set correctly based on EOS table in main)

};

class CoreThermalParams {
	public:
        gsl_spline* C_vCore; //Core-integrated specific heat capacity fit
        gsl_spline* Q_nuCore; //Core-integrated total neutrino emissivity fit
        gsl_interp_accel* C_vCore_acc;
        gsl_interp_accel* Q_nuCore_acc;
        double HeatFlux; //radial heat flux from core to crust. Uses temperature from previous time step/partial time step.
        double HeatFlux_pref; //prefactor (of temperature at next time step) of radial temperature derivative term in heat flux. Used for implicit timestepping.
        double HeatFlux_H; //thermal Hall contribution to heat flux. Uses temperature from previous time step/partial time step.
        double kappa_eCore; //core electron thermal conductivity at crust-core boundary in reduced units. Independent of temperature
        double kappa_nCore; //core neutron thermal conductivity at crust-core boundary in reduced units evaluated at T=T_0. Scales with temperature as 1/T.
        mutable double Q_nuCoreVal; //core-integrated total neutrino emissivity in reduced units. Updated

    ~CoreThermalParams() {
        gsl_spline_free (C_vCore);
        gsl_interp_accel_free (C_vCore_acc);
        gsl_spline_free (Q_nuCore);
        gsl_interp_accel_free (Q_nuCore_acc);
    }

};

/*
    TransCoeffs class

    Members are ScalarFields containing transport coefficients

    Class objects are instantiated with radial Nr and angular Nth sizes of the ScalarFields
*/
class TransCoeffs {
    public:
        size_t Nr, Nth;
        bool cond_aniso_In; //input variable for conductivity_anisotropy boolean
        bool conductivity_anisotropy; //permanent variable for conductivity_anisotropy boolean
        double C_hyp_In; //input variable for C_hyp
        double C_hyp; //toroidal B hyperdiffusivity dimensionless prefactor <1 (usually <~0.1)

    ScalarField eta_H; //
    ScalarField psi_H; //Hall stream function for calculation of advection velocities in toroidal B evolution
    ScalarField eta_O; //
    ScalarField eta_O_delta; //
    ScalarField eta_O_perp; //
    ScalarField kappa; //
    ScalarField kappa_delta; //
    ScalarField kappa_perp; //
    ScalarField kappa_H; //
    double max_eta_local;
    double kappa_eCore; //core electron thermal conductivity at crust-core boundary in reduced units. Independent of temperature
    double kappa_nCore; //core neutron thermal conductivity at crust-core boundary in reduced units evaluated at T=T_0. Scales with temperature as 1/T.

    TransCoeffs(size_t Nr, size_t Nth, bool cond_aniso_In, double C_hyp_In) :
        eta_H(boost::extents[Nr][Nth]), //
        psi_H(boost::extents[Nr][Nth]), //
        eta_O(boost::extents[Nr][Nth]), //
        eta_O_delta(boost::extents[Nr][Nth]), //
        eta_O_perp(boost::extents[Nr][Nth]), //
        kappa(boost::extents[Nr][Nth]), //
        kappa_delta(boost::extents[Nr][Nth]),//
        kappa_perp(boost::extents[Nr][Nth]),//
        kappa_H(boost::extents[Nr][Nth]) //
    {
        conductivity_anisotropy = cond_aniso_In;
        C_hyp = C_hyp_In;
        if(conductivity_anisotropy == true){
            //resize isotropic conductivities to zero if unused
            eta_O.resize(boost::extents[0][0]);
            kappa.resize(boost::extents[0][0]);
        }
        else{
            //resize anisotropic conductivities to zero if unused
            eta_O_delta.resize(boost::extents[0][0]);
            eta_O_perp.resize(boost::extents[0][0]);
            kappa_delta.resize(boost::extents[0][0]);
            kappa_perp.resize(boost::extents[0][0]);
            kappa_H.resize(boost::extents[0][0]);
        }
    }
};

/*
    ThermCoeffs class

    Members are ScalarFields containing thermodynamic coefficients

    Class objects are instantiated with radial Nr and angular Nth sizes of the ScalarFields
*/
class ThermCoeffs {
    public:
    size_t Nr, Nth;
    bool IMEX_In;
    bool IMEX; //if using an IMEX timestepper, need to compute dq_nuCrust/dT

    ScalarField c_v;
    ScalarField q_nuCrust;
    ScalarField dq_nuCrustdT;

    ThermCoeffs(size_t Nr, size_t Nth, bool IMEX_In) :
        c_v(boost::extents[Nr][Nth]), //specific heat capacity in reduced units
        q_nuCrust(boost::extents[Nr][Nth]), //neutrino emissivity in reduced units
        dq_nuCrustdT(boost::extents[Nr][Nth]) //temperature derivative of neutrino emissivity in reduced units
    {
        IMEX = IMEX_In;
        if(IMEX == false){
            dq_nuCrustdT.resize(boost::extents[0][0]);
        }
    }
};

/*
    MagCoeffs class

    Members are ScalarFields containing magnetic field magnitude and magnetization coefficients
    plus VectorField Bhat, the unit vector in the direction of the  magnetic field, and the
    electron chemical potential mu_e, which is a RadialScalarField

    Class objects are instantiated with radial Nr and angular Nth sizes of the ScalarFields
*/
struct MagCoeffs {
    public:
    size_t Nr, Nth;
    bool quantization_In; //input variable for magnetization boolean
    bool quantization; //permanent variable for magnetization boolean
    bool PolarAxisQuantization_In; //input variable for polar axis quantization boolean
    bool PolarAxisQuantization; //permanent variable for polar axis quantization boolean

    ScalarField Bmag; //magnetic field magnitude (reduced units)
    ScalarField vtheta_phi_lb_Lag;
    ScalarField vtheta_phi_lt_Lag;
    ScalarField vtheta_phi_rb_Lag;
    ScalarField vtheta_phi_rt_Lag;
    VectorField rBhat; //magnetic field unit vector at location of J_r
    VectorField thBhat; //magnetic field unit vector at location of J_theta
    VectorField phBhat; //magnetic field unit vector at location of J_phi
    VectorField J_B; //curl(B) (reduced units)
    VectorField E_H_Pol; //Cell-edge centered poloidal Hall electric field. Only contains E_r and E_theta.
    VectorField phE_H; //Poloidal Hall electric field at cell corners. Only contains E_r and E_theta.
    ScalarField M; //magnetization (reduced units)
    ScalarField chi; //magnetic susceptibility
    ScalarField M_mu; //mixed susceptibility (reduced units)
    ScalarField C_m; //magnetocaloric coefficient (reduced units)
    RadialScalarField mu_e; //electron chemical potential (MeV)

    MagCoeffs(size_t Nr, size_t Nth, bool quantization_In, bool PolarAxisQuantization_In) :
        Bmag(boost::extents[Nr][Nth]), //magnitude of magnetic field in reduced units
        vtheta_phi_lb_Lag(boost::extents[Nr][Nth]), //magnitude of magnetic field in reduced units
        vtheta_phi_lt_Lag(boost::extents[Nr][Nth]), //magnitude of magnetic field in reduced units
        vtheta_phi_rb_Lag(boost::extents[Nr][Nth]), //magnitude of magnetic field in reduced units
        vtheta_phi_rt_Lag(boost::extents[Nr][Nth]), //magnitude of magnetic field in reduced units
//        Bhat(boost::extents[3][Nr][Nth]), //magnetic field unit vector at cell centers (location of B_phi)
        rBhat(boost::extents[3][Nr][Nth]), //magnetic field unit vector at location of J_r
        thBhat(boost::extents[3][Nr][Nth]), //magnetic field unit vector at location of J_theta
        phBhat(boost::extents[3][Nr][Nth]), //magnetic field unit vector at location of J_phi
        J_B(boost::extents[3][Nr][Nth]), //curl(B) in reduced units
        E_H_Pol(boost::extents[3][Nr][Nth]), //Cell-edge centered poloidal Hall electric field. Only contains E_r and E_theta.
        phE_H(boost::extents[3][Nr][Nth]), //Poloidal Hall electric field at cell corners. Only contains E_r and E_theta.
        M(boost::extents[Nr][Nth]), //magnetization in reduced units
        chi(boost::extents[Nr][Nth]), //differential magnetic susceptibility (dimensionless)
        M_mu(boost::extents[Nr][Nth]), //mixed susceptibility (in MeV^{-1})
        C_m(boost::extents[Nr][Nth]), //magnetocaloric coefficient in reduced units
        mu_e(boost::extents[Nr]) //electron chemical potential in MeV
    {
        quantization = quantization_In;
        PolarAxisQuantization = PolarAxisQuantization_In;
        if(quantization == false){
            //resize magnetization and its derivatives to zero if unused
            J_B.resize(boost::extents[0][0][0]);
            M.resize(boost::extents[0][0]);
            //chi.resize(boost::extents[0][0]); don't resize this because it is always used, even if zero, to compute CFL-limited timestep
            M_mu.resize(boost::extents[0][0]);
            C_m.resize(boost::extents[0][0]);
        }
    }
};

/*
    Functions like Numpy's linspace
 */
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

/*
    Signum function- returns +1, -1, or zero based on sign of argument x
 */
template <typename T>
double sgn(T x, T eps = static_cast<T>(1e-12))
{
    static_assert(std::is_floating_point_v<T>, "T must be floating point");

    if (std::abs(x) <= eps) return 0.;
    return (x > T(0)) ? 1. : -1.;
}

template <typename T>
double sgnNonZero(T x, T eps = static_cast<T>(1e-12))
{
    static_assert(std::is_floating_point_v<T>, "T must be floating point");

    return (x > T(0)) ? 1. : -1.;
}

void load_params(SimParams & params, BParams & Bparams, TParams & Tparams, ConductParams & Cparams, int world_rank);
void MPE_Decomp1D(size_t N, int num_procs, int MyID, size_t & s, size_t & e);
void exchng2Vector(VectorField & A, size_t N_GC, const Process & process, MPI_Datatype stridetype);
void exchng2Scalar(ScalarField & A, size_t N_GC, const Process & process, MPI_Datatype stridetype);
void exchng2Array(std::vector<double> & A, size_t N_GC, const Process & process);
void save_timesCalc(std::vector<double> & t_next_Vec, double t_max, size_t saves_number);
// double MC(double a, double b, double c);
double superbee(double a, double b);
double vanLeer(double a, double b);
// double minmod2(double a, double b);
double minmodFlux(double r);
double superbeeFlux(double r);
double OuterEdgeLimiter(double r, double r_max, double cl);
double EdgeLimiter(double r, double r_max, double r_min, double cl);
double PoleEdgeLimiter(double theta, double cl);
double TrapezoidIntegrator(std::vector<double> & A, std::vector<double> & dx);
double dot(VectorField & A, VectorField & B, size_t i, size_t j);
double rdot(VectorField & A, VectorField & B, size_t i, size_t j);
double thdot(VectorField & A, VectorField & B, size_t i, size_t j);
double phdot(VectorField & A, VectorField & B, size_t i, size_t j);
double gradr(ScalarField & A, double Deltar, size_t i, size_t j);
double gradth(ScalarField & A, double rDeltatheta, size_t i, size_t j);
double gradrCentered(ScalarField & A, double Deltar, size_t i, size_t j);
double gradthCentered(ScalarField & A, double rDeltatheta, size_t i, size_t j);
void Bmag_and_BhatCalc(VectorField & B, MagCoeffs & mC, ScalarField & eta_H, ScalarField & eta_O, double & max_eta_local, const Domain & domain, size_t N_GC, const Process & process);

#endif // COMMON_H_INCLUDED
