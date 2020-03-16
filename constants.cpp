#include <map>
#include <string>
//#include "constants.hpp"

//using namespace CONST;

extern "C" {
std::map<std::string, double> element_masses = {
  {"+",     0.0},
  {"-",     0.0},
  {"Grain", 0.0},
  {"E",     5.45e-4},
  {"H",     1.0},
  {"D",     2.0},
  {"He",    4.0},
  {"Li",    6.9},
  {"C",    12.0},
  {"N",    14.0},
  {"O",    16.0},
  {"F",    19.0},
  {"Na",   23.0},
  {"Mg",   24.3},
  {"Al",   27.0},
  {"Si",   28.1},
  {"P",    31.0},
  {"S",    32.1},
  {"Cl",   35.5},
  {"Ar",   39.9},
  {"K",    39.1},
  {"Ca",   40.1},
  {"Fe",   55.8}
};


double PI = 3.1415926535897932384626433;
double PI_div_2 = 1.57079632679489661923132;
double PI_mul_2 = 6.283185307179586476925;
double SQRT2PI = 2.5066282746310005024;
double LN10 = 2.30258509299;

double phy_elementaryCharge_SI = 1.602176487e-19;
double phy_electronClassicalRadius_SI = 2.8179403267e-15;
double phy_electronClassicalRadius_CGS = 2.8179403267e-13;
double phy_CoulombConst_SI = 8.9875517873681764e9;
double phy_atomicMassUnit_SI  = 1.660539040e-27;
double phy_mProton_SI  = 1.67262158e-27;
double phy_mProton_CGS = 1.67262158e-24;
double phy_mElectron_SI  = 9.10938188E-31;
double phy_mElectron_CGS = 9.10938188e-28;
double phy_kBoltzmann_SI  = 1.3806503e-23;
double phy_kBoltzmann_CGS = 1.3806503e-16;
double phy_hPlanck_SI  = 6.62606896e-34;
double phy_hPlanck_CGS = 6.62606896e-27;
double phy_hbarPlanck_SI  = 1.054571628e-34;
double phy_hbarPlanck_CGS = 1.054571628e-27;
double phy_GravitationConst_SI  = 6.67428e-11;
double phy_GravitationConst_CGS = 6.67428e-8;
double phy_SpeedOfLight_SI  = 299792458e0;
double phy_SpeedOfLight_CGS = 299792458e2;
double phy_StefanBoltzmann_SI = 5.6704e-8;
double phy_StefanBoltzmann_CGS = 5.670373e-5;
double phy_IdealGasConst_SI = 8.314472e0;
double phy_ThomsonScatterCross_CGS = 6.6524574e-25;

double phy_Lsun_SI = 3.839e26; // J s-1
double phy_Lsun_CGS = 3.839e33; // erg s-1
double phy_Msun_SI = 1.9891e30; // kg
double phy_Msun_CGS = 1.9891e33; // g
double phy_Rsun_SI = 6.955e8; // m
double phy_Rsun_CGS = 6.955e10; // cm
double phy_Rearth_CGS = 6371e5; // cm
double phy_Mearth_CGS = 5.97219e27; // g
double phy_Rmoon_CGS = 1737.1e5; // cm
double phy_Mmoon_CGS = 7.34767309e25; // g
double phy_RJupiter_CGS = 7.1492e9; // cm
double phy_MJupiter_CGS = 1.8982e30; // g

double phy_SecondsPerYear = 3600e0*24e0*365e0;
double phy_Deg2Rad = PI/180e0;
double phy_erg2joule = 1e-7;
double phy_m2cm = 1e2;
double phy_kg2g = 1e3;
double phy_eV2erg = 1.60217657e-12;
double phy_cm_1_2erg = phy_hPlanck_CGS * phy_SpeedOfLight_CGS;
double phy_cm_1_2K = phy_cm_1_2erg/phy_kBoltzmann_CGS;
double phy_AvogadroConst = 6.02214179e23;
double phy_AU2cm = 1.49597871e13;
double phy_AU2m  = 1.49597871e11;
double phy_pc2m  = 3.08567758e16;
double phy_pc2cm = 3.08567758e18;
double phy_Angstrom2micron  = 1e-4;
double phy_Angstrom2cm  = 1e-8;
double phy_micron2cm  = 1e-4;

double phy_jansky2CGS = 1e-23;
double phy_jansky2SI  = 1e-26;

double phy_CMB_T = 2.72548e0;

double phy_ratioDust2GasMass_ISM = 0.01e0;
double phy_Habing_photon_energy_CGS = 1.99e-11;
double phy_LyAlpha_energy_CGS = 1.64e-11;
double phy_UV_cont_energy_CGS = phy_Habing_photon_energy_CGS;
double phy_Habing_energy_density_CGS = 5.29e-14; // Draine 2011 book, equation 12.6
double phy_Habing_photon_flux_CGS = 6e7; // cm-2 s-1
double phy_Habing_energy_flux_CGS = 1.194e-3; // erg cm-2 s-1
double phy_UVext2Av = 2.6e0; // Tielens 2005, eq 3.19

double phy_LyAlpha_nu0 = 2.4660718e15;
double phy_LyAlpha_l0 = 1215.668e0;
double phy_LyAlpha_dnul = 9.938e7;
double phy_LyAlpha_f12 = 0.4162e0;

double phy_LyAlpha_cross_H2O = 1.2e-17; // Van Dishoeck 2006, Table 1
double phy_LyAlpha_cross_OH  = 1.8e-18; // Van Dishoeck 2006, Table 1

double phy_cosmicray_ionization_rate = 1.36e-17; // s-1
double phy_cosmicray_desorption_factor = 3.1e-19; // s-1
double phy_cosmicray_desorption_T = 70.0; // s-1
double phy_cosmicray_attenuate_N_CGS = 5.75e25; // 96 g cm-2, Nomura 2007
double phy_cosmicray_attenuate_m_CGS = 96.0; // 96 g cm-2, Nomura 2007
double phy_PAH_abundance_0 = 1.6e-7;
double phy_SitesDensity_CGS = 1e15;
double phy_DiffBarrierWidth_CGS = 1.0e-8;
double phy_Diff2DesorRatio = 0.5;
double phy_DiffBarrierDefault = 1e4;
double phy_vibFreqDefault = 1e12;

double colDen2Av_coeff = 1e-21; // Sun Kwok, eq 10.21
double phy_colDen2Av_coeff = 5.3e-22; // Draine 2011, eq 21.7
}
