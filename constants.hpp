#ifndef CONST_H
#define CONST_H

#include <map>
#include <string>

namespace CONST {

extern "C" {
extern std::map<std::string, double> element_masses;

extern const double PI;
extern const double PI_div_2;
extern const double PI_mul_2;
extern const double SQRT2PI;
extern const double LN10;

extern const double phy_elementaryCharge_SI;
extern const double phy_electronClassicalRadius_SI;
extern const double phy_electronClassicalRadius_CGS;
extern const double phy_CoulombConst_SI;
extern const double phy_atomicMassUnit_SI ;
extern const double phy_mProton_SI ;
extern const double phy_mProton_CGS;
extern const double phy_mElectron_SI ;
extern const double phy_mElectron_CGS;
extern const double phy_kBoltzmann_SI ;
extern const double phy_kBoltzmann_CGS;
extern const double phy_hPlanck_SI ;
extern const double phy_hPlanck_CGS;
extern const double phy_hbarPlanck_SI ;
extern const double phy_hbarPlanck_CGS;
extern const double phy_GravitationConst_SI ;
extern const double phy_GravitationConst_CGS;
extern const double phy_SpeedOfLight_SI ;
extern const double phy_SpeedOfLight_CGS;
extern const double phy_StefanBoltzmann_SI;
extern const double phy_StefanBoltzmann_CGS;
extern const double phy_IdealGasConst_SI;
extern const double phy_ThomsonScatterCross_CGS;

extern const double phy_Lsun_SI;
extern const double phy_Lsun_CGS;
extern const double phy_Msun_SI;
extern const double phy_Msun_CGS;
extern const double phy_Rsun_SI;
extern const double phy_Rsun_CGS;
extern const double phy_Rearth_CGS;
extern const double phy_Mearth_CGS;
extern const double phy_Rmoon_CGS;
extern const double phy_Mmoon_CGS;
extern const double phy_RJupiter_CGS;
extern const double phy_MJupiter_CGS;

extern const double phy_SecondsPerYear;
extern const double phy_Deg2Rad;
extern const double phy_erg2joule;
extern const double phy_m2cm;
extern const double phy_kg2g;
extern const double phy_eV2erg;
extern const double phy_cm_1_2erg;
extern const double phy_cm_1_2K;
extern const double phy_AvogadroConst;
extern const double phy_AU2cm;
extern const double phy_AU2m ;
extern const double phy_pc2m ;
extern const double phy_pc2cm;
extern const double phy_Angstrom2micron ;
extern const double phy_Angstrom2cm ;
extern const double phy_micron2cm ;

extern const double phy_jansky2CGS;
extern const double phy_jansky2SI ;

extern const double phy_CMB_T;

extern const double phy_ratioDust2GasMass_ISM;
extern const double phy_Habing_photon_energy_CGS;
extern const double phy_LyAlpha_energy_CGS;
extern const double phy_UV_cont_energy_CGS;
extern const double phy_Habing_energy_density_CGS;
extern const double phy_Habing_photon_flux_CGS;
extern const double phy_Habing_energy_flux_CGS;
extern const double phy_UVext2Av;

extern const double phy_LyAlpha_nu0;
extern const double phy_LyAlpha_l0;
extern const double phy_LyAlpha_dnul;
extern const double phy_LyAlpha_f12;

extern const double phy_LyAlpha_cross_H2O;
extern const double phy_LyAlpha_cross_OH ;

extern const double phy_cosmicray_ionization_rate;
extern const double phy_cosmicray_desorption_factor;
extern const double phy_cosmicray_desorption_T;
extern const double phy_cosmicray_attenuate_N_CGS;
extern const double phy_cosmicray_attenuate_m_CGS;
extern const double phy_PAH_abundance_0;
extern const double phy_SitesDensity_CGS;
extern const double phy_DiffBarrierWidth_CGS;
extern const double phy_Diff2DesorRatio;
extern const double phy_DiffBarrierDefault;
extern const double phy_vibFreqDefault;

extern const double colDen2Av_coeff;
extern const double phy_colDen2Av_coeff;

}

}

#endif //CONST_H
