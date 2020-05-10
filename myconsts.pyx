# distutils: language = c++

"""
'<,'>s/\(\<\w\+\>\),\s*/  @property\r  def \1(self):\r    return \1\r\r/g
'<,'>s/\(\<\w\+\>\),\s*/  @property\r  def \1(self):\r    return myconsts.\1\r\r/g
"""

cimport myconsts

class Consts:
  @property
  def element_masses(self):
    return myconsts.element_masses

  @property
  def PI(self):
    return myconsts.PI

  @property
  def PI_div_2(self):
    return myconsts.PI_div_2

  @property
  def PI_mul_2(self):
    return myconsts.PI_mul_2

  @property
  def SQRT2PI(self):
    return myconsts.SQRT2PI

  @property
  def LN10(self):
    return myconsts.LN10

  @property
  def phy_elementaryCharge_SI(self):
    return myconsts.phy_elementaryCharge_SI

  @property
  def phy_electronClassicalRadius_SI(self):
    return myconsts.phy_electronClassicalRadius_SI

  @property
  def phy_electronClassicalRadius_CGS(self):
    return myconsts.phy_electronClassicalRadius_CGS

  @property
  def phy_CoulombConst_SI(self):
    return myconsts.phy_CoulombConst_SI

  @property
  def phy_atomicMassUnit_SI(self):
    return myconsts.phy_atomicMassUnit_SI

  @property
  def phy_mProton_SI(self):
    return myconsts.phy_mProton_SI

  @property
  def phy_mProton_CGS(self):
    return myconsts.phy_mProton_CGS

  @property
  def phy_mElectron_SI(self):
    return myconsts.phy_mElectron_SI

  @property
  def phy_mElectron_CGS(self):
    return myconsts.phy_mElectron_CGS

  @property
  def phy_kBoltzmann_SI(self):
    return myconsts.phy_kBoltzmann_SI

  @property
  def phy_kBoltzmann_CGS(self):
    return myconsts.phy_kBoltzmann_CGS

  @property
  def phy_hPlanck_SI(self):
    return myconsts.phy_hPlanck_SI

  @property
  def phy_hPlanck_CGS(self):
    return myconsts.phy_hPlanck_CGS

  @property
  def phy_hbarPlanck_SI(self):
    return myconsts.phy_hbarPlanck_SI

  @property
  def phy_hbarPlanck_CGS(self):
    return myconsts.phy_hbarPlanck_CGS

  @property
  def phy_GravitationConst_SI(self):
    return myconsts.phy_GravitationConst_SI

  @property
  def phy_GravitationConst_CGS(self):
    return myconsts.phy_GravitationConst_CGS

  @property
  def phy_SpeedOfLight_SI(self):
    return myconsts.phy_SpeedOfLight_SI

  @property
  def phy_SpeedOfLight_CGS(self):
    return myconsts.phy_SpeedOfLight_CGS

  @property
  def phy_StefanBoltzmann_SI(self):
    return myconsts.phy_StefanBoltzmann_SI

  @property
  def phy_StefanBoltzmann_CGS(self):
    return myconsts.phy_StefanBoltzmann_CGS

  @property
  def phy_IdealGasConst_SI(self):
    return myconsts.phy_IdealGasConst_SI

  @property
  def phy_ThomsonScatterCross_CGS(self):
    return myconsts.phy_ThomsonScatterCross_CGS

  @property
  def phy_Lsun_SI(self):
    return myconsts.phy_Lsun_SI

  @property
  def phy_Lsun_CGS(self):
    return myconsts.phy_Lsun_CGS

  @property
  def phy_Msun_SI(self):
    return myconsts.phy_Msun_SI

  @property
  def phy_Msun_CGS(self):
    return myconsts.phy_Msun_CGS

  @property
  def phy_Rsun_SI(self):
    return myconsts.phy_Rsun_SI

  @property
  def phy_Rsun_CGS(self):
    return myconsts.phy_Rsun_CGS

  @property
  def phy_Rearth_CGS(self):
    return myconsts.phy_Rearth_CGS

  @property
  def phy_Mearth_CGS(self):
    return myconsts.phy_Mearth_CGS

  @property
  def phy_Rmoon_CGS(self):
    return myconsts.phy_Rmoon_CGS

  @property
  def phy_Mmoon_CGS(self):
    return myconsts.phy_Mmoon_CGS

  @property
  def phy_RJupiter_CGS(self):
    return myconsts.phy_RJupiter_CGS

  @property
  def phy_MJupiter_CGS(self):
    return myconsts.phy_MJupiter_CGS

  @property
  def phy_SecondsPerYear(self):
    return myconsts.phy_SecondsPerYear

  @property
  def phy_Deg2Rad(self):
    return myconsts.phy_Deg2Rad

  @property
  def phy_erg2joule(self):
    return myconsts.phy_erg2joule

  @property
  def phy_m2cm(self):
    return myconsts.phy_m2cm

  @property
  def phy_kg2g(self):
    return myconsts.phy_kg2g

  @property
  def phy_eV2erg(self):
    return myconsts.phy_eV2erg

  @property
  def phy_cm_1_2erg(self):
    return myconsts.phy_cm_1_2erg

  @property
  def phy_cm_1_2K(self):
    return myconsts.phy_cm_1_2K

  @property
  def phy_AvogadroConst(self):
    return myconsts.phy_AvogadroConst

  @property
  def phy_AU2cm(self):
    return myconsts.phy_AU2cm

  @property
  def phy_AU2m(self):
    return myconsts.phy_AU2m

  @property
  def phy_pc2m(self):
    return myconsts.phy_pc2m

  @property
  def phy_pc2cm(self):
    return myconsts.phy_pc2cm

  @property
  def phy_Angstrom2micron(self):
    return myconsts.phy_Angstrom2micron

  @property
  def phy_Angstrom2cm(self):
    return myconsts.phy_Angstrom2cm

  @property
  def phy_micron2cm(self):
    return myconsts.phy_micron2cm

  @property
  def phy_jansky2CGS(self):
    return myconsts.phy_jansky2CGS

  @property
  def phy_jansky2SI(self):
    return myconsts.phy_jansky2SI

  @property
  def phy_CMB_T(self):
    return myconsts.phy_CMB_T

  @property
  def phy_ratioDust2GasMass_ISM(self):
    return myconsts.phy_ratioDust2GasMass_ISM

  @property
  def phy_Habing_photon_energy_CGS(self):
    return myconsts.phy_Habing_photon_energy_CGS

  @property
  def phy_LyAlpha_energy_CGS(self):
    return myconsts.phy_LyAlpha_energy_CGS

  @property
  def phy_UV_cont_energy_CGS(self):
    return myconsts.phy_UV_cont_energy_CGS

  @property
  def phy_Habing_energy_density_CGS(self):
    return myconsts.phy_Habing_energy_density_CGS

  @property
  def phy_Habing_photon_flux_CGS(self):
    return myconsts.phy_Habing_photon_flux_CGS

  @property
  def phy_Habing_energy_flux_CGS(self):
    return myconsts.phy_Habing_energy_flux_CGS

  @property
  def phy_UVext2Av(self):
    return myconsts.phy_UVext2Av

  @property
  def phy_LyAlpha_nu0(self):
    return myconsts.phy_LyAlpha_nu0

  @property
  def phy_LyAlpha_l0(self):
    return myconsts.phy_LyAlpha_l0

  @property
  def phy_LyAlpha_dnul(self):
    return myconsts.phy_LyAlpha_dnul

  @property
  def phy_LyAlpha_f12(self):
    return myconsts.phy_LyAlpha_f12

  @property
  def phy_LyAlpha_cross_H2O(self):
    return myconsts.phy_LyAlpha_cross_H2O

  @property
  def phy_LyAlpha_cross_OH(self):
    return myconsts.phy_LyAlpha_cross_OH

  @property
  def phy_cosmicray_ionization_rate(self):
    return myconsts.phy_cosmicray_ionization_rate

  @property
  def phy_cosmicray_desorption_factor(self):
    return myconsts.phy_cosmicray_desorption_factor

  @property
  def phy_cosmicray_desorption_T(self):
    return myconsts.phy_cosmicray_desorption_T

  @property
  def phy_cosmicray_attenuate_N_CGS(self):
    return myconsts.phy_cosmicray_attenuate_N_CGS

  @property
  def phy_cosmicray_attenuate_m_CGS(self):
    return myconsts.phy_cosmicray_attenuate_m_CGS

  @property
  def phy_PAH_abundance_0(self):
    return myconsts.phy_PAH_abundance_0

  @property
  def phy_SitesDensity_CGS(self):
    return myconsts.phy_SitesDensity_CGS

  @property
  def phy_DiffBarrierWidth_CGS(self):
    return myconsts.phy_DiffBarrierWidth_CGS

  @property
  def phy_Diff2DesorRatio(self):
    return myconsts.phy_Diff2DesorRatio

  @property
  def phy_DiffBarrierDefault(self):
    return myconsts.phy_DiffBarrierDefault

  @property
  def phy_vibFreqDefault(self):
    return myconsts.phy_vibFreqDefault

  @property
  def colDen2Av_coeff(self):
    return myconsts.colDen2Av_coeff

  @property
  def phy_colDen2Av_coeff(self):
    return myconsts.phy_colDen2Av_coeff

  @property
  def phy_colDen2AUV_1000A(self):
    return myconsts.phy_colDen2AUV_1000A
