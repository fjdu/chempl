#include "math.h"
#include "types.hpp"
#include "constants.hpp"

#define NV_Ith_S(x, i) (x[i])

namespace CALC_RATE {


inline TYPES::DTP_FLOAT thermal_velocity_CGS(
    const TYPES::DTP_FLOAT T_CGS,
    const TYPES::DTP_FLOAT massnum) {
  return sqrt((8.0/CONST::PI)
              * CONST::phy_kBoltzmann_CGS * T_CGS
              / (massnum * CONST::phy_mProton_CGS));
}


TYPES::DTP_FLOAT rate_adsorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = r.abc[0] // Act as sticking efficiency
       * p.n_gas
       * p.dust2gas_num
       * p.dust_crosssec
       * thermal_velocity_CGS(p.T_gas,
                              s.massSpecies.at(r.idxReactants[0]));
  r.rate = NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_desorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = s.vibFreqs.at(r.idxReactants[0])
    * r.abc[0]
    * (exp(-r.abc[2]/p.T_dust)
     + p.chi_cosmicray
     * CONST::phy_cosmicray_desorption_factor
     * exp(-r.abc[2] / CONST::phy_cosmicray_desorption_T));
  r.rate = r.drdy[0] * NV_Ith_S(y, r.idxReactants[0]);
  return r.rate;
}


double arrhenius(const double &T,
    const std::vector<double> &abc, const int &iS) {
  return abc[iS] * pow(T/3e2, abc[iS+1]) * exp(-abc[iS+2] / T);
}


inline double interval(const double t, const double t0, const double t1,
    const double width0=0.1, const double width1=0.1,
    const double expOverflow=50.0) {
    double x0 = std::min((t0-t) / width0, expOverflow);
    double x1 = std::min((t-t1) / width1, expOverflow);
    double y0 = exp(x0);
    double y1 = exp(x1);
    return 1.0 / ((1.0+y0)*(1.0+y1));
}


TYPES::DTP_FLOAT rateArrhenius(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double T = p.T_gas, tmp, dfTmin = 1e20, dfTmin2 = 1e20;
  int iIntvmin = 0, iIntvmin2, iTrgmin = 0, iTrgmin2 = 0;
  const double TTol = 5.0;

  for (std::size_t i=0; i<r.Trange.size(); ++i) {
    double dfthis = abs(T - r.Trange[i]);
    if (dfTmin > dfthis) {
      dfTmin2 = dfTmin;
      iTrgmin2 = iTrgmin;
      dfTmin = dfthis;
      iTrgmin = i;
    } else if (dfTmin2 > dfthis) {
      dfTmin2 = dfthis;
      iTrgmin2 = i;
    }
  }

  iIntvmin = iTrgmin / 2;
  iIntvmin2 = iTrgmin2 / 2;
  if ((iIntvmin == iIntvmin2) ||
      (dfTmin > TTol)) {
    if ((r.Trange[iIntvmin*2] <= T) && (T <= r.Trange[iIntvmin*2+1])) {
      tmp = p.n_gas * arrhenius(T, r.abc, iIntvmin*3);
    } else {
      tmp = p.n_gas * arrhenius(T, r.abc, iIntvmin2*3);
    }
  } else {
    double dT = 0.5 * abs(r.Trange[iTrgmin] - r.Trange[iTrgmin2]);
    double w = 0.1 + dT, f1, f2;
    if (r.Trange[iIntvmin*2+1] <= r.Trange[iIntvmin2*2]) {
      f1 = interval(T, r.Trange[iIntvmin*2], r.Trange[iIntvmin*2+1]+dT, w, w);
      f2 = interval(T, r.Trange[iIntvmin2*2]-dT, r.Trange[iIntvmin2*2+1], w, w);
    } else {
      f1 = interval(T, r.Trange[iIntvmin*2]-dT, r.Trange[iIntvmin*2+1], w, w);
      f2 = interval(T, r.Trange[iIntvmin2*2], r.Trange[iIntvmin2*2+1]+dT, w, w);
    }
    tmp = p.n_gas * (arrhenius(T, r.abc, iIntvmin*3) * f1 +
                     arrhenius(T, r.abc, iIntvmin2*3) * f2);
  }
  r.drdy[0] = tmp * NV_Ith_S(y, r.idxReactants[1]);
  r.drdy[1] = tmp * NV_Ith_S(y, r.idxReactants[0]);
  r.rate = NV_Ith_S(y, r.idxReactants[0])
       * NV_Ith_S(y, r.idxReactants[1])
       * tmp;
  return r.rate;
}


TYPES::DTP_FLOAT rate_cosmicray_ionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = r.abc[0] * p.chi_cosmicray;
  r.rate = NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_cosmicray_induced_ionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = r.abc[0] * p.chi_cosmicray
       * pow(p.T_gas/3e2, r.abc[1])
       * r.abc[2] / (1.0 - p.dust_albedo);
  r.rate = NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_photoionization(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = r.abc[0] * p.G0_UV * exp(-r.abc[2] * p.Av);
  r.rate = NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_photodissociation_H2(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{ // Draine, 31.25 -- 31.27
  double x = p.Ncol_H2/5e14;
  double t1 = 1.0 + x/p.dv_km_s, t2 = sqrt(1.0+x);
  double f_SS;
  if (s.idx2name[r.idxReactants[0]] == "H2") {
    f_SS = 0.965/(t1*t1)
         + 0.035/t2 * exp(-8.5e-4*t2);
  } else {f_SS = 1.0;}
  r.drdy[0] = f_SS * r.abc[0] * p.G0_UV * exp(-r.abc[2] * p.Av);
  r.rate = NV_Ith_S(y, r.idxReactants[0]) * r.drdy[0];
  return r.rate;
}


inline TYPES::DTP_FLOAT calc_cross_surf_barrier_prob(
    const double Tdust, TYPES::Reaction& r)
{
  TYPES::DTP_FLOAT tmp = 1.0;
  if ((r.abc[2] > 0.0) && (r.Trange[0] > 0.0)) {
    tmp =
      exp(std::max(-r.abc[2]/Tdust,
                   -2.0 * r.abc[1] * CONST::phy_DiffBarrierWidth_CGS
                   / CONST::phy_hbarPlanck_CGS
                   * sqrt(2.0 * r.Trange[0] * CONST::phy_mProton_CGS
                        * CONST::phy_kBoltzmann_CGS * r.abc[2])));
  }
  return tmp;
}


TYPES::DTP_FLOAT rate_surface_AA(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double Tdust = p.T_dust;
  int i0 = r.idxReactants[0];
  TYPES::DTP_FLOAT tmp =
         std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/Tdust),
                  s.quantMobilities.at(i0))
       / (4.0*p.dust_crosssec*p.dust_site_density)
       / p.dust2gas_num * r.abc[0] * calc_cross_surf_barrier_prob(Tdust, r);
  r.drdy[0] = tmp * NV_Ith_S(y, i0);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  r.rate = NV_Ith_S(y, i0) * NV_Ith_S(y, i0) * tmp;
  return r.rate;
}


TYPES::DTP_FLOAT rate_surface_AB(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double Tdust = p.T_dust;
  int i0 = r.idxReactants[0], i1 = r.idxReactants[1];
  TYPES::DTP_FLOAT tmp =
      (std::max(s.vibFreqs.at(i0) * exp(-s.diffBarriers.at(i0)/Tdust),
                s.quantMobilities.at(i0))
     + std::max(s.vibFreqs.at(i1) * exp(-s.diffBarriers.at(i1)/Tdust),
                s.quantMobilities.at(i1)))
    / (4.0*p.dust_crosssec*p.dust_site_density)
    / p.dust2gas_num * r.abc[0] * calc_cross_surf_barrier_prob(Tdust, r);
  r.drdy[0] = tmp * NV_Ith_S(y, i1);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  r.rate = NV_Ith_S(y, i0) * NV_Ith_S(y, i1) * tmp;
  return r.rate;
}


TYPES::DTP_FLOAT rate_surface_AA_desorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double f = p.chemdesorption_factor;
  rate_surface_AA(t, y, r, p, s, m);
  r.drdy[0] *= f;
  r.drdy[1] *= f;
  r.rate *= f;
  return r.rate;
}


TYPES::DTP_FLOAT rate_surface_AB_desorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double f = p.chemdesorption_factor;
  rate_surface_AB(t, y, r, p, s, m);
  r.drdy[0] *= f;
  r.drdy[1] *= f;
  r.rate *= f;
  return r.rate;
}


void update_surfmant(
    const TYPES::DTP_FLOAT& t,
    double *y,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m) {
  m.t_calc = t;
  m.k_ads_tot = 0.0;
  m.k_eva_tot = 0.0;
  m.surf_tot = 0.0;
  m.mant_tot = 0.0;
  for (auto & a: m.ads_reactions) {
    m.k_ads_tot += rate_adsorption(t, y, a, p, s, m);
  }
  for (auto & a: m.eva_reactions) {
    m.k_eva_tot += rate_desorption(t, y, a, p, s, m);
  }
  for (auto const& i: s.surfaceSpecies) {
    m.surf_tot += NV_Ith_S(y, i);
  }
  for (auto const& i: s.mantleSpecies) {
    m.mant_tot += NV_Ith_S(y, i);
  }
}


TYPES::DTP_FLOAT rate_surf2mant(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    update_surfmant(t, y, p, s, m);
  }
  r.drdy[0] = m.k_ads_tot / p.dust2gas_num
       / (4.0*p.dust_crosssec*p.dust_site_density);
  r.rate = NV_Ith_S(y, i0) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_mant2surf(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  int i0 = r.idxReactants[0];
  if (m.t_calc != t) {
    update_surfmant(t, y, p, s, m);
  }
  r.drdy[0] = m.k_eva_tot
       / (std::max(m.surf_tot, m.mant_tot) + 1e-60);
  r.rate = NV_Ith_S(y, i0) * r.drdy[0];
  return r.rate;
}


TYPES::DTP_FLOAT rate_iongrain(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double TemperatureReduced = CONST::phy_kBoltzmann_SI * p.T_gas
                            / (pow(CONST::phy_elementaryCharge_SI,2)
                             * CONST::phy_CoulombConst_SI
                             / (p.dust_radius*1e-2));
  double JNegaPosi = (1.0 + 1.0/TemperatureReduced)
                   * (1.0 + sqrt(2.0/(2.0+TemperatureReduced)));
  double JChargeNeut = 1.0 + sqrt(CONST::PI/2.0/TemperatureReduced);

  int i0 = r.idxReactants[0], i1 = r.idxReactants[1], i2;
  int charge = s.elementsSpecies.at(i0).at("+") * s.elementsSpecies.at(i1).at("-") +
               s.elementsSpecies.at(i0).at("-") * s.elementsSpecies.at(i1).at("+");
  if (s.elementsSpecies.at(i0).at("Grain") == 0) {
    i2 = i0;
  } else if (s.elementsSpecies.at(i1).at("Grain") == 0) {
    i2 = i1;
  } else {
    std::cout << "Error in rate_iongrain." << std::endl;
  }
  int tmp;
  if (charge < 0) {
    tmp = p.dust_crosssec * JNegaPosi
        * thermal_velocity_CGS(p.T_gas, s.massSpecies.at(i2));
  } else {
    tmp = p.dust_crosssec * JChargeNeut
        * thermal_velocity_CGS(p.T_gas, s.massSpecies.at(i2));
  }
  r.drdy[0] = tmp * NV_Ith_S(y, i1);
  r.drdy[1] = tmp * NV_Ith_S(y, i0);
  r.rate = NV_Ith_S(y, i0) * NV_Ith_S(y, i1) * tmp;
  return r.rate;
}


TYPES::DTP_FLOAT rate_photodesorption(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  double yield = r.abc[0] + r.abc[1] * p.T_gas;
  double uv_flux = p.G0_UV * exp(-CONST::phy_UVext2Av * p.Av)
                 * CONST::phy_Habing_photon_flux_CGS;
  double r0 = uv_flux * yield * p.dust2gas_num * p.dust_crosssec;
  double tmp1 = p.dust2gas_num
              * (4.0 * p.dust_site_density * p.dust_crosssec);
  double tmp = NV_Ith_S(y, r.idxReactants[0]) / tmp1;
  if (tmp < 1e-6) {
    r.drdy[0] = r0 / tmp1;
    r.rate = r0 * tmp;
  } else {
    double tmp2 = exp(-tmp);
    r.drdy[0] = r0 * tmp2 / tmp1;
    r.rate = r0 * (1.0 - tmp2);
  }
  return r.rate;
}


TYPES::DTP_FLOAT rate_CO_photodissociation(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{ // Adapted from: udfa_rate12/rate13_cse_code/cssubs.f
  // This function is for the case of expanding envelope of evolved stars.
  // It should NOT be used in the generic ISM cases.
  const double FRACE = 1.0/3.0, LAMDAE = 1000.0*1.0E-08,
    FOSCE = 0.017, BANDS = 1.0, smallnum=1e-6;
  double XCO = y[r.idxReactants[0]];
  double AUV = p.Ncol_H * CONST::phy_colDen2AUV_1000A;
  double GAMMAD = exp(-1.644*pow(AUV, 0.86));
  double c0 = 1.5 * 0.0265 * FRACE * FOSCE * LAMDAE
            / (p.v_km_s*1e5) * p.Ncol_H;
  double c1 = r.abc[0] * BANDS * GAMMAD;
  double TAUE = c0 * XCO;
  double BETAE, dBdX;
  if (abs(TAUE) >= smallnum) {
    BETAE = (1.0 - exp(-TAUE)) / TAUE;
    dBdX = (-BETAE + (1.0-BETAE)/TAUE) * c0;
  } else {
    BETAE = 1.0 - TAUE * 0.5;
    dBdX = -0.5 * c0;
  }
  r.rate = c1 * BETAE * XCO;
  r.drdy[0] = c1 * (BETAE + XCO * dBdX);
  return r.rate;
}


TYPES::DTP_FLOAT rate_CO_photodissociation_better(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  // This function can be used for generic ISM.
  // Reference: Morris & Jura 1983
  const double FRACE = 1.0/3.0, LAMDAE = 1000.0*1.0E-08,
    FOSCE = 0.017, BANDS = 1.0, smallnum=1e-6;
  double AUV = p.Ncol_H * CONST::phy_colDen2AUV_1000A;
  double GAMMAD = exp(-1.644*pow(AUV, 0.86));
  double c0 = 1.5 * 0.0265 * FRACE * FOSCE * LAMDAE
            / (p.dv_km_s*1e5) * p.Ncol_CO;
  double c1 = p.G0_UV * r.abc[0] * BANDS * GAMMAD;
  double TAUE = c0;
  double BETAE;
  if (abs(TAUE) >= smallnum) {
    BETAE = (1.0 - exp(-TAUE)) / TAUE;
  } else {
    BETAE = 1.0 - TAUE * 0.5;
  }
  r.drdy[0] = c1 * BETAE;
  r.rate = c1 * BETAE * y[r.idxReactants[0]];
  return r.rate;
}


TYPES::DTP_FLOAT rate_dummy(
    const TYPES::DTP_FLOAT& t,
    double *y,
    TYPES::Reaction& r,
    const TYPES::PhyParams& p,
    const TYPES::Species& s,
    TYPES::AuxData& m)
{
  r.drdy[0] = 0.0;
  r.drdy[1] = 0.0;
  r.rate = 0.0;
  return r.rate;
}


void assignAReactionHandler(TYPES::RateCalculators& rcs,
                            const TYPES::RateCalculator& rc,
                            const int& itype) {
  rcs[itype] = rc;
}


void assignReactionHandlers(TYPES::Chem_data& user_data) {
  (user_data.rate_calculators)[1]  = rate_cosmicray_ionization;
  (user_data.rate_calculators)[71] = rate_cosmicray_ionization;
  (user_data.rate_calculators)[2]  = rate_cosmicray_induced_ionization;
  (user_data.rate_calculators)[72] = rate_cosmicray_induced_ionization;
  (user_data.rate_calculators)[4]  = rate_photodissociation_H2;
  (user_data.rate_calculators)[5]  = rateArrhenius;
  (user_data.rate_calculators)[61] = rate_adsorption;
  (user_data.rate_calculators)[62] = rate_desorption;
  (user_data.rate_calculators)[63] = rate_surface_AA;
  (user_data.rate_calculators)[64] = rate_surface_AB;
  (user_data.rate_calculators)[65] = rate_mant2surf;
  (user_data.rate_calculators)[66] = rate_surf2mant;
  (user_data.rate_calculators)[68] = rate_surface_AA_desorption;
  (user_data.rate_calculators)[69] = rate_surface_AB_desorption;
  (user_data.rate_calculators)[3]  = rate_photoionization;
  (user_data.rate_calculators)[20] = rate_cosmicray_induced_ionization;
  (user_data.rate_calculators)[21] = rate_iongrain;
  (user_data.rate_calculators)[75] = rate_photodesorption;
  (user_data.rate_calculators)[6]  = rate_CO_photodissociation;
  (user_data.rate_calculators)[106]= rate_CO_photodissociation_better;
  (user_data.rate_calculators)[53] = rate_dummy;
  (user_data.rate_calculators)[13] = rate_dummy;
  (user_data.rate_calculators)[67] = rate_dummy;
  //(user_data.rate_calculators)[53] = rateArrhenius;
}

}
