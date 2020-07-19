#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include <utility>
#include "constants.hpp"
#include "types.hpp"
#include "utils.hpp"

namespace TYPES {


TimeDependency::TimeDependency(std::string name_,
      std::vector<double> ts_,
      std::vector<double> vs_) {
      name = name_;
      ts = ts_;
      vs = vs_;
}


inline DTP_FLOAT get_n_gas(PhyParams& p) {return p.n_gas;}
inline DTP_FLOAT get_T_gas(PhyParams& p) {return p.T_gas;}
inline DTP_FLOAT get_T_dust(PhyParams& p) {return p.T_dust;}
inline DTP_FLOAT get_Av(PhyParams& p) {return p.Av;}
inline DTP_FLOAT get_G0_UV(PhyParams& p) {return p.G0_UV;}
inline DTP_FLOAT get_chi_Xray(PhyParams& p) {return p.chi_Xray;}
inline DTP_FLOAT get_chi_cosmicray(PhyParams& p) {return p.chi_cosmicray;}
inline DTP_FLOAT get_dust2gas_num(PhyParams& p) {return p.dust2gas_num;}
inline DTP_FLOAT get_dust2gas_mass(PhyParams& p) {return p.dust2gas_mass;}
inline DTP_FLOAT get_dust_material_density(PhyParams& p) {return p.dust_material_density;}
inline DTP_FLOAT get_dust_site_density(PhyParams& p) {return p.dust_site_density;}
inline DTP_FLOAT get_dust_radius(PhyParams& p) {return p.dust_radius;}
inline DTP_FLOAT get_dust_crosssec(PhyParams& p) {return p.dust_crosssec;}
inline DTP_FLOAT get_dust_albedo(PhyParams& p) {return p.dust_albedo;}
inline DTP_FLOAT get_mean_mol_weight(PhyParams& p) {return p.mean_mol_weight;}
inline DTP_FLOAT get_chemdesorption_factor(PhyParams& p) {return p.chemdesorption_factor;}
inline DTP_FLOAT get_Ncol_CO(PhyParams& p) {return p.Ncol_CO;}
inline DTP_FLOAT get_Ncol_H2(PhyParams& p) {return p.Ncol_H2;}
inline DTP_FLOAT get_Ncol_H(PhyParams& p) {return p.Ncol_H;}
inline DTP_FLOAT get_dv_km_s(PhyParams& p) {return p.dv_km_s;}
inline DTP_FLOAT get_v_km_s(PhyParams& p) {return p.v_km_s;}
inline DTP_FLOAT get_t_max_year(PhyParams& p) {return p.t_max_year;}

void set_n_gas(PhyParams& p, DTP_FLOAT v) {p.n_gas=v;}
void set_T_gas(PhyParams& p, DTP_FLOAT v) {p.T_gas=v;}
void set_T_dust(PhyParams& p, DTP_FLOAT v) {p.T_dust=v;}
void set_Av(PhyParams& p, DTP_FLOAT v) {p.Av=v;}
void set_G0_UV(PhyParams& p, DTP_FLOAT v) {p.G0_UV=v;}
void set_chi_Xray(PhyParams& p, DTP_FLOAT v) {p.chi_Xray=v;}
void set_chi_cosmicray(PhyParams& p, DTP_FLOAT v) {p.chi_cosmicray=v;}
void set_dust2gas_num(PhyParams& p, DTP_FLOAT v) {p.dust2gas_num=v;}
void set_dust2gas_mass(PhyParams& p, DTP_FLOAT v) {p.dust2gas_mass=v;}
void set_dust_material_density(PhyParams& p, DTP_FLOAT v) {p.dust_material_density=v;}
void set_dust_site_density(PhyParams& p, DTP_FLOAT v) {p.dust_site_density=v;}
void set_dust_radius(PhyParams& p, DTP_FLOAT v) {p.dust_radius=v;}
void set_dust_crosssec(PhyParams& p, DTP_FLOAT v) {p.dust_crosssec=v;}
void set_dust_albedo(PhyParams& p, DTP_FLOAT v) {p.dust_albedo=v;}
void set_mean_mol_weight(PhyParams& p, DTP_FLOAT v) {p.mean_mol_weight=v;}
void set_chemdesorption_factor(PhyParams& p, DTP_FLOAT v) {p.chemdesorption_factor=v;}
void set_Ncol_CO(PhyParams& p, DTP_FLOAT v) {p.Ncol_CO=v;}
void set_Ncol_H2(PhyParams& p, DTP_FLOAT v) {p.Ncol_H2=v;}
void set_Ncol_H(PhyParams& p, DTP_FLOAT v) {p.Ncol_H=v;}
void set_dv_km_s(PhyParams& p, DTP_FLOAT v) {p.dv_km_s=v;}
void set_v_km_s(PhyParams& p, DTP_FLOAT v) {p.v_km_s=v;}
void set_t_max_year(PhyParams& p, DTP_FLOAT v) {p.t_max_year=v;}


std::map<std::string, void (*)(PhyParams&, DTP_FLOAT)> phySetterDict = {
  {"n_gas", set_n_gas},
  {"T_gas", set_T_gas},
  {"T_dust", set_T_dust},
  {"Av", set_Av},
  {"G0_UV", set_G0_UV},
  {"chi_Xray", set_chi_Xray},
  {"chi_cosmicray", set_chi_cosmicray},
  {"dust2gas_num", set_dust2gas_num},
  {"dust2gas_mass", set_dust2gas_mass},
  {"dust_material_density", set_dust_material_density},
  {"dust_site_density", set_dust_site_density},
  {"dust_radius", set_dust_radius},
  {"dust_crosssec", set_dust_crosssec},
  {"dust_albedo", set_dust_albedo},
  {"mean_mol_weight", set_mean_mol_weight},
  {"chemdesorption_factor", set_chemdesorption_factor},
  {"Ncol_CO", set_Ncol_CO},
  {"Ncol_H2", set_Ncol_H2},
  {"Ncol_H", set_Ncol_H},
  {"dv_km_s", set_dv_km_s},
  {"v_km_s", set_v_km_s},
  {"t_max_year", set_t_max_year}
};


std::map<std::string, DTP_FLOAT (*)(PhyParams&)> phyGetterDict = {
  {"n_gas", get_n_gas},
  {"T_gas", get_T_gas},
  {"T_dust", get_T_dust},
  {"Av", get_Av},
  {"G0_UV", get_G0_UV},
  {"chi_Xray", get_chi_Xray},
  {"chi_cosmicray", get_chi_cosmicray},
  {"dust2gas_num", get_dust2gas_num},
  {"dust2gas_mass", get_dust2gas_mass},
  {"dust_material_density", get_dust_material_density},
  {"dust_site_density", get_dust_site_density},
  {"dust_radius", get_dust_radius},
  {"dust_crosssec", get_dust_crosssec},
  {"dust_albedo", get_dust_albedo},
  {"mean_mol_weight", get_mean_mol_weight},
  {"chemdesorption_factor", get_chemdesorption_factor},
  {"Ncol_CO", get_Ncol_CO},
  {"Ncol_H2", get_Ncol_H2},
  {"Ncol_H", get_Ncol_H},
  {"dv_km_s", get_dv_km_s},
  {"v_km_s", get_v_km_s},
  {"t_max_year", get_t_max_year}
};


void PhyParams::prep_params() {
  dust2gas_num = dust2gas_mass * (CONST::phy_mProton_CGS * mean_mol_weight)
              / (4.0*CONST::PI/3.0 * dust_material_density *
                 dust_radius * dust_radius * dust_radius);
  dust_crosssec = CONST::PI * dust_radius * dust_radius;
}


int PhyParams::from_file(std::string fname)
{
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s*(=)\s*(\S+))");
  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        continue;
      }
      std::smatch matches;
      if (std::regex_search(line, matches, entry)) {
        std::string var_name = matches.str(1);
        if (phySetterDict.find(var_name) != phySetterDict.end()) {
          DTP_FLOAT var_val = std::stod(matches.str(3));
          phySetterDict[var_name](*this, var_val);
        } else {
          std::cout << "Invalid line: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error in PhyParams::from_file: " << fname << std::endl;
  }
  return 0;
}


void PhyParams::add_a_timedependency(std::string name,
    std::vector<double> ts, std::vector<double> vs) {
  if (phySetterDict.find(name) == phySetterDict.end()) {
    std::cout << "In add_a_timedependency: key does not exist: "
              << name << std::endl;
    return;
  }
  this->timeDependencies.push_back(TimeDependency(name, ts, vs));
}


void PhyParams::remove_a_timedependency(std::string name) {
  int ifound = -1;
  for (int i=0; i < this->timeDependencies.size(); ++i) {
    if (this->timeDependencies[i].name == name) {
      ifound = i;
      break;
    }
  }
  if (ifound >= 0) {
    this->timeDependencies.erase(this->timeDependencies.begin()+ifound);
  }
}


std::vector<std::string> PhyParams::get_timeDependency_names() {
  std::vector<std::string> names;
  for (auto t: this->timeDependencies) {
    names.push_back(t.name);
  }
  return names;
}


AuxData::AuxData():
  t_calc(NAN), k_eva_tot(0.0), k_ads_tot(0.0),
  mant_tot(0.0), surf_tot(0.0),
  n_surf2mant(0), n_mant2surf(0) {}


Reaction::Reaction() {
  drdy[0] = NAN;
  drdy[1] = NAN;
  rate = 0.0;
}


Reaction::Reaction(
  std::vector<std::string> sR,
  std::vector<std::string> sP,
  std::vector<DTP_FLOAT> a_,
  std::vector<DTP_FLOAT> Tr,
  int it) {
  sReactants = sR;
  sProducts = sP;
  abc = a_;
  Trange = Tr;
  itype = it;
  drdy[0] = NAN;
  drdy[1] = NAN;
  rate = 0.0;
}


bool operator== (const Reaction &r1, const Reaction &r2)
{
  return ((r1.sReactants == r2.sReactants) &&
          (r1.sProducts == r2.sProducts) &&
          (r1.abc == r2.abc) &&
          (r1.Trange == r2.Trange) &&
          (r1.itype == r2.itype));
}


void Chem_data::allocate_y() {
  if (y != nullptr) {
    delete [] y;
  }
  y = new DTP_FLOAT[species.idx2name.size()];
}


void Chem_data::allocate_y(int n=0) {
  if (y != nullptr) {
    delete [] y;
  }
  if (n>0) {
    y = new DTP_FLOAT[n];
  } else {
    y = new DTP_FLOAT[species.idx2name.size()];
  }
}


void Chem_data::deallocate_y() {
  if (y != nullptr) {
    delete [] y;
    y = nullptr;
  }
}


Chem_data::Chem_data() {
  ptr = this;
  y = nullptr;
}


Chem_data::~Chem_data() {
  deallocate_y();
}


bool Chem_data::updateIfIsDuplicate(const Reaction& rs) {
  bool isDuplicate = false;
  for (auto &r: reactions) {
    if ((r.sReactants == rs.sReactants) &&
        (r.sProducts == rs.sProducts) &&
        (r.itype == rs.itype)) {
      int iTr;
      for (iTr=r.Trange.size()-2; iTr>=0; iTr -= 2) {
        if (r.Trange[iTr] <= rs.Trange[0]) {
          break;
        }
      }
      r.Trange.insert(r.Trange.begin()+iTr+2, rs.Trange.begin(), rs.Trange.end());
      r.abc.insert(r.abc.begin() + (iTr*3)/2+3, rs.abc.begin(), rs.abc.end());
      isDuplicate = true;
      break;
    }
  }
  return isDuplicate;
}


void Chem_data::add_reaction(Reaction rs) {
  if (updateIfIsDuplicate(rs)) {
    return;
  }

  int nReactants=rs.sReactants.size(), nProducts=rs.sProducts.size();

  std::vector<std::string> stmp;
  stmp.reserve(nReactants+nProducts);
  stmp.insert(stmp.end(), rs.sReactants.begin(), rs.sReactants.end());
  stmp.insert(stmp.end(), rs.sProducts.begin(), rs.sProducts.end());

  for (int i=0; i<nReactants+nProducts; ++i) {
    std::string name = UTILS::trim(stmp[i]);

    if ((name == "PHOTON") || (name == "CRPHOT") ||
        (name == "CRP") || (name.size() == 0)) {
      continue;
    }

    if (species.name2idx.find(name) == species.name2idx.end()) {
      species.name2idx[name] = species.name2idx.size();
      species.idx2name.push_back(name);
    }

    if (i < nReactants) {
      rs.idxReactants.push_back(species.name2idx[name]);
    } else {
      rs.idxProducts.push_back(species.name2idx[name]);
    }
  }

  if (reaction_types.find(rs.itype) == reaction_types.end()) {
    reaction_types[rs.itype] = 1;
  } else {
    reaction_types[rs.itype] += 1;
  }

  reactions.push_back(rs);
}


void Chem_data::modify_reaction(const int& iReact, const std::map<std::string, std::vector<double> > &par) {
  if ((iReact < 0) || (iReact >= reactions.size())) {return;}
  for (auto &p: par) {
    if (p.first == "abc") {
      reactions[iReact].abc = p.second;
    }
    if (p.first == "Trange") {
      reactions[iReact].Trange = p.second;
    }
    if (p.first == "itype") {
      reaction_types[reactions[iReact].itype] -= 1;
      reactions[iReact].itype = (int)p.second[0];
      reaction_types[reactions[iReact].itype] += 1;
    }
  }
}


void Chem_data::find_duplicate_reactions() {
  for (int i=0; i < reactions.size(); ++i) {
    if (reactions[i].abc.size() > 3) {
      if (std::find(dupli.begin(), dupli.end(), i) == dupli.end()) {
        dupli.push_back(i);
      }
    }
  }
}


void Chem_data::clear_reactions() {
  this->reactions.clear();
  this->reaction_types.clear();
  this->auxdata.ads_reactions.clear();
  this->auxdata.eva_reactions.clear();
  this->auxdata.ads_species.clear();
  this->auxdata.eva_species.clear();
  this->auxdata.n_surf2mant = 0;
  this->auxdata.n_mant2surf = 0;
}


void Chem_data::set_phy_param(std::string pname, DTP_FLOAT val) {
  if (phySetterDict.find(pname) != phySetterDict.end()) {
    phySetterDict[pname](this->physical_params, val);
  } else {
    std::cout << "In set_phy_param:" << std::endl;
    std::cout << "Unknown param name: " << pname << std::endl;
  }
}

DTP_FLOAT Chem_data::get_phy_param(std::string pname) {
  if (phyGetterDict.find(pname) != phyGetterDict.end()) {
    return phyGetterDict[pname](this->physical_params);
  } else {
    std::cout << "In get_phy_param:" << std::endl;
    std::cout << "Unknown param name: " << pname << std::endl;
    return NAN;
  }
}

std::map<std::string, DTP_FLOAT> Chem_data::get_all_phy_params() {
  std::map<std::string, DTP_FLOAT> res;
  for (auto const& p: phyGetterDict) {
    res[p.first] = p.second(this->physical_params);
  }
  return res;
}

void Chem_data::assort_reactions()
{
  std::vector<Reaction>().swap(auxdata.ads_reactions);
  std::vector<Reaction>().swap(auxdata.eva_reactions);
  std::vector<int>().swap(this->auxdata.ads_species);
  std::vector<int>().swap(this->auxdata.eva_species);
  this->auxdata.n_surf2mant = 0;
  this->auxdata.n_mant2surf = 0;

  for (auto const& r: this->reactions) {
    // if ((std::find(auxdata.ads_reactions.begin(),
    //                auxdata.ads_reactions.end(), r)
    //             != auxdata.ads_reactions.end()) ||
    //     (std::find(auxdata.eva_reactions.begin(),
    //                auxdata.eva_reactions.end(), r)
    //             != auxdata.eva_reactions.end())) {
    //   return;
    // }

    if (r.itype == 61) {
      //if ((species.idx2name[r.idxReactants[0]] == "H2") ||
      //    (species.idx2name[r.idxReactants[0]] == "H")) {
      //  continue;
      //}
      this->auxdata.ads_reactions.push_back(r);
      this->auxdata.ads_species.push_back(r.idxReactants[0]);
    } else if (r.itype == 62) {
      //if ((species.idx2name[r.idxReactants[0]] == "gH2") ||
      //    (species.idx2name[r.idxReactants[0]] == "gH")) {
      //  continue;
      //}
      this->auxdata.eva_reactions.push_back(r);
      this->auxdata.eva_species.push_back(r.idxReactants[0]);
    } else if (r.itype == 65) {
      this->auxdata.n_mant2surf += 1;
    } else if (r.itype == 66) {
      this->auxdata.n_surf2mant += 1;
    }
  }
}

std::map<std::string, int> assignElementsToOneSpecies(
    const std::string& name, const Elements& elements)
{
  std::regex npat(R"(^(\d+))");

  std::map<std::string, int> eleDict;
  for (auto const& e: elements) {
    eleDict[e.first] = 0;
  }

  for (int iBg=0; iBg<name.size();) {
    bool found = false;
    for (int nlen=name.size()-iBg; nlen>0; --nlen) {
      auto q = elements.find(name.substr(iBg, nlen));
      if (q == elements.end()) {
        continue;
      }

      iBg += q->first.size();

      std::smatch matches;
      std::string subname = name.substr(iBg);
      if (std::regex_search(subname, matches, npat)) {
        try {
          eleDict[q->first] += std::stoi(matches.str(1));
        } catch (...) {
          std::cout << "Exception in assignElementsToOneSpecies: "
                    << " " << matches.str(1) << std::endl;
        }
        iBg += matches.str(1).size();
      } else {
        eleDict[q->first] += 1;
      }
      found = true;
      break;
    }
    if (!found) {
      ++iBg;
    }
  }
  return eleDict;
}


void Chem_data::assignElementsToSpecies(const Elements& elements) {
  for (auto const& s: this->species.name2idx) {
    this->species.elementsSpecies[s.second] =
    assignElementsToOneSpecies(s.first, elements);
  }
}


void Chem_data::assignElementsToSpecies() {
  assignElementsToSpecies(CONST::element_masses);
}


void Chem_data::calculateSpeciesMasses(const Elements& elements) {
  for (auto const& s: this->species.name2idx) {
    this->species.massSpecies[s.second] = 0.0;
    auto const& eleDict = this->species.elementsSpecies[s.second];
    for (auto const& e: eleDict) {
      this->species.massSpecies[s.second] +=
      elements.at(e.first) * ((double)e.second);
    }
  }
}

void Chem_data::calculateSpeciesMasses() {
  calculateSpeciesMasses(CONST::element_masses);
}

void Chem_data::calculateSpeciesVibFreqs() {
  for (auto const& r: this->reactions) {
    if (r.itype == 62) {
      this->species.vibFreqs[r.idxReactants[0]] =
        sqrt(2.0 * CONST::phy_SitesDensity_CGS
           * CONST::phy_kBoltzmann_CGS * r.abc[2]
           / (CONST::PI * CONST::PI)
           / (CONST::phy_mProton_CGS * this->species.massSpecies[r.idxReactants[0]]));
    }
  }
  for (auto const& r: this->reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (this->species.vibFreqs.find(i) == this->species.vibFreqs.end()) {
          this->species.vibFreqs[i] = CONST::phy_vibFreqDefault;
        }
      }
    }
  }
}

void Chem_data::calculateSpeciesDiffBarriers() {
  for (auto const& r: this->reactions) {
    if (r.itype == 62) {
      this->species.diffBarriers[r.idxReactants[0]] =
      r.abc[2] * CONST::phy_Diff2DesorRatio;
    }
  }
  for (auto const& r: this->reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (this->species.diffBarriers.find(i) ==
            this->species.diffBarriers.end()) {
          this->species.diffBarriers[i] = CONST::phy_DiffBarrierDefault;
        }
      }
    }
  }
}


void Chem_data::calculateSpeciesQuantumMobilities() {
  for (auto const& r: this->reactions) {
    if (r.itype == 62) {
      this->species.quantMobilities[r.idxReactants[0]] =
          this->species.vibFreqs[r.idxReactants[0]] *
          exp(-2.0 * CONST::phy_DiffBarrierWidth_CGS
             / CONST::phy_hbarPlanck_CGS
             * sqrt(2.0 * this->species.massSpecies[r.idxReactants[0]]
                  * CONST::phy_mProton_CGS
                  * CONST::phy_kBoltzmann_CGS
                  * this->species.diffBarriers[r.idxReactants[0]]));
    }
  }
  for (auto const& r: this->reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (this->species.quantMobilities.find(i) ==
            this->species.quantMobilities.end()) {
          this->species.quantMobilities[i] = 0.0;
        }
      }
    }
  }
}


void Chem_data::classifySpeciesByPhase() {
  for (auto const& s: this->species.name2idx) {
    if (s.first[0] == 'm') {
      this->species.mantleSpecies.insert(s.second);
    } else if (s.first[0] == 'g') {
      this->species.surfaceSpecies.insert(s.second);
    } else {
      this->species.gasSpecies.insert(s.second);
    }
  }
}


inline void calculateOneReactionHeat(Reaction& r,
    const std::map<int, DTP_FLOAT>& enths) {
  double h = 0.0;
  for (auto const& i: r.idxReactants) {
    if (enths.find(i) != enths.end()) {
      h += enths.at(i);
    } else {
      r.heat = 0.0;
      return;
    }
  }
  for (auto const& i: r.idxProducts) {
    if (enths.find(i) != enths.end()) {
      h -= enths.at(i);
    } else {
      r.heat = 0.0;
      return;
    }
  }
  r.heat = h;
}


void Chem_data::calculateReactionHeat() {
  for (auto& r: this->reactions) {
    calculateOneReactionHeat(r, this->species.enthalpies);
  }
}


double Chem_data::calculate_a_rate(
    double t,
    double *y,
    Reaction& r, bool updatePhyParams) {
  if (updatePhyParams) {
    update_phy_params(t, physical_params);
  }
  return rate_calculators[r.itype](t, y, r,
    physical_params, species, auxdata);
}


std::vector<int> Chem_data::getFormationReactions(int iSpecies) {
  std::vector<int> res;
  for (int i=0; i<reactions.size(); ++i) {
    if (std::find(reactions[i].idxProducts.begin(),
                  reactions[i].idxProducts.end(), iSpecies)
        != reactions[i].idxProducts.end()) {
      res.push_back(i);
    }
  }
  return res;
}


std::vector<int> Chem_data::getDestructionReactions(int iSpecies) {
  std::vector<int> res;
  for (int i=0; i<reactions.size(); ++i) {
    if (std::find(reactions[i].idxReactants.begin(),
                  reactions[i].idxReactants.end(), iSpecies)
        != reactions[i].idxReactants.end()) {
      res.push_back(i);
    }
  }
  return res;
}


std::vector<std::pair<int, double> >
Chem_data::getFormationReactionsWithRates(
  int iSpecies, double t, double* y) {
  update_phy_params(t, physical_params);

  std::vector<std::pair<int, double> > res;
  for (int i=0; i<reactions.size(); ++i) {
    if (std::find(reactions[i].idxProducts.begin(),
                  reactions[i].idxProducts.end(), iSpecies)
        != reactions[i].idxProducts.end()) {
      std::pair<int, double> p;
      p.first = i;
      p.second = calculate_a_rate(t, y, reactions[i]);
      res.push_back(p);
    }
  }
  std::sort(res.begin(), res.end(),
            [](std::pair<int, double> a,
               std::pair<int, double> b) {
                 return a.second > b.second;});
  return res;
}


std::vector<std::pair<int, double> >
Chem_data::getDestructionReactionsWithRates(
  int iSpecies, double t, double* y) {
  update_phy_params(t, physical_params);

  std::vector<std::pair<int, double> > res;
  for (int i=0; i<reactions.size(); ++i) {
    if (std::find(reactions[i].idxReactants.begin(),
                  reactions[i].idxReactants.end(), iSpecies)
        != reactions[i].idxReactants.end()) {
      std::pair<int, double> p;
      p.first = i;
      p.second = calculate_a_rate(t, y, reactions[i]);
      res.push_back(p);
    }
  }
  std::sort(res.begin(), res.end(),
            [](std::pair<int, double> a,
               std::pair<int, double> b) {
                 return a.second > b.second;});
  return res;
}


void update_phy_params(
    TYPES::DTP_FLOAT t,
    TYPES::PhyParams& p) {
  for (auto& s: p.timeDependencies) {
    TYPES::phySetterDict[s.name](p, UTILS::interpol(s.ts, s.vs, t));
  }
}


void Species::allocate_abundances() {
  this->abundances = std::vector<TYPES::DTP_FLOAT>(this->name2idx.size());
}


Recorder::Recorder(std::string fname_): fname(fname_) {
  ofs.open(fname);
  if (ofs.fail()) {
    std::cerr << "Fail to open file: "
              << fname << std::endl;
    throw std::runtime_error(std::strerror(errno));
  }
}


Recorder::~Recorder() {
  ofs.close();
}


int Recorder::write_header(std::vector<std::string> names,
                     int fwidth) {
  ofs << std::setw(fwidth) << std::left << "Time";
  ofs << std::setw(fwidth) << std::left << "T_gas";
  ofs << std::setw(fwidth) << std::left << "T_dust";
  ofs << std::setw(fwidth) << std::left << "n_gas";
  for (std::string const& n: names) {
    ofs << std::setw(fwidth) << std::left << n;
  }
  ofs << std::endl;
  return 0;
}


int Recorder::write_row(double t, int neq, double *y,
              const TYPES::PhyParams& p,
              int fwidth, int prec) {
  ofs << std::setw(fwidth) << std::left << std::scientific
      << std::setprecision(prec) << t / CONST::phy_SecondsPerYear;
  ofs << std::setw(fwidth) << std::left << std::scientific
    << std::setprecision(prec) << p.T_gas;
  ofs << std::setw(fwidth) << std::left << std::scientific
    << std::setprecision(prec) << p.T_dust;
  ofs << std::setw(fwidth) << std::left << std::scientific
    << std::setprecision(prec) << p.n_gas;
  for (int i=0; i<neq; ++i) {
    ofs << std::setw(fwidth) << std::left << std::scientific
      << std::setprecision(prec) << y[i];
  }
  ofs << std::endl;
  return 0;
}


PathsDict loadPathConfig(std::string fname) {
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s*=\s*([^!#]*[^ !#]+)\s*(?:([!#])|($)))");
  std::ifstream inputFile(fname);
  PathsDict pd;
  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        continue;
      }
      std::smatch matches;
      if (std::regex_search(line, matches, entry)) {
        pd[matches.str(1)] = matches.str(2);
      }
    }
  } else {
    std::cout << "Error in loadPathConfig: " << fname << std::endl;
  }
  return pd;
}


Reaction str2reaction(const std::string& str,
    int nReactants, int nProducts, int nABC,
    int lenSpeciesName, int lenABC, int nT, int lenT,
    int lenType, int rowlen_min) {
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  TYPES::Reaction reaction;

  if (std::regex_match(str, comment)) {
    return reaction;
  }
  if (std::regex_match(str, emptyline)) {
    return reaction;
  }
  if (str.size() < rowlen_min) {
    return reaction;
  }

  for (int i=0; i<nReactants+nProducts; ++i) {
    std::string name = str.substr(i*lenSpeciesName, lenSpeciesName);
    name = UTILS::trim(name);

    if ((name == "PHOTON") || (name == "CRPHOT") || (name == "CRP")) {
      continue;
    }

    if (name.size() > 0) {
      if (i < nReactants) {
        reaction.sReactants.push_back(name);
      } else {
        reaction.sProducts.push_back(name);
      }
    }
  }

  int iABC = (nReactants + nProducts) * lenSpeciesName;
  for (int i=0; i<nABC; ++i) {
    std::string tmp = str.substr(iABC+i*lenABC, lenABC);
    tmp = UTILS::trim(tmp);
    if (tmp.size() > 0) {
      reaction.abc.push_back(std::stod(tmp));
    } else {
      reaction.abc.push_back(0.0);
    }
  }

  int iT = iABC + nABC * lenABC;
  for (int i=0; i<nT; ++i) {
    std::string tmp = str.substr(iT+i*lenT, lenT);
    tmp = UTILS::trim(tmp);
    if (tmp.size() > 0) {
      reaction.Trange.push_back(std::stod(tmp));
    } else {
      reaction.Trange.push_back(0.0);
    }
  }

  int iType = iT + nT * lenT;
  std::string tmp = str.substr(iType, lenType);
  try {
    reaction.itype = std::stoi(tmp);
  } catch (...) {
    std::cerr << "Exception in str2reaction: "
              << tmp << std::endl;
  }

  return reaction;
}



void load_reactions(const std::string& fname, Chem_data& user_data,
    int nReactants, int nProducts, int nABC, int lenSpeciesName,
    int lenABC, int nT, int lenT, int lenType, int rowlen_min)
{
  std::ifstream inputFile(fname);
  std::string line;

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {

      Reaction rthis = str2reaction(line,
        nReactants, nProducts, nABC, lenSpeciesName,
        lenABC, nT, lenT, lenType, rowlen_min);
      if ((rthis.sReactants.size() == 0) ||
          (rthis.sProducts.size() == 0)) {
        continue;
      }

      user_data.add_reaction(rthis);
    }
  } else {
    std::cerr << "Error in load_reactions: "
              << fname << std::endl;
  }

  inputFile.close();
}


void loadInitialAbundances(TYPES::Species& species, std::string fname) {
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s+(\S+))");

  if (species.abundances.size() != species.name2idx.size()) {
    species.abundances = std::vector<TYPES::DTP_FLOAT>(species.name2idx.size());
  }
  for (int i=0; i<species.abundances.size(); ++i) {
    species.abundances[i] = 0.0;
  }

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        continue;
      }
      std::smatch matches;
      if (std::regex_search(line, matches, entry)) {
        std::string name = matches.str(1);
        if (species.name2idx.find(name) != species.name2idx.end()) {
          TYPES::DTP_FLOAT val = std::stod(matches.str(2));
          species.abundances[species.name2idx[name]] = val;
        } else {
          std::cout << "In loadInitialAbundances: Unrecognized row: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error in loadInitialAbundances: " << fname << std::endl;
  }
}


void loadSpeciesEnthalpies(TYPES::Species& species, std::string fname) {
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s+(\S+))");

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        continue;
      }
      std::smatch matches;
      if (std::regex_search(line, matches, entry)) {
        std::string name = matches.str(1);
        if (species.name2idx.find(name) != species.name2idx.end()) {
          TYPES::DTP_FLOAT val = std::stod(matches.str(2));
          species.enthalpies[species.name2idx[name]] = val;
        } else {
          //std::cout << "Invalid line: " << line << std::endl;
        }
      }
    }
  } else {
    std::cout << "Error loadSpeciesEnthalpies: " << fname << std::endl;
  }
}


}
