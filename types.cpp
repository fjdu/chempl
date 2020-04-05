#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include "constants.hpp"
#include "types.hpp"
#include "utils.hpp"

namespace TYPES {

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
inline DTP_FLOAT get_Ncol_H2(PhyParams& p) {return p.Ncol_H2;}
inline DTP_FLOAT get_dv_km_s(PhyParams& p) {return p.dv_km_s;}
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
void set_Ncol_H2(PhyParams& p, DTP_FLOAT v) {p.Ncol_H2=v;}
void set_dv_km_s(PhyParams& p, DTP_FLOAT v) {p.dv_km_s=v;}
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
  {"Ncol_H2", set_Ncol_H2},
  {"dv_km_s", set_dv_km_s},
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
  {"Ncol_H2", get_Ncol_H2},
  {"dv_km_s", get_dv_km_s},
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


void User_data::add_reaction(const Reaction& rs) {
  Species& species = this->species;
  ReactionTypes& r_types = this->reaction_types;
  Reaction reaction = rs;

  int nReactants=rs.sReactants.size(), nProducts=rs.sProducts.size();

  std::vector<std::string> stmp;
  stmp.reserve(nReactants+nProducts);

  stmp.insert(stmp.end(), rs.sReactants.begin(), rs.sReactants.end());
  stmp.insert(stmp.end(), rs.sProducts.begin(), rs.sProducts.end());

  for (int i=0; i<nReactants+nProducts; ++i) {
    std::string name = UTILS::trim(stmp[i]);

    if ((name == "PHOTON") || (name == "CRPHOT") || (name == "CRP")) {
      continue;
    }

    if ((name.size() > 0) &&
        (species.name2idx.find(name) == species.name2idx.end())) {
      int nSpecies = species.name2idx.size();
      species.name2idx[name] = nSpecies;
      species.idx2name.push_back(name);
    }

    if (name.size() > 0) {
      int iSpecies = species.name2idx[name];
      if (i < nReactants) {
        reaction.idxReactants.push_back(iSpecies);
      } else {
        reaction.idxProducts.push_back(iSpecies);
      }
    }
  }

  if (r_types.find(reaction.itype) == r_types.end()) {
    r_types[reaction.itype] = 1;
  } else {
    r_types[reaction.itype] += 1;
  }

  this->reactions.push_back(reaction);

  // int n=this->reactions.size();
  // if (n > 0) {
  //   std::cout << "Reaction " << n << " added: ";
  //   for (int i=0; i<nReactants; ++i) {
  //     std::cout << this->reactions[n-1].sReactants[i];
  //     if (i<nReactants-1) {
  //       std::cout << " + ";
  //     }
  //   }
  //   std::cout << " -> ";
  //   for (int i=0; i<nProducts; ++i) {
  //     std::cout << this->reactions[n-1].sProducts[i];
  //     if (i<nProducts-1) {
  //       std::cout << " + ";
  //     }
  //   }
  //   std::cout << std::endl;
  // }
}


void User_data::clear_reactions() {
  this->reactions.clear();
  this->reaction_types.clear();
  this->auxdata.ads_reactions.clear();
  this->auxdata.eva_reactions.clear();
}


void User_data::set_phy_param(std::string pname, DTP_FLOAT val) {
  if (phySetterDict.find(pname) != phySetterDict.end()) {
    phySetterDict[pname](this->physical_params, val);
  } else {
    std::cout << "In set_phy_param:" << std::endl;
    std::cout << "Unknown param name: " << pname << std::endl;
  }
}

DTP_FLOAT User_data::get_phy_param(std::string pname) {
  if (phyGetterDict.find(pname) != phyGetterDict.end()) {
    return phyGetterDict[pname](this->physical_params);
  } else {
    std::cout << "In get_phy_param:" << std::endl;
    std::cout << "Unknown param name: " << pname << std::endl;
    return NAN;
  }
}

std::map<std::string, DTP_FLOAT> User_data::get_all_phy_params() {
  std::map<std::string, DTP_FLOAT> res;
  for (auto const& p: phyGetterDict) {
    res[p.first] = p.second(this->physical_params);
  }
  return res;
}

void User_data::assort_reactions()
{
  for (auto const& r: this->reactions) {
    if (r.itype == 61) {
      //if ((species.idx2name[r.idxReactants[0]] == "H2") ||
      //    (species.idx2name[r.idxReactants[0]] == "H")) {
      //  continue;
      //}
      this->auxdata.ads_reactions.push_back(r);
    } else if (r.itype == 62) {
      //if ((species.idx2name[r.idxReactants[0]] == "gH2") ||
      //    (species.idx2name[r.idxReactants[0]] == "gH")) {
      //  continue;
      //}
      this->auxdata.eva_reactions.push_back(r);
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


void User_data::assignElementsToSpecies(const Elements& elements) {
  for (auto const& s: this->species.name2idx) {
    this->species.elementsSpecies[s.second] =
    assignElementsToOneSpecies(s.first, elements);
  }
}


void User_data::assignElementsToSpecies() {
  assignElementsToSpecies(CONST::element_masses);
}


void User_data::calculateSpeciesMasses(const Elements& elements) {
  for (auto const& s: this->species.name2idx) {
    this->species.massSpecies[s.second] = 0.0;
    auto const& eleDict = this->species.elementsSpecies[s.second];
    for (auto const& e: eleDict) {
      this->species.massSpecies[s.second] +=
      elements.at(e.first) * ((double)e.second);
    }
  }
}

void User_data::calculateSpeciesMasses() {
  calculateSpeciesMasses(CONST::element_masses);
}

void User_data::calculateSpeciesVibFreqs() {
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

void User_data::calculateSpeciesDiffBarriers() {
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


void User_data::calculateSpeciesQuantumMobilities() {
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


void User_data::classifySpeciesByPhase() {
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


void User_data::calculateReactionHeat() {
  for (auto& r: this->reactions) {
    calculateOneReactionHeat(r, this->species.enthalpies);
  }
}


}
