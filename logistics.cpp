#ifndef LOGIS_H
#define LOGIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <regex>
#include <set>
#include "types.hpp"
#include "utils.hpp"
#include "math.h"

namespace LOGIS {


void load_reactions(const std::string& fname, TYPES::User_data& user_data,
    int nReactants=3, int nProducts=4, int nABC=3, int lenSpeciesName=12,
    int lenABC=9, int nT=2, int lenT=6, int lenType=3, int rowlen_min=126)
{
  std::ifstream inputFile(fname);
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");

  TYPES::Species& species = user_data.species;
  TYPES::Reactions& reactions = user_data.reactions;
  TYPES::ReactionTypes& r_types = user_data.reaction_types;

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {
      if (std::regex_match(line, comment)) {
        // std::cout << "Comments!" << std::endl;
        continue;
      }
      if (std::regex_match(line, emptyline)) {
        // std::cout << "Empty line!" << std::endl;
        continue;
      }
      if (line.size() < rowlen_min) {
        // std::cout << "Too short line!" << std::endl;
        continue;
      }

      TYPES::Reaction reaction;

      for (int i=0; i<nReactants+nProducts; ++i) {
        std::string name = line.substr(i*lenSpeciesName, lenSpeciesName);
        name = UTILS::trim(name);

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
            reaction.sReactants.push_back(name);
          } else {
            reaction.idxProducts.push_back(iSpecies);
            reaction.sProducts.push_back(name);
          }
        }
      }

      int iABC = (nReactants + nProducts) * lenSpeciesName;
      for (int i=0; i<nABC; ++i) {
        std::string tmp = line.substr(iABC+i*lenABC, lenABC);
        tmp = UTILS::trim(tmp);
        if (tmp.size() > 0) {
          reaction.abc.push_back(std::stod(tmp));
        } else {
          reaction.abc.push_back(0.0);
        }
      }

      int iT = iABC + nABC * lenABC;
      for (int i=0; i<nT; ++i) {
        std::string tmp = line.substr(iT+i*lenT, lenT);
        tmp = UTILS::trim(tmp);
        if (tmp.size() > 0) {
          reaction.Trange.push_back(std::stod(tmp));
        } else {
          reaction.Trange.push_back(0.0);
        }
      }

      int iType = iT + nT * lenT;
      std::string tmp = line.substr(iType, lenType);
      try {
        reaction.itype = std::stoi(tmp);
      } catch (...) {
        std::cerr << "Exception in load_reactions: "
                  << tmp << std::endl;
      }
      if (r_types.find(reaction.itype) == r_types.end()) {
        r_types[reaction.itype] = 1;
      } else {
        r_types[reaction.itype] += 1;
      }

      reactions.push_back(reaction);
    }
  } else {
    std::cerr << "Error in load_reactions: "
              << fname << std::endl;
  }

  inputFile.close();
}


/*
void assort_reactions(const TYPES::Reactions& reactions,
                      const TYPES::Species& species,
                      TYPES::AuxData& m)
{
  for (auto const& r: reactions) {
    if (r.itype == 61) {
      //if ((species.idx2name[r.idxReactants[0]] == "H2") ||
      //    (species.idx2name[r.idxReactants[0]] == "H")) {
      //  continue;
      //}
      m.ads_reactions.push_back(r);
    } else if (r.itype == 62) {
      //if ((species.idx2name[r.idxReactants[0]] == "gH2") ||
      //    (species.idx2name[r.idxReactants[0]] == "gH")) {
      //  continue;
      //}
      m.eva_reactions.push_back(r);
    }
  }
}


std::map<std::string, int> assignElementsToOneSpecies(
    const std::string& name, const TYPES::Elements& elements)
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


int assignElementsToSpecies(TYPES::Species& species, const TYPES::Elements& elements) {
  for (auto const& s: species.name2idx) {
    species.elementsSpecies[s.second] = assignElementsToOneSpecies(s.first, elements);
  }
  return 0;
}


int classifySpeciesByPhase(TYPES::Species& species) {
  for (auto const& s: species.name2idx) {
    if (s.first[0] == 'm') {
      species.mantleSpecies.insert(s.second);
    } else if (s.first[0] == 'g') {
      species.surfaceSpecies.insert(s.second);
    } else {
      species.gasSpecies.insert(s.second);
    }
  }
  return 0;
}


int calculateSpeciesMasses(TYPES::Species& species, const TYPES::Elements& elements) {
  for (auto const& s: species.name2idx) {
    species.massSpecies[s.second] = 0.0;
    auto const& eleDict = species.elementsSpecies[s.second];
    for (auto const& e: eleDict) {
      species.massSpecies[s.second] += elements.at(e.first) * ((double)e.second);
    }
  }
  return 0;
}


int calculateSpeciesVibFreqs(TYPES::Species& species, TYPES::Reactions& reactions) {
  for (auto const& r: reactions) {
    if (r.itype == 62) {
      species.vibFreqs[r.idxReactants[0]] =
        sqrt(2.0 * CONST::phy_SitesDensity_CGS
           * CONST::phy_kBoltzmann_CGS * r.abc[2]
           / (CONST::PI * CONST::PI)
           / (CONST::phy_mProton_CGS * species.massSpecies[r.idxReactants[0]]));
    }
  }
  for (auto const& r: reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (species.vibFreqs.find(i) == species.vibFreqs.end()) {
          species.vibFreqs[i] = CONST::phy_vibFreqDefault;
        }
      }
    }
  }
  return 0;
}


int calculateSpeciesDiffBarriers(TYPES::Species& species, TYPES::Reactions& reactions) {
  for (auto const& r: reactions) {
    if (r.itype == 62) {
      species.diffBarriers[r.idxReactants[0]] = r.abc[2] * CONST::phy_Diff2DesorRatio;
    }
  }
  for (auto const& r: reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (species.diffBarriers.find(i) == species.diffBarriers.end()) {
          species.diffBarriers[i] = CONST::phy_DiffBarrierDefault;
        }
      }
    }
  }
  return 0;
}


int calculateSpeciesQuantumMobilities(TYPES::Species& species, TYPES::Reactions& reactions) {
  for (auto const& r: reactions) {
    if (r.itype == 62) {
      species.quantMobilities[r.idxReactants[0]] =
          species.vibFreqs[r.idxReactants[0]] *
          exp(-2.0 * CONST::phy_DiffBarrierWidth_CGS
             / CONST::phy_hbarPlanck_CGS
             * sqrt(2.0 * species.massSpecies[r.idxReactants[0]] * CONST::phy_mProton_CGS
                  * CONST::phy_kBoltzmann_CGS
                  * species.diffBarriers[r.idxReactants[0]]));
    }
  }
  for (auto const& r: reactions) {
    if ((r.itype == 63) || (r.itype == 64)) {
      for (auto const& i: r.idxReactants) {
        if (species.quantMobilities.find(i) == species.quantMobilities.end()) {
          species.quantMobilities[i] = 0.0;
        }
      }
    }
  }
  return 0;
}
*/


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


TYPES::PathsDict loadPathConfig(std::string fname) {
  std::string line;
  std::regex comment(R"(^[!#].*$)");
  std::regex emptyline(R"(^\s*$)");
  std::regex entry(R"(^\s*(\S+)\s*=\s*([^!#]*[^ !#]+)\s*(?:([!#])|($)))");
  std::ifstream inputFile(fname);
  TYPES::PathsDict pd;
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


}

#endif // LOGIS_H
