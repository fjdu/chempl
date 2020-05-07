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


TYPES::Reaction str2reaction(const std::string& str,
    int nReactants=3, int nProducts=4, int nABC=3, int lenSpeciesName=12,
    int lenABC=9, int nT=2, int lenT=6, int lenType=3, int rowlen_min=126) {
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
    std::cerr << "Exception in load_reactions: "
              << tmp << std::endl;
  }

  return reaction;
}


void load_reactions(const std::string& fname, TYPES::Chem_data& user_data,
    int nReactants=3, int nProducts=4, int nABC=3, int lenSpeciesName=12,
    int lenABC=9, int nT=2, int lenT=6, int lenType=3, int rowlen_min=126)
{
  std::ifstream inputFile(fname);
  std::string line;

  TYPES::Species& species = user_data.species;
  TYPES::Reactions& reactions = user_data.reactions;
  TYPES::ReactionTypes& r_types = user_data.reaction_types;

  if (inputFile.good()) {
    while (std::getline(inputFile, line)) {

      TYPES::Reaction rthis = str2reaction(line);
      if (rthis.sReactants.size() == 0) {
        continue;
      }

      bool isDuplicate = false;
      for (auto &r: reactions) {
        if ((r.sReactants == rthis.sReactants) &&
            (r.sProducts == rthis.sProducts) &&
            (r.itype == rthis.itype)) {
          int iTr;
          for (iTr=r.Trange.size()-2; iTr>=0; iTr -= 2) {
            if (r.Trange[iTr] <= rthis.Trange[0]) {
              break;
            }
          }
          r.Trange.insert(r.Trange.begin()+iTr+2, rthis.Trange.begin(), rthis.Trange.end());
          r.abc.insert(r.abc.begin() + (iTr*3)/2+3, rthis.abc.begin(), rthis.abc.end());
          isDuplicate = true;
          break;
        }
      }
      if (isDuplicate) {
        continue;
      }

      for (auto s: rthis.sReactants) {
        if (species.name2idx.find(s) == species.name2idx.end()) {
          species.name2idx[s] = species.name2idx.size();
          species.idx2name.push_back(s);
        }
        rthis.idxReactants.push_back(species.name2idx[s]);
      }
      for (auto s: rthis.sProducts) {
        if (species.name2idx.find(s) == species.name2idx.end()) {
          species.name2idx[s] = species.name2idx.size();
          species.idx2name.push_back(s);
        }
        rthis.idxProducts.push_back(species.name2idx[s]);
      }

      if (r_types.find(rthis.itype) == r_types.end()) {
        r_types[rthis.itype] = 1;
      } else {
        r_types[rthis.itype] += 1;
      }

      reactions.push_back(rthis);
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
