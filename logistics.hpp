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
#include "math.h"

namespace LOGIS {


extern void load_reactions(const std::string& fname,
           TYPES::User_data& user_data,
           int nReactants=3, int nProducts=4, int nABC=3,
           int lenSpeciesName=12, int lenABC=9, int nT=2, int lenT=6, int lenType=3,
           int rowlen_min=126);

/*
extern void assort_reactions(const TYPES::Reactions& reactions, const TYPES::Species& species, TYPES::AuxData& m);


extern std::map<std::string, int> assignElementsToOneSpecies(
           const std::string& name, const TYPES::Elements& elements);


extern int assignElementsToSpecies(TYPES::Species& species, const TYPES::Elements& elements);


extern int classifySpeciesByPhase(TYPES::Species& species);


extern int calculateSpeciesMasses(TYPES::Species& species, const TYPES::Elements& elements);


extern int calculateSpeciesVibFreqs(TYPES::Species& species, TYPES::Reactions& reactions);


extern int calculateSpeciesDiffBarriers(TYPES::Species& species, TYPES::Reactions& reactions);


extern int calculateSpeciesQuantumMobilities(TYPES::Species& species, TYPES::Reactions& reactions);
*/


extern void loadInitialAbundances(TYPES::Species& species, std::string fname);

extern void loadSpeciesEnthalpies(TYPES::Species& species, std::string fname);

extern TYPES::PathsDict loadPathConfig(std::string fname);


}

#endif // LOGIS_H
