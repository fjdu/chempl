#include <iostream>
#include <string>
#include <algorithm>
#include "logistics.hpp"
#include "types.hpp"
#include "rate_equation_lsode.hpp"
#include "constants.hpp"
#include "calculate_reaction_rate.hpp"
#include <ctime>

int main(int argc, char **argv)
{
  for (int i=0; i<argc; ++i) {
    std::cout << std::string(argv[i]) << " ";
  }
  std::cout << std::endl;

  std::string fname_path_config;
  if (argc < 2) {
    fname_path_config = "paths.dat";
  } else {
    fname_path_config = argv[1];
  }

  TYPES::PathsDict pdict = LOGIS::loadPathConfig(fname_path_config);
  for (auto const& p: pdict) {
    std::cout << p.first << " = " << p.second << std::endl;
  }

  TYPES::Chem_data user_data;

  user_data.physical_params.from_file(pdict["f_phy_params"]);
  user_data.physical_params.prep_params();
  for (auto const& p: TYPES::phySetterDict) {
    std::cout << p.first << " = "
              << TYPES::phyGetterDict[p.first](user_data.physical_params)
              << " " << std::endl;
  }

  LOGIS::load_reactions(pdict["f_reactions"], user_data);
  std::cout << "Number of reactions: "
            << user_data.reactions.size() << std::endl;
  std::cout << "Number of species: "
            << user_data.species.name2idx.size() << std::endl;
  std::cout << "Number of reaction types: "
            << user_data.reaction_types.size() << std::endl;
  for (auto const& i: user_data.reaction_types) {
    std::cout << i.first << ": " << i.second << "\n";
  }
  user_data.assort_reactions();
  std::cout << "Number of adsorption reactions: "
            << user_data.auxdata.ads_reactions.size() << std::endl;
  std::cout << "Number of evaporation reactions: "
            << user_data.auxdata.eva_reactions.size() << std::endl;

  user_data.assignElementsToSpecies();
  user_data.calculateSpeciesMasses();
  user_data.calculateSpeciesVibFreqs();
  user_data.calculateSpeciesDiffBarriers();
  user_data.calculateSpeciesQuantumMobilities();
  user_data.classifySpeciesByPhase();

  std::cout << "Number of gas species: "
            << user_data.species.gasSpecies.size() << std::endl;
  std::cout << "Number of surface species: "
            << user_data.species.surfaceSpecies.size() << std::endl;
  std::cout << "Number of mantle species: "
            << user_data.species.mantleSpecies.size() << std::endl;

  LOGIS::loadInitialAbundances(user_data.species, pdict["f_initial_abundances"]);
  std::cout << "Number of species with initial abundances: "
            << std::count_if(user_data.species.abundances.begin(),
                             user_data.species.abundances.end(),
                             [](TYPES::DTP_FLOAT v){return v>1e-40;})
            << std::endl;
  for (int i=0; i<user_data.species.idx2name.size(); ++i) {
    if (user_data.species.abundances[i] > 1.0e-40) {
      std::cout << user_data.species.idx2name[i] << " "
                << user_data.species.abundances[i] << std::endl;
    }
  }

  LOGIS::loadSpeciesEnthalpies(user_data.species, pdict["f_enthalpies"]);
  std::cout << "Number of species with enthalpies: "
            << user_data.species.enthalpies.size() << std::endl;
  //for (auto const& s: user_data.species->enthalpies) {
  //  std::cout << user_data.species->idx2name[s.first]
  //            << " " << s.second << std::endl;
  //}

  user_data.calculateReactionHeat();

  CALC_RATE::assignReactionHandlers(user_data);
  for (auto const& i: user_data.reaction_types) {
    if (user_data.rate_calculators.find(i.first) ==
        user_data.rate_calculators.end()) {
      std::cout << "Reaction type " << i.first
                << " has no handler!" << std::endl;
    }
  }

  TYPES::Recorder recorder(pdict["f_record"]);
  recorder.write_header(user_data.species.idx2name);

  RATE_EQ::Updater_RE updater_re;
  updater_re.set_user_data(&user_data);
  updater_re.set_sparse();
  updater_re.initialize_solver(1e-6, 1e-30);
  updater_re.set_solver_msg(1);
  //updater_re.set_solver_msg_lun(79);

  TYPES::DTP_FLOAT t=0.0, dt=1e-1, t_ratio=1.08;
  int NMAX = 5000;
  double t_max_seconds = user_data.physical_params.t_max_year * CONST::phy_SecondsPerYear;
  double *y = new double[updater_re.NEQ];
  for (int i=0; i<updater_re.NEQ; ++i) {
    y[i] = user_data.species.abundances[i];
  }

  recorder.write_row(t, updater_re.NEQ, y, user_data.physical_params);
  std::cout << std::endl;

  clock_t rt_begin = std::clock();
  for (int i=0; i<NMAX; ++i) {
    t = updater_re.update(t, dt, y);
    recorder.write_row(t, updater_re.NEQ,
                       y, user_data.physical_params);
    if (t >= t_max_seconds) {
      break;
    }
    dt *= t_ratio;
    if (t + dt > t_max_seconds) {
      dt = t_max_seconds - t;
    }
    if (updater_re.ISTATE != 2) {
      std::cout << "Failed. ISTATE=" << updater_re.ISTATE << "\n" << std::endl;
      if ((updater_re.ISTATE == -1) ||
          (updater_re.ISTATE == -4) ||
          (updater_re.ISTATE == -5)) {
        updater_re.ISTATE = 3;
      } else {
        break;
      }
    }
  }
  clock_t rt_end = std::clock();
  double elapsed_secs = double(rt_end - rt_begin) / CLOCKS_PER_SEC;
  std::cout << elapsed_secs << " seconds elapsed." << std::endl;

  delete [] y;

  return 0;
}
