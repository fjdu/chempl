#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <cstring>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include <iomanip>
#include <cerrno>
#include <stdexcept>
#include <utility>
#include "constants.hpp"
#include "math.h"

namespace TYPES {

typedef double DTP_FLOAT;
typedef double* DTP_Y;

class TimeDependency {
  public:
    std::string name;
    std::vector<double> ts, vs;
    TimeDependency(std::string name_,
      std::vector<double> ts_,
      std::vector<double> vs_);
};

typedef std::vector<TimeDependency> TimeDependencies;


class PhyParams {
public:
  void prep_params();

  int from_file(std::string fname);

  DTP_FLOAT n_gas, T_gas, T_dust, Av, G0_UV, chi_Xray, chi_cosmicray,
    dust2gas_num, dust2gas_mass, dust_material_density,
    dust_site_density, dust_radius, dust_crosssec, dust_albedo,
    mean_mol_weight, chemdesorption_factor, Ncol_H2, Ncol_H, Ncol_CO,
    dv_km_s, v_km_s;
  DTP_FLOAT t_max_year;
  TimeDependencies timeDependencies;

  void add_a_timedependency(std::string name, std::vector<double> ts,
    std::vector<double> vs);
  void remove_a_timedependency(std::string name);
  std::vector<std::string> get_timeDependency_names();
};


extern std::map<std::string, void (*)(PhyParams&, DTP_FLOAT)> phySetterDict;
extern std::map<std::string, DTP_FLOAT (*)(PhyParams&)> phyGetterDict;


class Reaction {
  public:
    int itype;
    std::vector<std::string> sReactants, sProducts;
    std::vector<int> idxReactants, idxProducts;
    std::vector<DTP_FLOAT> abc;
    std::vector<DTP_FLOAT> Trange;
    double drdy[2], rate, heat;
    Reaction();
    Reaction(
      std::vector<std::string> sR,
      std::vector<std::string> sP,
      std::vector<DTP_FLOAT> a_,
      std::vector<DTP_FLOAT> Tr,
      int it);
  friend bool operator== (const Reaction &r1, const Reaction &r2);
};


typedef std::vector<Reaction> Reactions;
typedef std::map<int, int> ReactionTypes;
typedef std::map<std::string, DTP_FLOAT> Elements;


class Species {
  public:
    std::map<std::string, int> name2idx;
    std::vector<std::string> idx2name;
    std::map<int, std::map<std::string, int> > elementsSpecies;
    std::map<int, DTP_FLOAT> massSpecies, enthalpies,
                             vibFreqs, diffBarriers, quantMobilities;
    std::set<int> gasSpecies, surfaceSpecies, mantleSpecies;
    std::vector<DTP_FLOAT> abundances;
    void allocate_abundances();
    Species();
    ~Species();
};


class AuxData {
  public:
    DTP_FLOAT t_calc, k_eva_tot, k_ads_tot, mant_tot, surf_tot;
    Reactions ads_reactions, eva_reactions;
    std::vector<int> ads_species, eva_species;
    int n_surf2mant, n_mant2surf;
    AuxData();
};


typedef DTP_FLOAT (*RateCalculator)(
    const DTP_FLOAT&, double *, Reaction&,
    const PhyParams&, const Species&, TYPES::AuxData&);
typedef std::vector<DTP_FLOAT> (*dRdyCalculator)(
    const DTP_FLOAT&, double *, Reaction&,
    const PhyParams&, const Species&, TYPES::AuxData&);


typedef std::map<int, RateCalculator> RateCalculators;
typedef std::map<int, dRdyCalculator> dRdyCalculators;

typedef std::map<std::string, std::string> PathsDict;

class Chem_data {
  public:
    PhyParams physical_params;
    Reactions reactions;
    ReactionTypes reaction_types;
    Species species;
    RateCalculators rate_calculators;
    AuxData auxdata;
    std::vector<int> dupli;
    Chem_data* ptr;
    DTP_FLOAT* y;

    bool updateIfIsDuplicate(const Reaction&);
    void add_reaction(Reaction);
    void modify_reaction(const int&,
      const std::map<std::string, std::vector<double> > &);
    void find_duplicate_reactions();
    void clear_reactions();
    void allocate_y();
    void allocate_y(int n);
    void deallocate_y();

    void set_phy_param(std::string name, DTP_FLOAT v);
    DTP_FLOAT get_phy_param(std::string name);
    std::map<std::string, DTP_FLOAT> get_all_phy_params();
    void assort_reactions();
    void assignElementsToSpecies(const Elements& elements);
    void assignElementsToSpecies();
    void calculateSpeciesMasses(const Elements& elements);
    void calculateSpeciesMasses();
    void calculateSpeciesVibFreqs();
    void calculateSpeciesDiffBarriers();
    void calculateSpeciesQuantumMobilities();
    void classifySpeciesByPhase();
    void calculateReactionHeat();
    double calculate_a_rate(double t, double *y, Reaction& r, bool updatePhyParams=false);
    std::vector<int> getFormationReactions(int iSpecies);
    std::vector<int> getDestructionReactions(int iSpecies);
    std::vector<std::pair<int, double> > getFormationReactionsWithRates(int iSpecies, double t, double* y);
    std::vector<std::pair<int, double> > getDestructionReactionsWithRates(int iSpecies, double t, double* y);
    Chem_data();
    ~Chem_data();
};


void update_phy_params(
    TYPES::DTP_FLOAT t,
    TYPES::PhyParams& p);


class Recorder {
  public:
    std::string fname;
    std::ofstream ofs;

    Recorder(std::string fname_);
    ~Recorder();

    int write_header(std::vector<std::string> names,
                     int fwidth=14);

    int write_row(double t, int neq, double *y,
                  const TYPES::PhyParams& p,
                  int fwidth=14, int prec=5);
};


std::map<std::string, int> assignElementsToOneSpecies(
    const std::string& name, const Elements& elements);

PathsDict loadPathConfig(std::string);


Reaction str2reaction(const std::string& str,
    int nReactants=3, int nProducts=4, int nABC=3,
    int lenSpeciesName=12, int lenABC=9, int nT=2, int lenT=6,
    int lenType=3, int rowlen_min=126);


void load_reactions(const std::string& fname,
           Chem_data& user_data,
           int nReactants=3, int nProducts=4, int nABC=3,
           int lenSpeciesName=12, int lenABC=9, int nT=2,
           int lenT=6, int lenType=3, int rowlen_min=126);

void loadInitialAbundances(Species& species, std::string fname);

void loadSpeciesEnthalpies(TYPES::Species& species, std::string fname);


}
#endif //TYPES_H
