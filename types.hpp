#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <regex>
#include "math.h"
#include "constants.hpp"
#include <iomanip>
#include <cerrno>
#include <stdexcept>

namespace TYPES {

typedef double DTP_FLOAT;
typedef double* DTP_Y;


class PhyParams {
public:
  inline DTP_FLOAT get_n_gas() const {return n_gas;}
  inline DTP_FLOAT get_T_gas() const {return T_gas;}
  inline DTP_FLOAT get_T_dust() const {return T_dust;}
  inline DTP_FLOAT get_Av() const {return Av;}
  inline DTP_FLOAT get_G0_UV() const {return G0_UV;}
  inline DTP_FLOAT get_chi_Xray() const {return chi_Xray;}
  inline DTP_FLOAT get_chi_cosmicray() const {return chi_cosmicray;}
  inline DTP_FLOAT get_dust2gas_num() const {return dust2gas_num;}
  inline DTP_FLOAT get_dust2gas_mass() const {return dust2gas_mass;}
  inline DTP_FLOAT get_dust_material_density() const {return dust_material_density;}
  inline DTP_FLOAT get_dust_site_density() const {return dust_site_density;}
  inline DTP_FLOAT get_dust_radius() const {return dust_radius;}
  inline DTP_FLOAT get_dust_crosssec() const {return dust_crosssec;}
  inline DTP_FLOAT get_dust_albedo() const {return dust_albedo;}
  inline DTP_FLOAT get_mean_mol_weight() const {return mean_mol_weight;}
  inline DTP_FLOAT get_chemdesorption_factor() const {return chemdesorption_factor;}
  inline DTP_FLOAT get_t_max_year() const {return t_max_year;}

  DTP_FLOAT get_n_gas(DTP_FLOAT t) const;
  DTP_FLOAT get_T_gas(DTP_FLOAT t) const;
  DTP_FLOAT get_T_dust(DTP_FLOAT t) const;
  inline DTP_FLOAT get_Av(DTP_FLOAT t) const;
  inline DTP_FLOAT get_G0_UV(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chi_Xray(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chi_cosmicray(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust2gas_num(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust2gas_mass(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_material_density(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_site_density(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_radius(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_crosssec(DTP_FLOAT t) const;
  inline DTP_FLOAT get_dust_albedo(DTP_FLOAT t) const;
  inline DTP_FLOAT get_mean_mol_weight(DTP_FLOAT t) const;
  inline DTP_FLOAT get_chemdesorption_factor(DTP_FLOAT t) const;

  inline void set_n_gas(DTP_FLOAT _v);
  inline void set_T_gas(DTP_FLOAT _v);
  inline void set_T_dust(DTP_FLOAT _v);
  inline void set_Av(DTP_FLOAT _v);
  inline void set_G0_UV(DTP_FLOAT _v);
  inline void set_chi_Xray(DTP_FLOAT _v);
  inline void set_chi_cosmicray(DTP_FLOAT _v);
  inline void set_dust2gas_num(DTP_FLOAT _v);
  inline void set_dust2gas_mass(DTP_FLOAT _v);
  inline void set_dust_material_density(DTP_FLOAT _v);
  inline void set_dust_site_density(DTP_FLOAT _v);
  inline void set_dust_radius(DTP_FLOAT _v);
  inline void set_dust_crosssec(DTP_FLOAT _v);
  inline void set_dust_albedo(DTP_FLOAT _v);
  inline void set_mean_mol_weight(DTP_FLOAT _v);
  inline void set_chemdesorption_factor(DTP_FLOAT _v);
  inline void set_t_max_year(DTP_FLOAT _v);

  void prep_params();

  int from_file(std::string fname);

private:
  DTP_FLOAT n_gas, T_gas, T_dust, Av, G0_UV, chi_Xray, chi_cosmicray,
    dust2gas_num, dust2gas_mass, dust_material_density,
    dust_site_density, dust_radius, dust_crosssec, dust_albedo,
    mean_mol_weight, chemdesorption_factor;
  double t_max_year;
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
    double drdy[2], rate;
    Reaction() {
      drdy[0] = NAN;
      drdy[1] = NAN;
      rate = 0.0;
    }
    Reaction(
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
    Species() {
    }
    ~Species() {
    }
};


class AuxData {
  public:
    DTP_FLOAT t_calc, k_eva_tot, k_ads_tot, mant_tot, surf_tot;
    Reactions ads_reactions, eva_reactions;
    AuxData(): t_calc(NAN), k_eva_tot(0.0), k_ads_tot(0.0),
                 mant_tot(0.0), surf_tot(0.0) {}
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

class User_data {
  public:
    PhyParams physical_params;
    Reactions reactions;
    ReactionTypes reaction_types;
    Species species;
    RateCalculators rate_calculators;
    AuxData auxdata;
    User_data* ptr;
    DTP_FLOAT* y;

    void add_reaction(const Reaction&);
    void remove_reaction(const Reaction&);
    void clear_reactions();
    void allocate_y() {
      if (y != nullptr) {
        delete [] y;
      }
      y = new DTP_FLOAT[species.idx2name.size()];
    }
    void deallocate_y() {
      if (y != nullptr) {
        delete [] y;
        y = nullptr;
      }
    }

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
    User_data() {
      ptr = this;
      y = nullptr;
    }
    ~User_data() {
      deallocate_y();
    }
};


class Recorder {
  public:
    std::string fname;
    std::ofstream ofs;

    Recorder(std::string fname_): fname(fname_) {
      ofs.open(fname);
      if (ofs.fail()) {
        std::cerr << "Fail to open file: "
                  << fname << std::endl;
        throw std::runtime_error(std::strerror(errno));
      }
    }

    ~Recorder() {
      ofs.close();
    }

    int write_header(std::vector<std::string> names,
                     int fwidth=14) {
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

    int write_row(double t, int neq, double *y,
                  const TYPES::PhyParams& p,
                  int fwidth=14, int prec=5) {
      ofs << std::setw(fwidth) << std::left << std::scientific
          << std::setprecision(prec) << t / CONST::phy_SecondsPerYear;
      ofs << std::setw(fwidth) << std::left << std::scientific
        << std::setprecision(prec) << p.get_T_gas(t);
      ofs << std::setw(fwidth) << std::left << std::scientific
        << std::setprecision(prec) << p.get_T_dust(t);
      ofs << std::setw(fwidth) << std::left << std::scientific
        << std::setprecision(prec) << p.get_n_gas(t);
      for (int i=0; i<neq; ++i) {
        ofs << std::setw(fwidth) << std::left << std::scientific
          << std::setprecision(prec) << y[i];
      }
      ofs << std::endl;
      return 0;
    }
};


}
#endif //TYPES_H
