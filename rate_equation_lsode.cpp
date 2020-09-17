#include <iostream>
#include <vector>
#include <algorithm>
#include "types.hpp"
#include "rate_equation_lsode.hpp"
#include "calculate_reaction_rate.hpp"

namespace RATE_EQ {

int Updater_RE::makeSparse(
    const TYPES::Reactions& reactions,
    std::vector<std::vector<bool> >& sps) {
  for (std::size_t i=0; i<sps.size(); ++i) {
    for (std::size_t j=0; j<sps[i].size(); ++j) {
      sps[i][j] = false;
    }
  }
  for (auto const& r: reactions) {
    for (auto const& i: r.idxReactants) {
      for (auto const& j: r.idxReactants) {
        sps[j][i] = true;
      }
      for (auto const& j: r.idxProducts) {
        sps[j][i] = true;
      }
    }
  }
  int nnz = 0;
  for (int i=0; i<NEQ; ++i) {
    for (int j=0; j<NEQ; ++j) {
      if (sps[i][j]) {
        ++nnz;
      }
    }
  }
  return nnz;
}


Updater_RE::Updater_RE() {
  LRW = 0;
  LIW = 0;
  data = nullptr;
  RWORK = nullptr;
  IWORK = nullptr;
  RSAV = nullptr;
  ISAV = nullptr;
}


Updater_RE::~Updater_RE() {
  data = nullptr;
  if (IWORK != nullptr) {delete [] IWORK; IWORK = nullptr;}
  if (RWORK != nullptr) {delete [] RWORK; RWORK = nullptr;}
  if (RSAV != nullptr) {delete [] RSAV; RSAV = nullptr;}
  if (ISAV != nullptr) {delete [] ISAV; ISAV = nullptr;}
}


TYPES::DTP_FLOAT Updater_RE::update(double t, double dt, double *y)
{
  double t0 = t;
  double tout = t + dt;

  dlsodes_w(f, NEQ, y, &t, tout, ITOL, RTOL, ATOL, ITASK,
            &ISTATE, IOPT, RWORK, LRW, IWORK, LIW, jac, MF);

  std::cout << (char)27 << "[A"
            << std::setw(9) << std::setprecision(3) << std::left
            << t0 << " -> "
            << std::setw(9) << std::setprecision(3) << std::left
            << t
            << " ISTATE=" << ISTATE
            << " #f="   << IWORK[11]
            << " #jac=" << IWORK[12]
            << " LRW=" << IWORK[16]
            << " LIW=" << IWORK[17]
            << " NNZ=" << IWORK[18]
            << " SID=" << SOLVER_ID
            << std::endl;
  return t;
}

void Updater_RE::set_user_data(TYPES::Chem_data *data_) {
  data = data_;
  NEQ = data->species.idx2name.size();
  // std::cout << "Solver data pointer -> " << data << "\n" << std::endl;
}

void Updater_RE::allocate_sparse() {
  if (sparseMaskJac.size() != (std::size_t)NEQ) {
    for (auto& s: sparseMaskJac) {std::vector<bool>().swap(s);}
    std::vector<std::vector<bool> >().swap(sparseMaskJac);
    for (int i=0; i<NEQ; ++i) {
      sparseMaskJac.push_back(std::vector<bool>(NEQ, false));
    }
  }
}

int Updater_RE::initialize_solver(
    double reltol,
    double abstol,
    int mf, //021: use Jac; 022: not.
    int LRW_F,
    int solver_id)
{
  MF = mf;
  IOPT = 1;
  ITOL = 1;
  RTOL = reltol;
  ATOL = abstol;
  ITASK = 1;
  ISTATE = 1;
  SOLVER_ID = solver_id;

  std::cout << "Making sparse..." << std::endl;
  NNZ = makeSparse(data->reactions, sparseMaskJac);
  std::cout << "NNZ = " << NNZ << " ("
            << (double)NNZ / (double)(NEQ*NEQ) << ")" << std::endl;

  int lrw = 20 + 20 * NEQ + LRW_F * NNZ;
  int liw = 31 + NEQ + NNZ;

  if ((IWORK == nullptr) || (RWORK == nullptr) ||
      (LRW != lrw) || (LIW != liw)) {
    LRW = lrw;
    LIW = liw;

    if (IWORK != nullptr) {delete [] IWORK; IWORK = nullptr;}
    if (RWORK != nullptr) {delete [] RWORK; RWORK = nullptr;}
    IWORK = new int[LIW];
    RWORK = new double[LRW];

    std::cout << "RWORK size = " << LRW << std::endl;
    std::cout << "IWORK size = " << LIW << std::endl;
  }

  std::fill(IWORK, IWORK+LIW, 0);
  std::fill(RWORK, RWORK+LRW, 0.0);
  IWORK[4] = 5;
  IWORK[5] = 5000;
  IWORK[6] = 10;

  int k = 1;
  IWORK[30] = k;
  for (int i=0; i<NEQ; ++i) {
    for (int j=0; j<NEQ; ++j) {
      if (sparseMaskJac[j][i]) {
        IWORK[30 + NEQ + k] = j+1;
        ++k;
      }
    }
    IWORK[31+i] = k;
  }
  if (NNZ != (k-1)) {
    std::cerr << "NNZ != (k-1)" << k-1 << std::endl;
  }

  return 0;
}

void Updater_RE::set_solver_msg(int mflag) {
  xsetf_w(&mflag);
}

void Updater_RE::set_solver_msg_lun(int lun) {
  xsetun_w(&lun);
}

void Updater_RE::allocate_rsav_isav() {
  lrsav = 256;
  lisav = 128;
  if (RSAV == nullptr) {
    RSAV = new double[lrsav];
  }
  if (ISAV == nullptr) {
    ISAV = new int[lisav];
  }
}

void Updater_RE::save_restore_common_block(int job) {
  dsrcms_w(RSAV, &lrsav, ISAV, &lisav, &job);
}

void Updater_RE::f(int *neq, double *t, double *y, double *ydot)
{
  for (int i=0; i<*neq; ++i) {
    ydot[i] = 0.0;
  }

  TYPES::update_phy_params(*t, data->physical_params);
  CALC_RATE::update_surfmant(
      *t, y, data->physical_params,
      data->species, data->auxdata);

  for (auto &reaction: data->reactions) {
    TYPES::DTP_FLOAT r =
      (((data->rate_calculators))[reaction.itype])(*t, y,
        reaction, data->physical_params,
        data->species, data->auxdata);
    for (auto const& i: reaction.idxReactants) {
      ydot[i] -= r;
    }
    for (auto const& i: reaction.idxProducts) {
      ydot[i] += r;
    }
  }
}


void Updater_RE::jac(int *neq, double *t, double *y, int *j, double *ian, double *jan, double *pdj)
{
  int jc = (*j) - 1;
  for (auto const& r: data->reactions) {
    for (std::size_t i=0; i < r.idxReactants.size(); ++i) {
      if (r.idxReactants[i] != jc) {
        continue;
      }
      for (auto const& k: r.idxReactants) {
        pdj[k] -= r.drdy[i];
      }
      for (auto const& k: r.idxProducts) {
        pdj[k] += r.drdy[i];
      }
    }
  }
}


// The have to be static because they are used by static functions f and jac,
// which themselves have to be static because they are used by the external
// dlsodes_w.
// The static variables have to be initialized here.
TYPES::Chem_data *Updater_RE::data;


}
