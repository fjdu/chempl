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


TYPES::DTP_FLOAT Updater_RE::update(double t, double dt, double *y)
{
  double t0 = t;
  double tout = t + dt;

  dlsodes_w(f, &NEQ, y, &t, &tout, &ITOL, &RTOL, &ATOL, &ITASK,
            &ISTATE, &IOPT, RWORK, &LRW, IWORK, &LIW, jac, &MF);

  std::cout << (char)27 << "[A"
            << std::setw(10) << std::setprecision(4) << std::left
            << t0 << " -> "
            << std::setw(10) << std::setprecision(4) << std::left
            << t
            << " ISTATE=" << ISTATE
            << " #f="   << IWORK[11]
            << " #jac=" << IWORK[12]
            << " LRW=" << IWORK[16]
            << " LIW=" << IWORK[17]
            << " NNZ=" << IWORK[18]
            << std::endl;
  return t;
}

void Updater_RE::set_user_data(TYPES::User_data *data_) {
  data = data_;
  NEQ = data->species.idx2name.size();
  if (sparseMaskJac.size() != NEQ) {
    for (auto& s: sparseMaskJac) {std::vector<bool>().swap(s);}
    std::vector<std::vector<bool> >().swap(sparseMaskJac);
  }
  for (int i=0; i<NEQ; ++i) {
    sparseMaskJac.push_back(std::vector<bool>(NEQ, false));
  }
}

int Updater_RE::initialize_solver(
    double reltol,
    double abstol,
    int mf, //021: use Jac; 022: not.
    int LRW_F)
{
  MF = mf;
  IOPT = 1;
  ITOL = 1;
  RTOL = reltol;
  ATOL = abstol;
  ITASK = 1;
  ISTATE = 1;

  NNZ = makeSparse(data->reactions, sparseMaskJac);
  std::cout << "NNZ = " << NNZ << " ("
            << (double)NNZ / (double)(NEQ*NEQ) << ")" << std::endl;

  LRW = 20 + 20 * NEQ + LRW_F * NNZ;
  LIW = 31 + NEQ + NNZ;

  if (IWORK != nullptr) {delete [] IWORK; IWORK = nullptr;}
  if (RWORK != nullptr) {delete [] RWORK; RWORK = nullptr;}
  IWORK = new int[LIW];
  RWORK = new double[LRW];

  std::cout << "RWORK size = " << LRW << std::endl;
  std::cout << "IWORK size = " << LIW << std::endl;

  // for (int i=0; i<LRW; ++i) {RWORK[i] = 0.0;}
  // for (int i=0; i<LIW; ++i) {IWORK[i] = 0;}
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


void Updater_RE::f(int *neq, double *t, double *y, double *ydot)
{
  for (int i=0; i<*neq; ++i) {
    ydot[i] = 0.0;
  }

  CALC_RATE::update_surfmant(
      *t, y, data->physical_params,
      data->species, data->auxdata);

  for (auto &reaction: data->reactions) {
    if (reaction.itype == 67) {
      continue;
    }
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
    if (r.itype == 67) {
      continue;
    }
    for (int i=0; i < r.idxReactants.size(); ++i) {
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
TYPES::User_data *Updater_RE::data;
int *Updater_RE::IWORK;
double *Updater_RE::RWORK;


}
