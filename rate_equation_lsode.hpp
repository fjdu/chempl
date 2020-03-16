#include <iostream>
#include <vector>
#include <algorithm>
#include "types.hpp"

#ifndef RATE_EQ_H
#define RATE_EQ_H

extern "C" {

void dlsodes_w(
  void (*)(int *, double *, double *, double *), //f
  int *NEQ, double *y, double *t, double *tout, int *ITOL,
  double *RTOL, double *ATOL, int *ITASK, int *ISTATE, int *IOPT,
  double *RWORK, int *LRW, int *IWORK, int *LIW,
  void (*)(int *, double *, double *, int *, double *, double *, double *),
  int *MF);

void xsetf_w(int* mflag);
void xsetun_w(int* lun);

}

namespace RATE_EQ {

class Updater_RE {
  public:
    static void f(int *neq, double *t, double *y, double *ydot);
    static void jac(int *neq, double *t, double *y, int *j, double *ian, double *jan, double *pdj);
    TYPES::DTP_FLOAT update(double t, double dt, double *y);
    int initialize_solver(double reltol=1e-6, double abstol=1e-30, int mf=21, int LRW_F=6);
    int makeSparse(const TYPES::Reactions& reactions, std::vector<std::vector<bool> >& sps);
    void set_user_data(TYPES::User_data *data_);
    void set_solver_msg(int mflag);
    void set_solver_msg_lun(int lun);

    static TYPES::User_data *data;
    static int *IWORK;
    static double *RWORK;
    int NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF, NNZ;
    double RTOL, ATOL;
    std::vector<std::vector<bool> > sparseMaskJac;

    void set_MF(int i_) {MF = i_;}
    void set_IOPT(int i_) {IOPT = i_;}
    void set_ITOL(int i_) {ITOL = i_;}
    void set_ITASK(int i_) {ITASK = i_;}
    void set_ISTATE(int i_) {ISTATE = i_;}
    void set_RTOL(double i_) {RTOL = i_;}
    void set_ATOL(double i_) {ATOL = i_;}

    Updater_RE() {RWORK = nullptr; IWORK = nullptr;}

    ~Updater_RE() {
      data = nullptr;
      if (IWORK != nullptr) {delete [] IWORK; IWORK = nullptr;}
      if (RWORK != nullptr) {delete [] RWORK; RWORK = nullptr;}
    }

};


}

#endif //RATE_EQ_H
