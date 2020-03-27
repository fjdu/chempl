# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from myconsts import Consts

cdef extern from "types.hpp" namespace "TYPES":

  cdef cppclass Reaction:
    int itype
    vector[string] sReactants, sProducts
    vector[double] abc, Trange
    double drdy[2]
    double rate
    double heat
    Reaction() except +
    Reaction(vector[string] sReactants,
                vector[string] sProducts,
                vector[double] abc,
                vector[double] Trange,
                int itype)
    
  cdef cppclass Species:
    cppmap[string, int] name2idx
    vector[string] idx2name
    cppmap[int, cppmap[string, int]] elementsSpecies
    cppmap[int, double] massSpecies, enthalpies,
    cppmap[int, double] vibFreqs, diffBarriers, quantMobilities
    cppset[int] gasSpecies, surfaceSpecies, mantleSpecies
    vector[double] abundances

  cdef cppclass AuxData:
    pass

  cdef cppclass PhyParams:
    void prep_params()
    int from_file(string fname)

  ctypedef vector[Reaction] Reactions
  ctypedef cppmap[string, double] Elements
  ctypedef cppmap[int, int] ReactionTypes
  ctypedef double (*RateCalculator)(const double&, double *,
           Reaction&, const PhyParams&, const Species&, AuxData&)
  ctypedef cppmap[int, RateCalculator] RateCalculators

  cdef cppclass User_data:
    void add_reaction(Reaction& rs)
    void clear_reactions()
    void set_phy_param(string name, double v)
    double get_phy_param(string name)
    cppmap[string, double] get_all_phy_params()
    void assort_reactions()
    void assignElementsToSpecies()
    void assignElementsToSpecies(Elements& elements)
    void calculateSpeciesMasses()
    void calculateSpeciesMasses(Elements& elements)
    void calculateSpeciesVibFreqs()
    void calculateSpeciesDiffBarriers()
    void calculateSpeciesQuantumMobilities()
    void calculateReactionHeat()
    void classifySpeciesByPhase()
    void allocate_y()
    void deallocate_y()
    Reactions reactions
    PhyParams physical_params
    Species species
    ReactionTypes reaction_types
    User_data* ptr
    RateCalculators rate_calculators
    double* y



cdef extern from "logistics.hpp" namespace "LOGIS":
  void load_reactions(const string& fname, User_data& user_data,
    int nReactants, int nProducts, int nABC, int lenSpeciesName,
    int lenABC, int nT, int lenT, int lenType, int rowlen_min)
  int loadInitialAbundances(Species& species, string fname)
  int loadSpeciesEnthalpies(Species& species, string fname)


cdef extern from "calculate_reaction_rate.hpp" namespace "CALC_RATE":
  void assignReactionHandlers(User_data&)
  void assignAReactionHandler(RateCalculators& rcs,
                              const RateCalculator& rc,
                              const int& itype)
  double rate_photodissociation_H2(
    const double& t,
    const double* y,
    Reaction& r,
    const PhyParams& p,
    const Species& s,
    AuxData& m)


cdef extern from "rate_equation_lsode.hpp" namespace "RATE_EQ":
  cdef cppclass Updater_RE:
    void set_user_data(User_data* udata)
    void set_sparse()
    int initialize_solver(double reltol, double abstol, int mf, int LRW_F, int solver_id)
    void allocate_rsav_isav()
    void save_restore_common_block(int job)
    double update(double t, double dt, double* y)
    void set_solver_msg(int mflag)
    void set_solver_msg_lun(int lun)
    void set_MF(int i_)
    void set_IOPT(int i_)
    void set_ITOL(int i_)
    void set_ITASK(int i_)
    void set_ISTATE(int i_)
    void set_RTOL(double i_)
    void set_ATOL(double i_)

    int NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF, NNZ
    double RTOL, ATOL


cdef class pyUserData:

  cdef User_data user_data
  cdef Updater_RE updater_re

  def set_solver(self, rtol=1e-6, atol=1e-30, mf=21, LRW_F=6,
                 showmsg=1, msglun=6, solver_id=0):
    self.updater_re.set_user_data(self.user_data.ptr)
    self.updater_re.set_sparse()
    self.updater_re.initialize_solver(rtol, atol, mf, LRW_F, solver_id)
    self.updater_re.set_solver_msg(showmsg);
    self.updater_re.set_solver_msg_lun(msglun);
    self.updater_re.allocate_rsav_isav()
    self.user_data.allocate_y()

  def deallocate_y(self):
    self.user_data.deallocate_y()

  def get_solver_internals(self):
    return {
      'NEQ':    self.updater_re.NEQ,
      'ITOL':   self.updater_re.ITOL,
      'ITASK':  self.updater_re.ITASK,
      'ISTATE': self.updater_re.ISTATE,
      'IOPT':   self.updater_re.IOPT,
      'LRW':    self.updater_re.LRW,
      'LIW':    self.updater_re.LIW,
      'MF':     self.updater_re.MF,
      'NNZ':    self.updater_re.NNZ,
      'RTOL':   self.updater_re.RTOL,
      'ATOL':   self.updater_re.ATOL}

  def save_common_block(self):
    self.updater_re.save_restore_common_block(job=1)
    
  def restore_common_block(self):
    self.updater_re.save_restore_common_block(job=2)
    
  def update(self, vector[double] y, double t, double dt, int istate=0):
    cdef int i
    cdef double t1

    self.updater_re.set_user_data(self.user_data.ptr)

    for i in range(self.updater_re.NEQ):
      self.user_data.y[i] = y[i]

    if istate != 0:
      self.updater_re.set_ISTATE(istate)
    if self.updater_re.ISTATE not in [0,1]:
      self.restore_common_block()
    if self.updater_re.ISTATE < 0:
      if self.updater_re.ISTATE in [-1, -4, -5]:
        self.updater_re.set_ISTATE(3)
      else:
        print('Unrecoverable error: ISTATE = ', self.updater_re.ISTATE)
        return

    t1 = self.updater_re.update(t, dt, self.user_data.y)

    if istate != 1:
      self.save_common_block()
    return t1, [self.user_data.y[i] for i in range(self.updater_re.NEQ)]

  cdef _get_all_reactions(self):
    return [{'reactants': _.sReactants,
             'products': _.sProducts,
             'abc': _.abc,
             'Trange': _.Trange,
             'itype': _.itype,
             'drdy': _.drdy,
             'rate': _.rate,
             'heat': _.heat
            }
            for _ in self.user_data.reactions]

  def add_reaction(self,
                   vector[string] sReactants,
                   vector[string] sProducts,
                   vector[double] abc,
                   vector[double] Trange,
                   int itype):
    cdef Reaction rs
    rs = Reaction(sReactants, sProducts, abc, Trange, itype)
    self.user_data.add_reaction(rs)

  def clear_reactions(self):
    self.user_data.clear_reactions()

  def set_phy_param(self, string name, double val):
    self.user_data.set_phy_param(name, val)
    self.user_data.physical_params.prep_params()

  def set_phy_param_from_file(self, string fname):
    self.user_data.physical_params.from_file(fname)
    self.user_data.physical_params.prep_params()

  def get_phy_param(self, string name):
    return self.user_data.get_phy_param(name)

  def set_phy_params_by_dict(self, d):
    for k in d:
      self.set_phy_param(k, d[k])
    self.user_data.physical_params.prep_params()

  def get_all_phy_params(self):
    return self.user_data.get_all_phy_params()

  def get_all_reactions(self):
    return self._get_all_reactions()

  def assort_reactions(self):
    return self.user_data.assort_reactions()

  def assignElementsToSpecies(self, elements=None):
    if elements:
      self.user_data.assignElementsToSpecies(elements)
    else:
      self.user_data.assignElementsToSpecies()

  def calculateSpeciesMasses(self, elements=None):
    if elements:
      self.user_data.calculateSpeciesMasses(elements)
    else:
      self.user_data.calculateSpeciesMasses()

  def calculateSpeciesVibFreqs(self):
    self.user_data.calculateSpeciesVibFreqs()

  def calculateSpeciesDiffBarriers(self):
    self.user_data.calculateSpeciesDiffBarriers()

  def calculateSpeciesQuantumMobilities(self):
    self.user_data.calculateSpeciesQuantumMobilities()

  def calculateReactionHeat(self):
    self.user_data.calculateReactionHeat()

  def classifySpeciesByPhase(self):
    self.user_data.classifySpeciesByPhase()

  @property
  def reactions(self):
    return self._get_all_reactions()

  @property
  def reaction_types(self):
    return self.user_data.reaction_types

  @property
  def physical_params(self):
    return self.get_all_phy_params()

  @property
  def name2idx(self):
    return self.user_data.species.name2idx

  @property
  def idx2name(self):
    return self.user_data.species.idx2name

  @property
  def elementsSpecies(self):
    return self.user_data.species.elementsSpecies

  @property
  def massSpecies(self):
    return self.user_data.species.massSpecies

  @property
  def enthalpies(self):
    return self.user_data.species.enthalpies

  @property
  def vibFreqs(self):
    return self.user_data.species.vibFreqs

  @property
  def diffBarriers(self):
    return self.user_data.species.diffBarriers

  @property
  def quantMobilities(self):
    return self.user_data.species.quantMobilities

  @property
  def gasSpecies(self):
    return self.user_data.species.gasSpecies

  @property
  def surfaceSpecies(self):
    return self.user_data.species.surfaceSpecies

  @property
  def mantleSpecies(self):
    return self.user_data.species.mantleSpecies

  @property
  def abundances(self):
    return self.user_data.species.abundances

  def load_reactions(self, fname, nReactants=3, nProducts=4, nABC=3,
    lenSpeciesName=12, lenABC=9, nT=2, lenT=6, lenType=3, rowlen_min=126):
    load_reactions(fname, self.user_data, nReactants, nProducts,
    nABC, lenSpeciesName, lenABC, nT, lenT, lenType, rowlen_min)

  def loadInitialAbundances(self, fname):
    loadInitialAbundances(self.user_data.species, fname)

  def loadSpeciesEnthalpies(self, fname):
    loadSpeciesEnthalpies(self.user_data.species, fname)

  def assignReactionHandlers(self):
    assignReactionHandlers(self.user_data)

  def setAbundanceByName(self, name, val):
    idx = self.user_data.species.name2idx[name]
    self.user_data.species.abundances[idx] = val
