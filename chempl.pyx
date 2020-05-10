# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cppmap
from libcpp.set cimport set as cppset
from libcpp.utility cimport pair
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
    void allocate_abundances()

  cdef cppclass AuxData:
    pass

  cdef cppclass PhyParams:
    void prep_params()
    int from_file(string fname)
    void add_a_timedependency(string name,
        vector[double] ts, vector[double] vs)
    void remove_a_timedependency(string name)
    vector[string] get_timeDependency_names()

  ctypedef vector[Reaction] Reactions
  ctypedef cppmap[string, double] Elements
  ctypedef cppmap[int, int] ReactionTypes
  ctypedef double (*RateCalculator)(const double&, double *,
           Reaction&, const PhyParams&, const Species&, AuxData&)
  ctypedef cppmap[int, RateCalculator] RateCalculators

  cdef cppclass Chem_data:
    void add_reaction(Reaction& rs)
    void modify_reaction(const int& iReact, const cppmap[string, vector[double]] &par)
    void clear_reactions()
    void find_duplicate_reactions()
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
    double calculate_a_rate(double t, double *y, Reaction& r, int updatePhyParams)
    vector[pair[int, double]] getFormationReactionsWithRates(int iSpecies, double t, double*y)
    vector[pair[int, double]] getDestructionReactionsWithRates(int iSpecies, double t, double*y)
    Reactions reactions
    PhyParams physical_params
    Species species
    ReactionTypes reaction_types
    Chem_data* ptr
    RateCalculators rate_calculators
    vector[int] dupli
    double* y

  void update_phy_params(double t, PhyParams& p)
  cppmap[string, int] assignElementsToOneSpecies(string name, const Elements& elements)


cdef extern from "utils.hpp" namespace "UTILS":
  double interpol(const vector[double]& ts, const vector[double]& vs, double t)


cdef extern from "logistics.hpp" namespace "LOGIS":
  void load_reactions(const string& fname, Chem_data& cdata,
    int nReactants, int nProducts, int nABC, int lenSpeciesName,
    int lenABC, int nT, int lenT, int lenType, int rowlen_min)
  int loadInitialAbundances(Species& species, string fname)
  int loadSpeciesEnthalpies(Species& species, string fname)


cdef extern from "calculate_reaction_rate.hpp" namespace "CALC_RATE":
  void assignReactionHandlers(Chem_data&)
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
  double arrhenius(const double &T, const vector[double] &abc,
                       const int &iS)

cdef extern from "rate_equation_lsode.hpp" namespace "RATE_EQ":
  cdef cppclass Updater_RE:
    void set_user_data(Chem_data* udata)
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


cdef class ChemModel:

  cdef Chem_data cdata
  cdef Updater_RE updater_re
  cdef public all_reactions

  def set_solver(self, rtol=1e-6, atol=1e-30, mf=21, LRW_F=6,
                 showmsg=1, msglun=6, solver_id=0):
    self.updater_re.set_user_data(self.cdata.ptr)
    self.updater_re.set_sparse()
    self.updater_re.initialize_solver(rtol, atol, mf, LRW_F, solver_id)
    self.updater_re.set_solver_msg(showmsg);
    self.updater_re.set_solver_msg_lun(msglun);
    self.updater_re.allocate_rsav_isav()
    self.cdata.allocate_y()

  def allocate_y(self):
    self.cdata.allocate_y()

  def deallocate_y(self):
    self.cdata.deallocate_y()

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
    
  def update(self, vector[double] y, double t, double dt, int istate=0, interruptMode=False):
    cdef int i
    cdef double t1

    self.updater_re.set_user_data(self.cdata.ptr)

    for i in range(self.updater_re.NEQ):
      self.cdata.y[i] = y[i]

    if istate != 0:
      self.updater_re.set_ISTATE(istate)
    if interruptMode and self.updater_re.ISTATE not in [0,1]:
      self.restore_common_block()
    if self.updater_re.ISTATE < 0:
      if self.updater_re.ISTATE in [-1, -4, -5]:
        self.updater_re.set_ISTATE(3)
      else:
        print('Unrecoverable error: ISTATE = ', self.updater_re.ISTATE)
        return

    t1 = self.updater_re.update(t, dt, self.cdata.y)

    if interruptMode and istate != 1:
      self.save_common_block()
    return t1, [self.cdata.y[i] for i in range(self.updater_re.NEQ)]

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
            for _ in self.cdata.reactions]

  def add_reaction(self,
                   vector[string] sReactants,
                   vector[string] sProducts,
                   vector[double] abc,
                   vector[double] Trange,
                   int itype):
    cdef Reaction rs
    rs = Reaction(sReactants, sProducts, abc, Trange, itype)
    self.cdata.add_reaction(rs)

  def modify_reaction(self, const int& iReact, const cppmap[string, vector[double]] &par):
    """modify_reaction(iReact, map[string, vector[double]])
    string: b"abc" or b"Trange"
    """
    self.cdata.modify_reaction(iReact, par)

  def add_reaction_by_dict(self, r):
    cdef Reaction rs
    rs = Reaction(r['reactants'], r['products'], r['abc'], r['Trange'], r['itype'])
    self.cdata.add_reaction(rs)

  def clear_reactions(self):
    self.cdata.clear_reactions()

  def find_duplicate_reactions(self):
    self.cdata.find_duplicate_reactions()

  def calculate_a_rate(self, double t, vector[double] y, int iReac, int updatePhyParams=False):
    """calculate_a_rate(t, y, iReac, updatePhyParams=False)"""
    for i in range(len(y)):
      self.cdata.y[i] = y[i]
    return self.cdata.calculate_a_rate(t, self.cdata.y, self.cdata.reactions[iReac], updatePhyParams)

  def getFormationReactionsWithRates(self, int iSpecies, double t, vector[double] y):
    for i in range(len(y)):
      self.cdata.y[i] = y[i]
    return self.cdata.getFormationReactionsWithRates(iSpecies, t, self.cdata.y)

  def getDestructionReactionsWithRates(self, int iSpecies, double t, vector[double] y):
    for i in range(len(y)):
      self.cdata.y[i] = y[i]
    return self.cdata.getDestructionReactionsWithRates(iSpecies, t, self.cdata.y)

  def set_phy_param(self, string name, double val):
    self.cdata.set_phy_param(name, val)
    self.cdata.physical_params.prep_params()

  def set_phy_param_from_file(self, string fname):
    self.cdata.physical_params.from_file(fname)
    self.cdata.physical_params.prep_params()

  def get_phy_param(self, string name):
    return self.cdata.get_phy_param(name)

  def set_phy_params_by_dict(self, d):
    for k in d:
      self.set_phy_param(k, d[k])
    self.cdata.physical_params.prep_params()

  def get_all_phy_params(self):
    return self.cdata.get_all_phy_params()

  def add_time_dependency(self, name, ts, vs):
    """add_phy_param_time_dependency(name, ts, vs)"""
    self.cdata.physical_params.add_a_timedependency(name, ts, vs)

  def remove_time_dependency(self, name):
    """remove_phy_param_time_dependency(name)"""
    self.cdata.physical_params.remove_a_timedependency(name)

  def get_time_dependency_names(self):
    """get_timeDependency_names()"""
    return self.cdata.physical_params.get_timeDependency_names()

  def update_phy_params(self, t):
    update_phy_params(t, self.cdata.physical_params)
    
  def assort_reactions(self):
    return self.cdata.assort_reactions()

  def assignElementsToSpecies(self, elements=None):
    if elements:
      self.cdata.assignElementsToSpecies(elements)
    else:
      self.cdata.assignElementsToSpecies()

  def calculateSpeciesMasses(self, elements=None):
    if elements:
      self.cdata.calculateSpeciesMasses(elements)
    else:
      self.cdata.calculateSpeciesMasses()

  def calculateSpeciesVibFreqs(self):
    self.cdata.calculateSpeciesVibFreqs()

  def calculateSpeciesDiffBarriers(self):
    self.cdata.calculateSpeciesDiffBarriers()

  def calculateSpeciesQuantumMobilities(self):
    self.cdata.calculateSpeciesQuantumMobilities()

  def calculateReactionHeat(self):
    self.cdata.calculateReactionHeat()

  def classifySpeciesByPhase(self):
    self.cdata.classifySpeciesByPhase()

  def get_all_reactions(self):
    self.all_reactions = self._get_all_reactions()
    return self.all_reactions

  @property
  def reactions(self):
    if self.all_reactions is not None:
      return self.all_reactions
    self.all_reactions = self._get_all_reactions()
    return self.all_reactions

  @property
  def reaction_types(self):
    return self.cdata.reaction_types

  @property
  def physical_params(self):
    return self.get_all_phy_params()

  @property
  def name2idx(self):
    return self.cdata.species.name2idx

  @property
  def idx2name(self):
    return self.cdata.species.idx2name

  @property
  def elementsSpecies(self):
    return self.cdata.species.elementsSpecies

  @property
  def massSpecies(self):
    return self.cdata.species.massSpecies

  @property
  def enthalpies(self):
    return self.cdata.species.enthalpies

  @property
  def vibFreqs(self):
    return self.cdata.species.vibFreqs

  @property
  def diffBarriers(self):
    return self.cdata.species.diffBarriers

  @property
  def quantMobilities(self):
    return self.cdata.species.quantMobilities

  @property
  def gasSpecies(self):
    return self.cdata.species.gasSpecies

  @property
  def surfaceSpecies(self):
    return self.cdata.species.surfaceSpecies

  @property
  def mantleSpecies(self):
    return self.cdata.species.mantleSpecies

  @property
  def abundances(self):
    return self.cdata.species.abundances

  @property
  def duplicate_reactions(self):
    return self.cdata.dupli

  def load_reactions(self, fname, nReactants=3, nProducts=4, nABC=3,
    lenSpeciesName=12, lenABC=9, nT=2, lenT=6, lenType=3, rowlen_min=126):
    load_reactions(fname, self.cdata, nReactants, nProducts,
    nABC, lenSpeciesName, lenABC, nT, lenT, lenType, rowlen_min)

  def loadInitialAbundances(self, fname):
    loadInitialAbundances(self.cdata.species, fname)

  def loadSpeciesEnthalpies(self, fname):
    loadSpeciesEnthalpies(self.cdata.species, fname)

  def assignReactionHandlers(self):
    assignReactionHandlers(self.cdata)

  def setAbundanceByName(self, name, val):
    if len(self.cdata.species.abundances) == 0:
      self.cdata.species.allocate_abundances()
    idx = self.cdata.species.name2idx[name]
    self.cdata.species.abundances[idx] = val

  def setAbundanceByDict(self, a):
    if len(self.cdata.species.abundances) == 0:
      self.cdata.species.allocate_abundances()
    for k in a:
      idx = self.cdata.species.name2idx[k]
      self.cdata.species.abundances[idx] = a[k]

  def setAbundances(self, vals):
    if len(self.cdata.species.abundances) == 0:
      self.cdata.species.allocate_abundances()
    for i in range(len(self.cdata.species.idx2name)):
      self.cdata.species.abundances[i] = vals[i]

  def assignElementsToOneSpecies(self, name, elements):
    return assignElementsToOneSpecies(name, elements)

  def __init__(self, fReactions=None, fInitialAbundances=None,
               fSpeciesEnthalpies=None):
    """
  __init__(self, fReactions=None, fInitialAbundances=None,
           fSpeciesEnthalpies=None)
    """
    self.all_reactions = None

    if fReactions is not None:
      self.load_reactions(fReactions)
    if fInitialAbundances is not None:
      self.loadInitialAbundances(fInitialAbundances)
    if fSpeciesEnthalpies is not None:
      self.loadSpeciesEnthalpies(fSpeciesEnthalpies)

  def prepare(self):
    self.assort_reactions()
    self.assignElementsToSpecies()
    self.calculateSpeciesMasses()
    self.calculateSpeciesVibFreqs()
    self.calculateSpeciesDiffBarriers()
    self.calculateSpeciesQuantumMobilities()
    self.calculateReactionHeat()
    self.classifySpeciesByPhase()
    self.assignReactionHandlers()


def simpleInterpol(ts, vs, t):
    """simpleInterpol(ts, vs, t)
    ts: t values
    vs: values
    t: t value to be interpolated"""
    return interpol(ts, vs, t)

def rate_Arrhenius(T, abc, iS=0):
    return arrhenius(T, abc, iS)
