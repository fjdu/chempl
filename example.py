import chempl
from chempl import (cs, run_one_model, get_total_charge,
                    get_phy_params_default, N_H_to_Av, Av_to_N_H,
                    N_H_to_ngas, Td_from_Av, chem2tex, printFormationDestruction)

model = chempl.ChemModel(
    fReactions=b'rate12_noSurf_only_H2_formation.dat',
    fInitialAbundances=b'initial_abundances.dat',
    phy_params=
    {b'Av': 10.0,
     b'G0_UV': 1.0,
     b'Ncol_H2': 1e22,
     b'T_dust': 10.0,
     b'T_gas': 10.0,
     b'chemdesorption_factor': 0.05,
     b'chi_cosmicray': 1.0,
     b'dust2gas_mass': 1e-2,
     b'dust_albedo': 0.6,
     b'dust_material_density': 2.0,
     b'dust_radius': 0.1e-4,
     b'dust_site_density': 1e15,
     b'dv_km_s': 1.0,
     b'mean_mol_weight': 1.4,
     b'n_gas': 1e5}).prepare().set_solver()

print('Number of reactions', len(model.reactions))
print('Number of species', len(model.name2idx))

t = 0.0 # Initial time, in seconds
dt = 0.1 # Initial time step
tRatio = 1.1 # Ratio between time steps
tMax = 1e6 * cs.phy_SecondsPerYear # Targeted final time
nMax = 2000 # Number of time steps (estimated)
store = {'t': [], 'y': []} # Storage for the intermediate steps
y = model.abundances # Initial abundances as loaded in the setup

for i in range(nMax):
    t, y = model.update(y, t=t, dt=dt)
    store['t'].append(t)
    store['y'].append(y)
    if t >= tMax:
        break
    dt *= tRatio
    if t+dt > tMax:
        dt = tMax - t

'The actual number of time steps is:', len(store['t'])
