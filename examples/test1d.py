#!/usr/bin/env python

import cantera as ct
import ctmicro as cm

tinlet = 300.
p = ct.one_atm
mdot = .5
reactants = 'H2:.5, O2:1, AR:7'

profile = cm.ErfProfile((0, .1), (300., 1300.), .04, .02)

tol_ss = [1.0e-9, 1.0e-13]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-8, 1.0e-13]  # [rtol atol] for time stepping
loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = 1  # 1 to enable refinement, 0 to disable

gas = ct.Solution('h2o2.yaml')
gas.TPX = tinlet, p, reactants

f = cm.ChannelFlame(gas, profile)

# inlet
f.inlet.T = tinlet
f.inlet.X = reactants
f.inlet.mdot = mdot

# channel
f.flame.Nu = 4.36  # approximation
f.flame.diameter = .002  # 1mm tube

print(f.grid)
print(f.flame.get_wall_profile())

f.set_initial_guess()
# print(f.Tw())
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)
f.show_solution()

f.energy_enabled = False
f.set_max_jac_age(16, 16)
f.set_time_step(1.e-7, [1, 2, 4, 8, 16])
f.set_max_time_step(1.e-3)
f.solve(loglevel, refine_grid=False)
f.save('lossyflame_test.xml', 'no_energy',
       'solution with the energy equation disabled')

f.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.12)
# f.set_refine_criteria(ratio=2.0, slope=0.04, curve=0.08)
# f.set_refine_criteria(ratio=4.0, slope=0.08, curve=0.16)
# f.set_refine_criteria(ratio=4.5, slope=0.09, curve=0.18)
f.energy_enabled = True
f.solve(loglevel, refine_grid)
f.save('test1d.xml', 'energy',
       'solution with the energy equation enabled')

f.write_csv('test1d.csv', quiet=False)
