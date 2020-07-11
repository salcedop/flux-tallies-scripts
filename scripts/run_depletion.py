#!/usr/bin/env python3
import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import os
from fluxtallies import FluxTallies
import os
from energy_struc import energy_struc500
from data_nuc import DEPLETION_NUCLIDES,O16Peaks
###############################################################################
#                      Simulation Input File Parameters
###############################################################################

# OpenMC simulation parameters
xs_path1  = "/gpfs/jlse-fs0/users/salcedop/data/nndc_hdf5/cross_sections.xml"
my_vars = {'OPENMC_CROSS_SECTIONS':xs_path1}
os.environ.update(my_vars)

xs_path2  = "/home/salcedop/lib-piecewise-hybrid-updated-500/second-test/cross_sections_hybrids.xml"
my_vars2 = {'OPENMC_HYBRID_CROSS_SECTIONS':xs_path2}
os.environ.update(my_vars2)
batches = 70
inactive = 15
height = 14.4
particles = 500000

temperature = 300
tolerance = temperature - 290
# Depletion simulation parameters

chain_file = './chain_endfb71.xml'

power = 174 # W/cm, for 2D simulations only (use W for 3D)

'''
time_steps = np.array([0.25, 6.00, 6.25, 12.50, 25.00, 25.00, 25.00, # 100.00 EFPD 
             25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, # 300.00 EFPD
             62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50  # 1500.0 EFPD
             ])*24*60*60 
'''

time_steps = np.array([0.25, 6.00, 6.25, 12.50, 25.00, 25.00, 25.00, # 100.00 EFPD 
             25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, # 300.00 EFPD
             25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, 25.00, # 500.00 EFPD
             62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50, # 1000.0 EFPD
             62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50, 62.50  # 1500.0 EFPD
             ])*24*60*60

###############################################################################
#                              Define materials
###############################################################################

# Instantiate some Materials and register the appropriate Nuclides
uo2 = openmc.Material(material_id=1, name='UO2 fuel at 4.5% wt enrichment')
uo2.set_density('g/cm3', 10.29769)
uo2.add_element('U', 1., enrichment=4.5)
uo2.add_element('O', 2.)
uo2.depletable = True
uo2.temperature = temperature


helium = openmc.Material(material_id=2, name='Helium for gap')
helium.set_density('g/cm3', 0.001598)
helium.add_element('He', 2.4044e-4)
helium.temperature = temperature

zircaloy = openmc.Material(material_id=3, name='Zircaloy 4')
zircaloy.set_density('g/cm3', 6.55)
zircaloy.add_element('Sn', 0.014  , 'wo')
zircaloy.add_element('Fe', 0.00165, 'wo')
zircaloy.add_element('Cr', 0.001  , 'wo')
zircaloy.add_element('Zr', 0.98335, 'wo')
zircaloy.temperature = temperature

borated_water = openmc.Material(material_id=4, name='Borated water')
borated_water.set_density('g/cm3', 0.740582)
borated_water.add_element('B', 4.0e-5)
borated_water.add_element('H', 5.0e-2)
borated_water.add_element('O', 2.4e-2)
borated_water.add_s_alpha_beta('c_H_in_H2O')
borated_water.temperature = temperature
###############################################################################
#                             Create geometry
###############################################################################

# Instantiate ZCylinder surfaces
fuel_or = openmc.ZCylinder(surface_id=1, x0=0, y0=0, R=0.39218, name='Fuel OR')
clad_ir = openmc.ZCylinder(surface_id=2, x0=0, y0=0, R=0.40005, name='Clad IR')
clad_or = openmc.ZCylinder(surface_id=3, x0=0, y0=0, R=0.45720, name='Clad OR')
left = openmc.XPlane(surface_id=4, x0=-0.62992, name='left')
right = openmc.XPlane(surface_id=5, x0=0.62992, name='right')
bottom = openmc.YPlane(surface_id=6, y0=-0.62992, name='bottom')
top = openmc.YPlane(surface_id=7, y0=0.62992, name='top')

left.boundary_type = 'reflective'
right.boundary_type = 'reflective'
top.boundary_type = 'reflective'
bottom.boundary_type = 'reflective'

# Instantiate Cells
fuel = openmc.Cell(cell_id=1, name='cell 1')
gap = openmc.Cell(cell_id=2, name='cell 2')
clad = openmc.Cell(cell_id=3, name='cell 3')
water = openmc.Cell(cell_id=4, name='cell 4')

# Use surface half-spaces to define regions
fuel.region = -fuel_or
gap.region = +fuel_or & -clad_ir
clad.region = +clad_ir & -clad_or
water.region = +clad_or & +left & -right & +bottom & -top

# Register Materials with Cells
fuel.fill = uo2
gap.fill = helium
clad.fill = zircaloy
water.fill = borated_water

# Instantiate Universe
root = openmc.Universe(universe_id=0, name='root universe')

# Register Cells with Universe
root.add_cells([fuel, gap, clad, water])

# Instantiate a Geometry, register the root Universe
geometry = openmc.Geometry(root)

###############################################################################
#                     Set volumes of depletable materials
###############################################################################

# Compute cell areas
area = {}
area[fuel] = np.pi * fuel_or.r** 2

# Set materials volume for depletion. Set to an area for 2D simulations
uo2.volume = area[fuel]

###############################################################################
#                     Transport calculation settings
###############################################################################

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.hybrid = True
# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-0.62992, -0.62992, -1, 1,0.62992, 0.62992]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.source.Source(space=uniform_dist)

#settings_file.temperature = {'multipole': True, 'tolerance': tolerance,'default':temperature}

#settings_file.temperature = {'multipole': True, 'tolerance': tolerance,'default':temperature}
settings_file.output = {'tallies': False,'summary':False}

###############################################################################
#                   Initialize and run depletion calculation
###############################################################################

op = openmc.deplete.Operator(geometry,settings_file,hybrid=True,chain_file=chain_file,Energy_Struc=energy_struc500)

# Perform simulation using the predictor algorithm
integrator = openmc.deplete.PredictorIntegrator(op, time_steps, power)
integrator.integrate()

###############################################################################
#                    Read depletion calculation results
###############################################################################

# Open results file
results = openmc.deplete.ResultsList("depletion_results.h5")
'''
# Obtain K_eff as a function of time
time, keff = results.get_eigenvalue()

# Obtain U235 concentration as a function of time
time, n_U235 = results.get_atoms('1', 'U235')

# Obtain Xe135 absorption as a function of time
time, Xe_gam = results.get_reaction_rate('1', 'Xe135', '(n,gamma)')

###############################################################################
#                            Generate plots
###############################################################################

plt.figure()
plt.plot(time/(24*60*60), keff, label="K-effective")
plt.xlabel("Time (days)")
plt.ylabel("Keff")
plt.show()

plt.figure()
plt.plot(time/(24*60*60), n_U235, label="U 235")
plt.xlabel("Time (days)")
plt.ylabel("n U5 (-)")
plt.show()

plt.figure()
plt.plot(time/(24*60*60), Xe_gam, label="Xe135 absorption")
plt.xlabel("Time (days)")
plt.ylabel("RR (-)")
plt.show()
plt.close('all')
'''
