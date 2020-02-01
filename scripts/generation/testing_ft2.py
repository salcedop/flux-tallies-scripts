#!/usr/bin/env python3
import numpy as np
from fluxtallies import FluxTallies
import os
from data_nuc import DEPLETION_NUCLIDES,mypeaks
import openmc
import openmc.data

new_pwd = os.getcwd()

DEPLETION_NUCLIDES = ['H1', 'H2', 'H3', 'He3', 'He4', 'Li6', 'Li7']

nuclide_depletion_tally = DEPLETION_NUCLIDES
score_depletion_tally = ['fission','(n,2n)','(n,3n)','(n,4n)','(n,p)','(n,a)','(n,gamma)'] #['(n,gamma)','fission','(n,2n)']

a = {}

a['tally1'] = {}

a['tally1']['id'] = '1'

a['tally1']['name'] = 'ngamma-depletion-tally'

a['tally1']['filter'] = ['1']

a['tally1']['nuclides'] = nuclide_depletion_tally

a['tally1']['scores'] = score_depletion_tally


a['tally4'] = {}

a['tally4']['id'] = '2'

a['tally4']['name'] = 'single-group-flux'

a['tally4']['filter'] = ['1']

a['tally4']['nuclides'] = ''

a['tally4']['scores'] = ['flux']

a['tally5'] = {}

a['tally5']['id'] = '3'

a['tally5']['name'] = 'flux_tally'

a['tally5']['scores'] = ['group-flux']

a['tally5']['nuclides'] = ''

a['tally5']['filter'] = ['1','2']

b = {}

b[0] = {'id':'1','name':'reaction-rate','filter':['1'],'nuclides':nuclide_depletion_tally,'scores':score_depletion_tally}
#b[1] = {'id':'2','name':'single-group-flux','filter':['1'],'nuclides':'','scores':['flux']}
b[2] = {'id':'2','name':'flux-tally','filter':['1','2'],'nuclides':'','scores':['group-flux']}

#bps_matrix = np.matrix([[1.E-5,0.],[1.E+6,101.],[3.E+6,51.],[8.11E+6,551.],[2.E+7,26]])

#nuc = openmc.data.IncidentNeutron.from_hdf5('O16.h5')

#peaks,_ = find_peaks(nuc[107].xs['294K'].y)

#peak_energy = nuc[107].xs['294K'].x[mypeaks][0:27]

bps = np.ones(len(mypeaks)) * 11.

bps[0] = 11.

bps_matrix = np.matrix([[1.E-5,0.],[1.E+6,51.]])

last_row = [[2.E+7,21.]]

for irow,row in enumerate(mypeaks):
  new_row = [[row,bps[irow]]]
  bps_matrix = np.concatenate((bps_matrix,new_row))

bps_matrix = np.concatenate((bps_matrix,last_row))

#print(bps_matrix)

len_bps = (len(bps_matrix))

#print(len_bps)

#bps_matrix = np.matrix([[1.E-5,0.],[4.E+4,51.],[1.1E+6,[3.E+6,51.],[8.11E+6,551.],[2.E+7,26]])
pa = 'tests-results'
prueba = FluxTallies(nuclide_list=nuclide_depletion_tally,score_list=score_depletion_tally,tallies_dict=b,path_to_lib=pa)

#prueba.path_to_lib = 'shem-prac'

prueba.total_groups = 320


prueba.temperatures = [300.,1000.]
#prueba.background_xs = [1.E+10,1000.]

prueba.generate_tallies_xml('/home/salcedop/scripts/tallies.xml',bps_matrix,default_group_struc=False)#shem_path='/home/salcedop/scripts/SHEM361.hdf5')

#prueba.plot_xs(bps_matrix)

prueba.download_endf() 

prueba.run_njoy_arbitrary_struc(bps_matrix,default_group_struc=True)

#prueba.parse_groupr()

#prueba.generating_cross_section_xml()

#prueba.collapse_from_statepoint(sp)

#prueba.check_accuracy_from_summary(reaction_rate_filename='reaction-rates-sort.csv')
