#!/usr/bin/env python3
import numpy as np
from fluxtallies import FluxTallies
import os
from data_nuc import DEPLETION_NUCLIDES

new_pwd = os.getcwd()

nuclide_depletion_tally = DEPLETION_NUCLIDES
score_depletion_tally = ['fission','(n,2n)','(n,3n)','(n,4n)','(n,p)','(n,a)','(n,gamma)'] #['(n,gamma)','fission','(n,2n)']

a = {}

a['tally1'] = {}

a['tally1']['id'] = '1'

a['tally1']['name'] = 'ngamma-depletion-tally'

a['tally1']['filter'] = ['1']

a['tally1']['nuclides'] = nuclide_depletion_tally

a['tally1']['scores'] = score_depletion_tally

'''
a['tally2'] = {}

a['tally2']['id'] = '2'

a['tally2']['name'] = 'nfission-depletion-tally'

a['tally2']['filter'] = ['1']

a['tally2']['nuclides'] = ['U236','Pu240','Cm242','Cf252'] 

a['tally2']['scores'] = ['fission']

a['tally3'] = {}

a['tally3']['id'] = '3'

a['tally3']['name'] = 'np-depletion-tally'

a['tally3']['filter'] = ['1']

a['tally3']['nuclides'] = ['Co58','Re185','Ta180'] 

a['tally3']['scores'] = ['(n,p)']
'''

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

#bps_matrix = np.matrix([[1.E-5,0.],[1.E-1,200.],[1.E+2,400.],[1.E+4,800.],[1.E+6,100.],[8.E+6,450.],[8.12E+6,450.]])
#bps_matrix = np.matrix([[1.E-5,0.],[1.,100.],[1.E+2,1550.],[1.E+3,350.],[6.E+4,20],[1.E+5,600],[3.E+5,50],[1.E+6,1200],[2.E+6,100.],[8.E+6,330.],[8.12E+6,700.]])
#bps_matrix = np.matrix([[1.E-5,0.],[1.,100.],[1.E+2,550.],[1.E+3,350.],[6.E+4,20],[1.E+5,100],[3.E+5,50],[1.E+6,200],[2.E+6,100.],[8.E+6,330.],[8.12E+6,700.]])
#bps_matrix = np.matrix([[1.E-5,0.],[5.E+4,2500.],[8.11E+6,10000.]])
#bps_matrix = np.matrix([[1.E-5,0.],[1.,201.],[1.E+2,201.],[1.E+3,201.],[1.E+4,201.],[1.E+5,201.],[1.E+6,201.],[8.E+6,301.],[8.12E+6,1001.]])
#bps_matrix = np.matrix([[1.0000000000000000E-5,0.],[1.0000000000000000E+2,400.]])#,[1.E+3,200],[1.E+4,200],[1.E+5,200.],[1.E+6,200.],[8.E+6,300.],[8.12E+6,1000.]])

#bps_matrix = np.matrix([[2.E+6,0.],[4.0E+6,108.],[6.E+6,51.],[8.11E+6,101.]])
#bps_matrix = np.matrix([[1.E-5,0.],[1.E+6,101.],[2.E+6,51.],[4.E+6,201.],[6.E+6,51.],[8.11E+6,201.]])
bps_matrix = np.matrix([[1.E-5,0.],[1.E+6,101.],[3.E+6,51.],[8.11E+6,551.],[2.E+7,26]])


prueba = FluxTallies(nuclide_list=nuclide_depletion_tally,score_list=score_depletion_tally,tallies_dict=a)

prueba.path_to_lib = new_pwd

prueba.total_groups = 726  #362 + 250#2500

prueba.lower_lim = 1.E-5
prueba.epth_lim = 1.0
prueba.res_lim = 4.E+6
prueba.upper_lim = 9.E+6

prueba.generate_tallies_xml('tallies.xml',bps_matrix,SHEM361_only=False)

#prueba.plot_xs(bps_matrix)

#prueba.download_endf() 

#prueba.run_njoy_arbitrary_struc(bps_matrix,SHEM361_only=False)

#prueba.parse_groupr()

#prueba.collapse_from_statepoint(sp)

#prueba.check_accuracy_from_summary(reaction_rate_filename='reaction-rates-sort.csv')
