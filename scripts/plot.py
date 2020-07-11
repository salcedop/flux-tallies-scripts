#!/usr/bin/env python3
from pathlib import Path
import urllib.request
import re
import os
import openmc
#import openmc.data.njoy
#import openmc.data.njoy_groupr
import numpy as np
import h5py

try:
  from openmc.clean_xml import clean_xml_indentation
  old = True
except:
  from openmc._xml import clean_indentation
  old = False

import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

nuclide = 'Co58'
#ranges = {107:list(range(800,850))}
#map_rxn = {107:'(n,a)'}

ranges = {103:list(range(600,650))}
map_rxn = {103:'(n,p)'}

pattern= re.compile("([a-zA-Z]+)([0-9]+)")
components = pattern.match(nuclide)
isotope_name = components.group(1)
isotope_A = components.group(2)
length_A = len(isotope_A)

if (length_A < 3) and (length_A > 1):
  isotope_A = '0' + isotope_A
elif(length_A < 2):
  isotope_A = '00' + isotope_A

filename = '/home/salcedop/openmc/data/nndc/293.6K/{}_{}_293.6K.ace'.format(isotope_name,isotope_A)
inuc_endf=openmc.data.IncidentNeutron.from_ace(filename)

h5file = h5py.File('{}.h5'.format(nuclide),'r')
group = list(h5file.values())[0]
plt.figure()
for irxn,items in ranges.items():
 labels = []
 full_list = [irxn]
 full_list += items
 for mts in full_list:
  if (mts > 110):
   continue
  try:
   if (mts == irxn): 
      energy_values = inuc_endf[ranges[irxn][0]].xs['294K'].x
   else:
      energy_values = inuc_endf[mts].xs['294K'].x
   group_nuc = group['reactions/reaction_'+str(mts)]
   xs = group_nuc['294K/xs'].value
   plt.loglog(energy_values,xs)
   labels.append('{}-reconstructed'.format(mts))
  except:
   continue
 '''
 h5file2 = h5py.File('/home/salcedop/library-nuclide-bins-2500-groups-422-nucs-SHEM/{}/{}.h5'.format(nuclide,nuclide),'r')
 #group = list(h5file.values())[0]
 group_nuc = h5file2['{}/reactions/{}'.format(nuclide,map_rxn[irxn])]
 xs = group_nuc['groupr'][:]
 plt.loglog(inuc_endf[ranges[irxn][0]].xs['294K'].x,xs)
 labels.append('MG-{}'.format(irxn))
 '''
 
 energy_values = inuc_endf[irxn].xs['294K'].x
 xs = inuc_endf[irxn].xs['294K'].y
 plt.loglog(energy_values,xs)
 labels.append('{}-ACE'.format(irxn))
 
 plt.xlabel('Energy (eV)')
 plt.ylabel('Cross sections (b)')
 plt.title('{}-{} Cross-section versus Energy'.format(nuclide,irxn))
 plt.legend(labels,ncol=2,loc='upper left',prop={'size':8})

 plt.rcParams['axes.labelsize'] = 16

 plt.rcParams['axes.xmargin'] = 0

 plt.rcParams['axes.ymargin'] = 0

 plt.rcParams['font.family'] = 'serif'

 plt.rcParams['font.serif'] = ['Times']

 plt.rcParams['font.size'] = 14

 plt.rcParams['mathtext.fontset'] = 'stix'

 plt.rcParams['pdf.use14corefonts'] = True

 plt.rcParams['savefig.bbox'] = 'tight'

 #plt.rcParams['text.usetex'] = True

 plt.savefig('{}_{}.pdf'.format(nuclide,irxn))
 #plt.show()
