#!/usr/bin/env python3
import numpy as np
import openmc
import openmc.capi
import time
import pandas as pd
import os
import h5py
from mpi4py import MPI

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def get_fuel_mats():
  fuel_mats = []
  mats = openmc.capi.materials
  for key,value in mats.items():
    imat = openmc.capi.materials[key]
    nucs = imat.nuclides
    len_nucs = len(nucs)
    if (len_nucs > 200):
      fuel_mats.append(key)
  return fuel_mats


def prep_data(item):
        nuc_data[item] = {}
        filename = '{}.h5'.format(item,1500)
        path_to_xs = '/home/salcedop/groupr_xs/quick_prac/updating_file/smr-prac/dir/' + item +'/'+filename
        h5file = h5py.File(path_to_xs,'r')
        group = list(h5file.values())[0]
        for irx,rx in enumerate(depletion_rx_list):
          nuc_data[item][rx] = {}
          index = 'reactions/'+rx
          igroup = group[index]
          nuc_data[item][rx]['start_point'] = igroup.attrs['start_point']
          nuc_data[item][rx]['xs'] = igroup['groupr']
        #h5file.close()

def collapse(item):
      for inuc,nucs in enumerate(DEPLETION_NUCLIDES):
        res_dict[item][nucs] = {}
        dens = dens_data[item][nucs]['dens']
        for irx,rx in enumerate(depletion_rx_list):
          running_sum = 0.
          xs = nuc_data[nucs][rx]['xs'] * dens
          threshold = nuc_data[nucs][rx]['start_point']
          end = len(xs)
          for ielement in range(end):
            if (threshold == 0):
              group_xs = 0.
              ielement = 1
            else:
              group_xs = xs[ielement]
            running_sum += flux_tally[ielement+threshold+shift-1]*group_xs

          if (rx == '(n,fission)'):
            rx =  'fission'
          if (rx == '(n,g)'):
            rx = '(n,gamma)'
          res_dict[item][nucs][rx] = running_sum

DEPLETION_NUCLIDES = ['O16', 'O17', 'U234', 'U235', 'U238', 'U236', 'U239', 'U240', 'Np234', 'Np235', 'Np236', 'Np237', 'Np238', 'Np239', 
'Pu236', 'Pu237', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'B11', 'N14', 'N15', 'Fe57', 'Fe58', 'Co59', 'Ni60', 'Ni61', 'Ni62', 'Cu63', 
'Ni64', 'Zn64', 'Cu65', 'Zn65', 'Zn66', 'Zn67', 'Zn68', 'Ga69', 'Zn70', 'Ge70', 'Ga71', 'Ge72', 'Ge73', 'Ge74', 'As74', 'Se74', 'As75', 
'Ge76', 'Se76', 'Se77', 'Se78', 'Se79', 'Br79', 'Se80', 'Kr80', 'Br81', 'Se82', 'Kr82', 'Kr83', 'Kr84', 'Sr84', 'Kr85', 'Rb85', 'Kr86', 
'Rb86', 'Sr86', 'Rb87', 'Sr87', 'Sr88', 'Sr89', 'Y89', 'Sr90', 'Y90', 'Zr90', 'Y91', 'Zr91', 'Zr92', 'Zr93', 'Nb93', 'Zr94', 'Nb94', 
'Mo94', 'Zr95', 'Nb95', 'Mo95', 'Zr96', 'Mo96', 'Mo97', 'Mo98', 'Ru98', 'Mo99', 'Tc99', 'Ru99', 'Mo100', 'Ru100', 'Ru101', 'Ru102', 
'Pd102', 'Ru103', 'Rh103', 'Ru104', 'Pd104', 'Ru105', 'Rh105', 'Pd105', 'Ru106', 'Pd106', 'Pd107', 'Ag107', 'Pd108', 'Cd108', 'Ag109', 
'Pd110', 'Cd110', 'Ag111', 'Cd111', 'Cd112', 'Sn112', 'Cd113', 'In113', 'Sn113', 'Cd114', 'Sn114', 'In115', 'Sn115', 'Cd116', 'Sn116', 
'Sn117', 'Sn118', 'Sn119', 'Sn120', 'Te120', 'Sb121', 'Sn122', 'Te122', 'Sn123', 'Sb123', 'Te123', 'Sn124', 'Sb124', 'Te124', 'Sn125', 
'Sb125', 'Te125', 'Sn126', 'Sb126', 'Te126', 'Xe126', 'I127', 'Te128', 'Xe128', 'I129', 'Xe129', 'Te130', 'I130', 'Xe130', 
'I131', 'Xe131', 'Te132', 'Xe132', 'Ba132', 'Xe133', 'Cs133', 'Ba133', 'Xe134', 'Cs134', 'Ba134', 'I135', 'Xe135', 'Cs135', 
'Ba135', 'Xe136', 'Cs136', 'Ba136', 'Cs137', 'Ba137', 'Ba138', 'La138', 'Ce138', 'La139', 'Ce139', 'Ba140', 'La140', 'Ce140', 
'Ce141', 'Pr141', 'Ce142', 'Pr142', 'Nd142', 'Ce143', 'Pr143', 'Nd143', 'Ce144', 'Nd144', 'Nd145', 'Nd146', 'Nd147', 'Pm147', 
'Sm147', 'Nd148', 'Pm148', 'Sm148', 'Pm149', 'Sm149', 'Nd150', 'Sm150', 'Pm151', 'Sm151', 'Eu151', 'Sm152', 'Eu152', 'Gd152', 
'Sm153', 'Eu153', 'Gd153', 'Sm154', 'Eu154', 'Gd154', 'Eu155', 'Gd155', 'Eu156', 'Gd156', 'Eu157', 'Gd157', 'Gd158', 'Dy158', 'Tb159', 
'Gd160', 'Tb160', 'Dy160', 'Dy161', 'Dy162', 'Dy163', 'Dy164', 'Er164', 'Ho165', 'Er166', 'Er167', 'Er168', 'Tm168', 'Tm169', 'Er170', 'Tm170']
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size
nuclides = {}
nuc_data = {} #AutoVivification()
path_to_xs = '/home/salcedop/groupr_xs/quick_prac/updating_file/smr-prac/dir/'
cols = ['material','nuclide','score','reaction-rate-MG']
depletion_rx_list = ['(n,fission)','(n,2n)','(n,3n)','(n,4n)','(n,g)','(n,p)','(n,a)']
pwd = os.getcwd()

for inuc,item in enumerate(DEPLETION_NUCLIDES):
  prep_data(item)

total_groups = 1500
comm.barrier()
openmc.capi.init(args=None,intracomm=comm)

openmc.capi.run()


time_start1 = time.time()
fuel_mats = get_fuel_mats()
res_dict = {}
dens_data = {}
mat_keys = []
start_dens = time.time()

for region,mat_id in enumerate(fuel_mats):
  if (region % size != rank):
    continue
  dens_data[mat_id] = {}
  #print("Task number {} ({}) being done by processor {} of {}".format(inuc,item,rank,size-1))
  for inuc,nucs in enumerate(openmc.capi.materials[mat_id].nuclides):
    dens_data[mat_id][nucs] = {}
    dens_data[mat_id][nucs]['dens'] = openmc.capi.materials[mat_id].densities[inuc]

end_dens = time.time()
time_to_dens = end_dens - start_dens
#print("Time to process density: {}".format(time_to_dens))
start_openmc = time.time()
openmc.capi.run()

if (rank == 0):
  flux_tally = openmc.capi.tallies[2].mean
  comm.Send([flux_tally,MPI.FLOAT],dest=1,tag=13)
elif (rank==1):
  flux_tally = np.empty(1499*len(fuel_mats),dtype='float')
  comm.Recv([flux_tally,MPI.FLOAT],source=0,tag=13)

#flux_tally = openmc.capi.tallies[2].mean

#comm.barrier()  
openmc.capi.finalize()
for region,item in enumerate(fuel_mats):
  if (region % size != rank):
    continue
  #print("Task number {} ({}) being done by processor {} of {}".format(region,item,rank,size-1))
  res_dict[item] = {}
  shift = region * 1499
  collapse(item)

print(len(res_dict))
#if (rank==0):
end_openmc = time.time()
time_to_openmc = end_openmc-start_openmc
print("Time to run OpenMC + collapse: {}".format(time_to_openmc))
#print(res_dict)
#print(res_dict)
