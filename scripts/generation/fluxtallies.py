#!/usr/bin/env python3

from pathlib import Path
import urllib.request
import re
import os
import openmc
#import openmc.data.njoy
import openmc.data.njoy_groupr
import numpy as np
import h5py
import copy
import pandas as pd
import glob
from data_nuc import DEPLETION_NUCLIDES,depletion_rx_list,cols,error_cols
from xml.etree import ElementTree as ET
try:

  from openmc.clean_xml import clean_xml_indentation
  old = True
except:
  from openmc._xml import clean_indentation
  old = False

import matplotlib
matplotlib.use('TkAgg')
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

class FluxTallies:

  def __init__(self,path_to_lib=None,total_groups=500,temperatures=None,background_xs=None,lower_lim=1.e-5,epth_lim=1.,res_lim=950.e+3,upper_lim=2.e+7,nuclide_list=None,score_list=None,tallies_dict=None):
  #check if user specified any directory. If this isn't
  #case, create one in current directory
    if path_to_lib is None:
      self.path_to_lib = Path('MG-lib')
    else:
      self.path_to_lib = Path(path_to_lib)
    if nuclide_list is None:
      nuclide_list = ['U238']
    if score_list is None:
      score_list = ['(n,gamma)']
    
    self.total_groups = total_groups
    self.lower_lim = lower_lim
    self.epth_lim = epth_lim
    self.res_lim = res_lim
    self.upper_lim = upper_lim
    self.temperatures = temperatures
    self.background_xs = background_xs
    self.nuclide_list = nuclide_list
    self.score_list = score_list
    self.tallies_dict = tallies_dict

  @property
  def total_groups(self):
	  return self._total_groups

  @property
  def lower_lim(self):
	  return self._lower_lim

  @property
  def epth_lim(self):
	  return self._epth_lim

  @property
  def res_lim(self):
	  return self._res_lim

  @property
  def upper_lim(self):
	  return self._upper_lim

  @property
  def temperatures(self):
	  return self._temperatures

  @property
  def background_xs(self):
	  return self._background_xs

  @property
  def path_to_lib(self):
	  return self._path_to_lib
  @property
  def nuclide_list(self):
    return self._nuclide_list

  @property
  def score_list(self):
    return self._score_list
  @property
  def tallies_dict(self):
    return self._tallies_dict

  @total_groups.setter
  def total_groups(self,total_groups):
	  self._total_groups = total_groups

  @lower_lim.setter
  def lower_lim(self,lower_lim):
	  self._lower_lim = lower_lim

  @upper_lim.setter
  def upper_lim(self,upper_lim):
    self._upper_lim = upper_lim

  @epth_lim.setter
  def epth_lim(self,epth_lim):
    self._epth_lim = epth_lim

  @res_lim.setter
  def res_lim(self,res_lim):
	  self._res_lim = res_lim

  @temperatures.setter
  def temperatures(self,temperatures):
	  self._temperatures = temperatures

  @background_xs.setter
  def background_xs(self,background_xs):
    self._background_xs = background_xs
  
  @nuclide_list.setter
  def nuclide_list(self,nuclide_list):
    self._nuclide_list = nuclide_list 

  @score_list.setter
  def score_list(self,score_list):
    self._score_list = score_list
  
  @tallies_dict.setter
  def tallies_dict(self,tallies_dict):
    self._tallies_dict = tallies_dict
 
  @path_to_lib.setter
  def path_to_lib(self,path_to_lib):
    self._path_to_lib = Path(path_to_lib)
    if not os.path.isdir(self._path_to_lib):
  	  self._path_to_lib.mkdir()
  
  def download_endf(self):

    for nuc in DEPLETION_NUCLIDES:
      directory = str(self.path_to_lib)+'/'+nuc
      if not os.path.isdir(directory):
        os.mkdir(directory)
      os.chdir(directory) 
      #establishing 'regular-expression' (RE) to split isotope string.
      #for example 'U238' would be split into 'U' and '238'. This is necessary
      #to write the url we have to use to download the ENDFB-VII files.
      pattern= re.compile("([a-zA-Z]+)([0-9]+)")
      components = pattern.match(nuc)
      isotope_name = components.group(1)
      isotope_A = components.group(2)
      url = 'https://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/{}/{}'.format(isotope_name,isotope_A)  
      #calling it 'tape20' to make automation across all nuclides easier.
      urllib.request.urlretrieve(url,'tape20')
      os.chdir('..')

  def SHEM_struc(self):

    SHEM_lib = h5py.File('SHEM361.hdf5','r') 
    energy_bounds = SHEM_lib["energy_bounds"].value.flatten()

    return energy_bounds

  def gen_energy_struc(self,bps_matrix,SHEM361_only=False):
    
   len_matrix = len(bps_matrix)
   shem_lib = self.SHEM_struc()
   if (SHEM361_only == False):
    for i in range(1,len_matrix):
      if (i == 1):
        group_struc_accumulate = np.logspace(np.log10(bps_matrix[i-1,0]),np.log10(bps_matrix[i,0]),int(bps_matrix[i,1]))
      
      else:    
        #cutoff_shift = 1.e-1
        #print("group_accumulate: "+str(len(group_struc_accumulate)))
        current_struct = np.logspace(np.log10(bps_matrix[i-1,0]),np.log10(bps_matrix[i,0]),int(bps_matrix[i,1]))
       # print("current_struct: "+str(len(current_struct)))
        group_struc_accumulate = np.concatenate((group_struc_accumulate,current_struct[1:]),axis=0)
    #print(len(group_struc_accumulate))
                                                                                 
    #group_struc_accumulate = np.concatenate((group_struc_accumulate,shem_lib[8:362]),axis=0)
    group_struc_accumulate = sorted(group_struc_accumulate) 
    print(len(group_struc_accumulate))
   else:
    group_struc_accumulate = sorted(shem_lib)
 
   '''
   SHEM_lib = h5py.File('SHEM361.h5','r') 
   energy_bounds = SHEM_lib["energy_bounds"].value
   group_struc_accumulate = np.concatenate((group_struc_accumulate,energy_bounds[0,:]),axis=0)
   group_struc_accumulate = sorted(group_struc_accumulate)
   print(group_struc_accumulate)
   '''
   '''
   #resolve_groups = int(self.total_groups * 0.4) 
   thermal_fast_groups = 300#int(self.total_groups * 0.1)
   resolve_groups = 500
   fast_groups = 1200#int(self.total_groups * 0.5)
   #int((self.total_groups - resolve_groups)/2)
   #Need to specify a small shift to separate the end-point of the thermal_struc
   #and the beginning of the resolve_struc. Otherwise, NJOY will complain that the structure
   #is 'out-of order' because it doesn't increase monotonically. 
   #the same issue is addressed for the resolve_struc and the fast_struc.
   cutoff_shift = 1.e-1 #eV
   thermal_struc = np.logspace(np.log10(self.lower_lim),np.log10(self.epth_lim),thermal_fast_groups)
   resolve_struc = np.logspace(np.log10(self.epth_lim+cutoff_shift),np.log10(self.res_lim),resolve_groups)
   fast_struc = np.logspace(np.log10(self.res_lim+cutoff_shift),np.log10(self.upper_lim),fast_groups)
  
   group_struc = np.concatenate((thermal_struc,resolve_struc),axis=0)
   group_struc = np.concatenate((group_struc,fast_struc),axis=0)
   '''
   return (group_struc_accumulate)

  def plot_xs(self,bps_matrix,SHEM361_only=False):

    #results_file = open('minimum_threshold.txt','w')
    mydir = os.getcwd()
    rx_list = {}
    mt_list = ['(n,p)','(n,a)']
    for imt in mt_list:
      rx_list[imt] = list()
    txt = 'nuclide,rxn,minimum_energy\n'
    DEPLETION_NUCLIDES = ["As74"]
    for irx in mt_list:
          plt.figure()
          labels = []
          for nuc in DEPLETION_NUCLIDES:
            dire = mydir + '/' + nuc
            #dire = "/home/salcedop/library-nuclide-bins-2500-groups-422-nucs-fission/"+nuc#'/home/salcedop/library-nuclide-2500-
            os.chdir(dire)
            h5file = h5py.File('{}.h5'.format(nuc),'r')
            group = list(h5file.values())[0]
            group_nuc = group['reactions/'+irx]
            arr = group_nuc['groupr'].value
            len_arr = len(arr[:,1])
            g_struc = self.gen_energy_struc(bps_matrix,SHEM361_only)
            energies = g_struc[self.total_groups-len_arr:]
            plt.loglog(energies,arr[:,1])
            labels.append('{}-{}'.format(nuc,irx))
            #except:
            # continue
          plt.xlabel('Energy (eV)')
          plt.ylabel('Cross sections (b)')
          #plt.title('107 and levels for Gd154')
          plt.title('{}, {}-MG versus energy'.format(nuc,str(irx)))
          plt.legend(labels,ncol=9,loc='lower left',prop={'size':5})
          plt.show()


  def generate_tallies_xml(self,filename_string,bps_matrix,SHEM361_only=False):

    tree = ET.parse(filename_string)
    root = tree.getroot()
    mat_filter = root[0][0].text
    energy_filter = self.gen_energy_struc(bps_matrix,SHEM361_only)

    new_data = ET.Element('tallies')

    attrib_mat = {'id':'1', 'type':'material'}
    attrib_energy = {'id':'2', 'type':'energy'}
    empty_attrib = {}

    filter_mat = ET.SubElement(new_data,'filter',attrib_mat)
    filter_ene = ET.SubElement(new_data,'filter',attrib_energy)
    ET.SubElement(filter_mat,'bins',empty_attrib)
    ET.SubElement(filter_ene,'bins',empty_attrib)
    new_data[0][0].text = mat_filter
    new_data[1][0].text = ' '.join(str(i) for i in energy_filter)
    for key,value in enumerate(self.tallies_dict.items()):
      tally_id = value[1]['id']
      tally_name = value[1]['name']
      attrib_tallies = {'id':tally_id,'name':tally_name}
      tally_level = ET.SubElement(new_data,'tally',attrib_tallies)

      ET.SubElement(tally_level,'filters',empty_attrib)
      if (value[1]['nuclides'] != ''):
         ET.SubElement(tally_level,'nuclides',empty_attrib)
      ET.SubElement(tally_level,'scores',empty_attrib)

      last_index = 1
      filter_list = value[1]['filter']
      score_list = value[1]['scores']
      new_data[key+2][0].text = ' '.join(i for i in filter_list)
      
      if (value[1]['nuclides'] != ''):
        nuclide_list = value[1]['nuclides']
        new_data[key+2][1].text = ' '.join(nuc for nuc in nuclide_list)
        last_index = 2
      new_data[key+2][last_index].text = ' '.join(i for i in score_list)

    if (old is True):
      clean_xml_indentation(new_data)
    else:
      clean_indentation(new_data)

    tree = ET.ElementTree(new_data)
    os.rename(filename_string,'tallies-original.xml')
    tree.write('./tallies.xml', xml_declaration=True,
                   encoding='utf-8', method="xml")
                                                     
  def run_njoy_arbitrary_struc(self,bps_matrix,SHEM361_only=False):

    filename = './tape20'
    temperatures = self.temperatures
    background_xs = self.background_xs
    ace='./ace'
    xsdir='./xs_info'
    moder=True
    unresr=True
    stdout=True
    save_files=True
    groupr=True
    purr=False
    acer=False
 
    one_over_e=False
    

    group_struc = self.gen_energy_struc(bps_matrix,SHEM361_only)
    ''' 
    plt.figure()
    plt.plot(group_struc)
    plt.show()
    quit()
    '''
    #In this implementation we are selecting to enter a group-structure to NJOY
    #the structure is organize into 5 columns. 'energy_update' will contain that text.
    #'format_counter' helps organize the group-structure into 5 columns when 
    #writing the njoy input. This improves the aesthetics of the NJOY input and prevents 
	  #NJOY from crashing due to any problem when reading the structure. 
    #As soon as it gets to five, we reset 'format_counter'
    #to 0 and start counting again. 'total_counter' keeps going until we hit
    #the last element in the loop, which is important because 
    #The last line needs to be treated differently (no space character after the upper limit,
    #which is 20 MeV in the default case).

    energy_text_input=''
    format_counter = 0
    total_counter = 0
 
    for energy in group_struc:
      format_counter += 1
      total_counter += 1
      if (format_counter == 5):
        if (total_counter == self.total_groups):
          extra = ''
        else:
          extra = ' ' + '\n'
          format_counter = 0
      else:
        extra = ' '
      energy_text_input += str(energy) + extra
    for nuc in DEPLETION_NUCLIDES:
      directory = str(self.path_to_lib)+'/'+nuc
      os.chdir(directory)

      openmc.data.njoy_groupr.make_xs(filename,temperatures=self.temperatures,background_xs=self.background_xs,
                                    group_struct=energy_text_input,num_groups=self.total_groups-1,ace=ace,xsdir=xsdir,moder=moder,
                                    unresr=unresr,purr=purr,groupr=groupr,acer=acer,stdout=stdout,one_over_e=one_over_e)
   
      os.chdir('..')

  def parse_groupr(self): 
    key2='for'
    key3 = 'dilution\n \n'
    key4 = '\n \n'
    depletion_rx_list = ['(n,p)','(n,a)','(n,g)','(n,fission)','(n,2n)','(n,3n)','(n,4n)']
    
    for nuc in DEPLETION_NUCLIDES:
      directory = str(self.path_to_lib) + '/'+ nuc
      os.chdir(directory)
      njoy_file = open('njoy_output_groupr_{}'.format(self.total_groups),'r').read()
      h5name = '{}_{}_groupr.hdf5'.format(nuc,self.total_groups)
      if os.path.isfile(h5name):
        os.remove(h5name)
      f = h5py.File(h5name,'w')
      rxs_group = f.create_group(nuc+'/reactions')
      for rxn in depletion_rx_list:
        key1=rxn
        try:
          #I guess you can do the parsing in one full/long
          #line but this is probably more readable.
          parse1 = njoy_file.split(key1)[-1].split(key2)[0]
          parse2 = parse1.split(key3)[-1].split(key4)[0]
          parse3 = parse2.split()
          #since we are parsing threshold-rxns, we need to know
          #the threshold index and the array with the actual data.
          #'parse3' is a 1-D array, however it holds information that can be readily
          #mapped to a 2-D array. It contains both the indexes and cross-sections (xs)
          #adjacent to each other. For example, first element is an index, second one will 
          #be the xs for that index, third one will be an index, fourth one will be 
          #its xs value and so on. So really, all even elements will be indexes, while all 
          #odd elements will correspond to xs values (remember: pyhton slicing starts at 0). 
          start_point = int(parse3[0])
          len_xs = int((len(parse3)/2))
          dummy_array = np.zeros((self.total_groups-1))
          for j in range(len_xs):
            index = 2*j + 1
            index_even = int(parse3[index - 1])
            #print(index_even)
            #one of the issues with parsing the NJOY output file is that
            #the xs values have the following format '1.5294-9' instead of '1.5294e-9'
            #therefore, before saving the xs values in the *h5 file we need to write it correctly. 
            try:
              #print(start_point)
              xs_val=parse3[index].split('+')
              #print(parse3[index_even])
              #dummy_array[j,0] = float(parse3[index_even])
              dummy_array[index_even-1] = float(xs_val[0]+'e+'+xs_val[-1])
            except:
              xs_val=parse3[index].split('-')
              #dummy_array[j,0] = float(parse3[index_even])
              dummy_array[index_even-1] = float(xs_val[0]+'e-'+xs_val[-1])
 
        except:
          start_point = 0
          dummy_array = np.zeros((1))
        
        #matching OpenMC score strings. Unfortunely instead of
        #'fission' NJOY uses '(n,fission)' and instead '(n,gamma)'
        #NJOY uses '(n,g)'
        
        #if (rxn == '(n,fission)'):
        #  rxn =  'fission'
        if (rxn == '(n,g)'):
          rxn = '(n,gamma)'
        
        irxs_group = rxs_group.create_group('{}'.format(rxn))
        irxs_group.create_dataset('groupr',data=dummy_array)
        irxs_group.attrs['start_point'] = start_point
    f.close()

  def gen_materials(self,geometry,export_xml=True):

    materials = openmc.Materials(geometry.get_all_materials().values())
    if export_xml:
      materials.export_to_xml('material.xml')
      
  def gen_tallies(self,geometry,export_xml=True):

    # Extract all fuel materials and define tally-filters.
    materials = geometry.get_materials_by_name(name='Fuel', matching=False)
    nuclides = ['U238'] #user should specify this.
    mat_filter = [openmc.MaterialFilter(materials)]
    #keep in mind that the energy filter will use the default parameters
    #to create the energy filter if nothing was specify when the FluxTallies object
    #was instantiated. 
    energy_filter = self.gen_energy_struc()

    #reaction-rate-tallies
    tally_rr = openmc.Tally(name="depletion reaction-rate tally")
    tally_rr.scores = ['(n,gamma)']
    tally_rr.nuclides = nuclides
    tally_rr.filters =  mat_filter 

    #flux-tallies
    tally_flux = openmc.Tally(name="depletion flux-tally")
    tally_flux.scores=['flux']
    fine_energy_group = openmc.EnergyFilter(energy_filter)
    tally_flux.filters = mat_filter
    tally_flux.filters.append(fine_energy_group)
    tallies = openmc.Tallies([tally_rr,tally_flux])
    if export_xml:
      tallies.export_to_xml('tallies.xml')

  def collapse_from_statepoint(self,statepoint_filename):
    depletion_rx_list = ['(n,fission)','(n,2n)','(n,3n)','(n,4n)','(n,g)','(n,p)','(n,a)']
    summary_filename = 'summary.h5'
    su = openmc.Summary(summary_filename)
    geom = su.geometry
    mat_fuel = geom.get_materials_by_name(name='Fuel', matching=False)
    #mat = geom.get_all_materials()
    #mat_fuel = mat[10015]
    sp = openmc.StatePoint(statepoint_filename)
    tallies = sp.tallies
    tally_name = tallies[2].name
    #make this more general?
    flux = sp.get_tally(name=tally_name).get_pandas_dataframe()
    df = pd.DataFrame(columns=cols)
    num_groups = self.total_groups - 1
    regions = int(len(flux) / num_groups)
    running_sum = np.zeros([7])
    counter = -1 
    
    for region in range(regions): 
      shift = region * num_groups
      mat_id = flux['material'][shift]     
      nuclides_dens = mat_fuel[region].get_nuclide_densities()
      for nuc in DEPLETION_NUCLIDES:

        dens = nuclides_dens[nuc][1]
        directory = str(self.path_to_lib) + '/' + nuc
        os.chdir(directory)
        h5file = glob.glob('{}_{}_groupr.hdf5'.format(nuc,self.total_groups))
        h5file = h5py.File(h5file[0],'r')
        group = list(h5file.values())[0]
        for irx,rx in enumerate(depletion_rx_list):
          counter += 1
          running_sum[irx] = 0.
          index = 'reactions/'+rx
          igroup = group[index]
          threshold = igroup.attrs['start_point']
          end = len(igroup['groupr'])
          for ielement in range(end):
            if (threshold == 0):
              group_xs = 0.
              ielement = 1
            else:
              group_xs = igroup['groupr'][ielement]
              running_sum[irx] += flux['mean'][ielement+threshold+shift-1]*group_xs
          
          if (rx == '(n,fission)'):
            rx =  'fission'
          if (rx == '(n,g)'):
            rx = '(n,gamma)'
          
          df.loc[counter,'reaction-rate-MG'] = running_sum[irx] * dens
          df.loc[counter,'material'] = mat_id
          df.loc[counter,'nuclide'] = nuc
          df.loc[counter,'score'] = rx
      os.chdir('..')
    df.to_csv('reaction-rates-MG.csv',sep=' ') 

  def check_accuracy_from_summary(self,flux_tally_filename=None,reaction_rate_filename=None):
    summary_filename = 'summary.h5'
    if flux_tally_filename is None:
      flux_tally_filename = 'reaction-rates-MG.csv'
    if reaction_rate_filename is None:
      reaction_rate_filename = 'reaction-rates.csv'
    summary_file = '/home/salcedop/RESULTS/profiling/Throughput/block17_testing_direct-tallies/simulations/smr-reaction-rate-dens-included/summary.h5'    
    error_df = pd.DataFrame(columns=error_cols)
    rrates_df = pd.read_table(reaction_rate_filename,sep=' ')  
    fluxtally_df = pd.read_table(flux_tally_filename,sep=' ')
    #total absorption across entire system.
    total_abs = rrates_df['mean'].sum(axis=0) 
    #need to know number-density of each nuclide in each
    #material to compute macroscopic rate.
    su = openmc.Summary(summary_filename)
    geom = su.geometry
    cnt = -1
    for imat,mat in enumerate(rrates_df['material']):  
      cnt += 1
      nuc = rrates_df.loc[imat,'nuclide']
      
      dummy_ft = fluxtally_df[fluxtally_df['nuclide']==nuc]
      dummy_rr = rrates_df[rrates_df['nuclide']==nuc]
      #sort 'nuclide' and 'score' columns to make sure both dataframes
      #are in the same order.
      dummy_ft = dummy_ft.sort_values(by=['nuclide']).sort_values(by=['score'])
      dummy_rr = dummy_rr.sort_values(by=['nuclide']).sort_values(by=['score'])
      #new_index_rr =  dummy_rr[dummy_rr['score']==irxn].index[0]
      #new_index_ft = dummy_ft[dummy_ft['score']==irxn].index[0]
      if cnt > 6:
        cnt = 0
      new_index_rr =  dummy_rr['score'].index[cnt]
      new_index_ft = dummy_ft['score'].index[cnt]

      irow_ft = dummy_ft.loc[new_index_ft,'reaction-rate-MG']
      irow_rr = dummy_rr.loc[new_index_rr,'mean']
      error_df.loc[imat,'reaction-rate-MG'] = irow_ft
      error_df.loc[imat,'reaction-rate-MC'] = irow_rr
      error_df.loc[imat,'Error[pcm]'] = abs(irow_ft - irow_rr) / (total_abs * 1.E-5)
      
      error_df.loc[imat,'material'] = mat
      error_df.loc[imat,'nuclide'] = nuc
      error_df.loc[imat,'score'] = dummy_rr.loc[new_index_rr,'score']
      
    error_df.to_csv('reaction-rates-error.csv',sep=' ')
