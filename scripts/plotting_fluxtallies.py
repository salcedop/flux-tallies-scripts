#!/usr/bin/env python3

from pathlib import Path
import urllib.request
import re
import os
import openmc
import openmc.data.njoy
#import openmc.data.njoy_groupr
import numpy as np
import h5py
import copy
import pandas as pd
import glob
from data_nuc import DEPLETION_NUCLIDES,depletion_rx_list,cols,error_cols,mt_list,source

import matplotlib
matplotlib.use('Agg')
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

class FluxTallies:


  def __init__(self,path_to_lib=None,total_groups=500,temperatures=None,background_xs=None,lower_lim=1.e-5,epth_lim=1.,res_lim=950.e+3,upper_lim=2.e+7):
  #check if user specified any directory. If this isn't
  #case, create one in current directory
    if path_to_lib is None:
      self.path_to_lib = Path('MG-lib')
    else:
      self.path_to_lib = Path(path_to_lib)
    
    self.total_groups = total_groups
    self.lower_lim = lower_lim
    self.epth_lim = epth_lim
    self.res_lim = res_lim
    self.upper_lim = upper_lim
    self.temperatures = temperatures
    self.background_xs = background_xs

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
  
  @path_to_lib.setter
  def path_to_lib(self,path_to_lib):
    self._path_to_lib = Path(path_to_lib)
    if not os.path.isdir(self._path_to_lib):
  	  self._path_to_lib.mkdir()
  def p(self):
    '''
    plt.figure()
    plt.plot([1,2,3,4])
    plt.ylabel('some numbers')
    plt.savefig('p3.png')
    #plt.show()
    '''
    os.chdir(source)
    plt.figure()
    labels = []
    filename = 'B_010_293.6K.ace'
    inuc_endf = openmc.data.IncidentNeutron.from_ace(filename)
    #try:
    energy_values = inuc_endf[103].xs['294K'].x
    threshold = energy_values[0]
    xs_values = inuc_endf[103].xs['294K'].y
    plt.plot(energy_values,xs_values)
    ituple = (threshold,'B10')
    #rx_list[103].append(ituple)
    #labels.append('B10-103')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Cross sections (b)')
    #plt.legend(labels,ncol=4,loc='lower center')
    #plt.savefig('p3.png')
    plt.savefig('BBB.png')
    #plt.show()
  def plot_xs(self):

    #results_file = open('minimum_threshold.txt','w')
    mydir = os.getcwd()
    os.chdir(source)
    rx_list = {}
    #mt_list = [107]
    for imt in mt_list:
      rx_list[imt] = list()
    txt = 'nuclide,rxn,minimum_energy\n'
    #DEPLETION_NUCLIDES = ["O16"]
    for irx in mt_list:
          plt.figure()
          labels = []
          for nuc in DEPLETION_NUCLIDES:
            temp_len = len(nuc.split("_"))
            if (temp_len == 2):
              pattern = re.compile("([a-zA-Z]+)([0-9]+)([_]+)([a-z0-9]+)")
              components = pattern.match(nuc)
              isotope_name = components.group(1)
              length_A = len(components.group(2))
              isotope_A = components.group(2) + components.group(4)
            elif (temp_len == 1):
              pattern= re.compile("([a-zA-Z]+)([0-9]+)")
              components = pattern.match(nuc)
              length_A = len(components.group(2))
              isotope_name = components.group(1)
              isotope_A = components.group(2)
        
            else:
              raise ValueError("{} not in valid string format".format(nuc))
            if (length_A < 3) and (length_A > 1):
              isotope_A = '0' + isotope_A
            elif(length_A < 2):
              isotope_A = '00' + isotope_A
            filename = '{}_{}_293.6K.ace'.format(isotope_name,isotope_A)
            inuc_endf = openmc.data.IncidentNeutron.from_ace(filename)
            try:
             energy_values = inuc_endf[irx].xs['294K'].x
             threshold = energy_values[0]
             xs_values = inuc_endf[irx].xs['294K'].y
             plt.loglog(energy_values,xs_values)
             ituple = (threshold,nuc)
             rx_list[irx].append(ituple)
             labels.append('{}-{}'.format(nuc,irx))
            except:
             continue
          plt.rcParams['axes.labelsize'] = 16

          plt.rcParams['axes.xmargin'] = 0

          plt.rcParams['axes.ymargin'] = 0

          plt.rcParams['font.family'] = 'serif'

          plt.rcParams['font.serif'] = ['Times']

          plt.rcParams['font.size'] = 10

          plt.rcParams['mathtext.fontset'] = 'stix'

          plt.rcParams['pdf.use14corefonts'] = True

          plt.rcParams['savefig.bbox'] = 'tight'
          plt.xlabel('Energy (eV)')
          plt.ylabel('Cross section (b)')
          plt.title('{} (MT={}) Cross-section versus Energy'.format(nuc,str(irx)))
          plt.savefig('{}-{}-updated.pdf'.format(nuc,irx))
          #plt.savefig('{}/{}.png'.format(mydir,irx))
          #plt.show()
            #except:
              #continue
    '''
    os.chdir(mydir)
    for imt in mt_list:
     simt = str(imt)
     f = open('{}.txt'.format(simt),'w')
     txt = '{}\n----------------------\n'.format(simt)
     try:
      s = sorted(rx_list[imt])
      for itule in s:
        print(itule[1])
        print(itule[0])
        txt += itule[1]+','+str(itule[0])+'\n'
      #print(s)
      #itxt = s[0][1]+','+str(imt)+','+str(s[0][0])+'\n'
      f.write(txt)
      f.close()
     except:
      continue
     '''
    #results_file.write(txt)
    #results_file.close()
       
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

  def gen_energy_struc(self):
 
  
    resolve_groups = int(self.total_groups * 0.7)
    thermal_fast_groups = int((self.total_groups - resolve_groups)/2)
    #Need to specify a small shift to separate the end-point of the thermal_struc
    #and the beginning of the resolve_struc. Otherwise, NJOY will complain that the structure
    #is 'out-of order' because it doesn't increase monotonically. 
    #the same issue is addressed for the resolve_struc and the fast_struc.
    cutoff_shift = 1.e-1 #eV
    thermal_struc = np.logspace(np.log10(self.lower_lim),np.log10(self.epth_lim),thermal_fast_groups)
    resolve_struc = np.logspace(np.log10(self.epth_lim+cutoff_shift),np.log10(self.res_lim),resolve_groups)
    fast_struc = np.logspace(np.log10(self.res_lim+cutoff_shift),np.log10(self.upper_lim),thermal_fast_groups)
  
    group_struc = np.concatenate((thermal_struc,resolve_struc),axis=0)
    group_struc = np.concatenate((group_struc,fast_struc),axis=0)
 
    return group_struc

  def run_njoy_arbitrary_struc(self):

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
    

    group_struc = self.gen_energy_struc()
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

  def threshold(self):


    results_file = open('minimum_threshold.txt','w')
    rx_list = {}
    os.chdir(source)
    for imt in mt_list:
      rx_list[imt] = list()
    txt = 'nuclide,rxn,minimum_energy\n'
    for nuc in DEPLETION_NUCLIDES:
      filename = '{}.h5'.format(nuc)
      file15 = h5py.File(filename,'r')
      group_name = '{}/reactions'.format(nuc)
      main_group = file15[group_name]
      for irx in mt_list:
        try:
          imt_name = 'reaction_{}/294K/xs'.format(irx)
          threshold = main_group[imt_name][1]
        
          ituple = (threshold,nuc)
          rx_list[irx].append(ituple)
        except:
          continue

    for imt in mt_list:
     try:
      s = sorted(rx_list[imt])
      itxt = s[0][1]+','+str(imt)+','+str(s[0][0])+'\n'
      txt += itxt
     except:
      continue
    print(rx_list)
    results_file.write(txt)
    results_file.close()

  def threshold2(self):
    results_file = open('minimum_threshold.txt','w')
    rx_list = {}
    os.chdir(source)
    for imt in mt_list:
      rx_list[imt] = list()
    txt = 'nuclide,rxn,minimum_energy\n'
    for nuc in DEPLETION_NUCLIDES:
      pattern= re.compile("([a-zA-Z]+)([0-9]+)")
      components = pattern.match(nuc)
      isotope_name = components.group(1)
      isotope_A = components.group(2)
      length_A = len(isotope_A)
      if (length_A < 3) and (length_A > 1):
        isotope_A = '0' + isotope_A
      elif(length_A < 2):
         isotope_A = '00' + isotope_A
      filename = '{}_{}_293.6K.ace'.format(isotope_name,isotope_A)
      inuc_endf=openmc.data.IncidentNeutron.from_ace(filename)
      labels = []
      for irx in mt_list:
        #try:
          energy_values = inuc_endf[irx].xs['294K'].x
          threshold = energy_values[0]
          xs_values = inuc_endf[irx].xs['294K'].y
          plt.figure()
          plt.loglog(energy_values,xs_values)
          #plt.xlabel('Energy (eV)')
          #plt.ylabel('Cross sections (b)')
          ituple = (threshold,nuc)
          rx_list[irx].append(ituple)
          labels.append('{}-{}'.format(nuc,irx))
          #plt.legend(labels,ncol=4,loc='lower center')
          plt.savefig(nuc+'_'+str(irx)+'.png')
        #except:
          #continue
    for imt in mt_list:
     try:
      s = sorted(rx_list[imt])
      print(s)
      itxt = s[0][1]+','+str(imt)+','+str(s[0][0])+'\n'
      txt += itxt
     except:
      continue
    
    results_file.write(txt)
    results_file.close()


  def parse_groupr(self): 
    key2='for'
    key3 = 'dilution'
    key4 = 'group'
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
          print(start_point)
          len_xs = int((len(parse3)/2))
          dummy_array = np.zeros(len_xs)
          for j in range(len_xs):
            index = 2*j + 1
            #one of the issues with parsing the NJOY output file is that
            #the xs values have the following format '1.5294-9' instead of '1.5294e-9'
            #therefore, before saving the xs values in the *h5 file we need to write it correctly. 
            try:
              xs_val=parse3[index].split('+')
              dummy_array[j] = float(xs_val[0]+'e+'+xs_val[-1])
            except:
              xs_val=parse3[index].split('-')
              dummy_array[j] = float(xs_val[0]+'e-'+xs_val[-1])

        except:
          start_point = 0
          dummy_array = np.zeros(1)
      
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
