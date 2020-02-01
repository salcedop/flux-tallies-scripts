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
from data_nuc import O16Peaks,DEPLETION_NUCLIDES,depletion_rx_list,def_300group_struc
from xml.etree import ElementTree as ET

#DEPLETION_NUCLIDES = ['H1', 'H2','H3','Na23']
DEPLETION_NUCLIDES = ['U235']
try:

  from openmc.clean_xml import clean_xml_indentation
  old = True
except:
  from openmc._xml import clean_indentation
  old = False

#in case user decides to set another directory
original_dir = os.getcwd()
class FluxTallies:
  """This class handles the logic of fluxtallies within OpenMC..

  """
  def __init__(self,path_to_lib=None,total_groups=300,temperatures=[300.],nuclide_list=['U238'],score_list=['(n,gamma)'],tallies_dict=None):
    
    if path_to_lib is None:
      working_dir = original_dir
    else:
      working_dir = path_to_lib

    self.path_to_lib = Path(working_dir)
    #if not os.path.isdir(self.path_to_lib):
  	#  self.path_to_lib.mkdir()
    #os.chdir(working_dir)

    #self.path_to_lib = Path('MG-lib')  
    self.total_groups = total_groups
    self.temperatures = temperatures
    self.nuclide_list = nuclide_list
    self.score_list = score_list
    self.tallies_dict = tallies_dict

  @property
  def total_groups(self):
	  return self._total_groups

  @property
  def temperatures(self):
	  return self._temperatures

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

  @temperatures.setter
  def temperatures(self,temperatures):
	  self._temperatures = temperatures
 
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
    os.chdir(original_dir)
    self._path_to_lib = Path(path_to_lib) 
    if not os.path.isdir(self._path_to_lib):
  	  self._path_to_lib.mkdir()
    os.chdir(self._path_to_lib)
    
  def download_endf(self):
    """This method downloads endf files
       for all nuclides on the DEPLETION_NUCLIDES
       string.
    """
    for nuc in DEPLETION_NUCLIDES:
      print(os.getcwd())
      directory = nuc
      if not os.path.isdir(directory):
        os.mkdir(directory)
      os.chdir(directory) 
      #establishing 'regular-expression' (RE) to split isotope string.
      #for example 'U238' would be split into 'U' and '238'. This is necessary
      #to write the url we have to use to download the ENDFB-VII files.
      #check if the string splits to verify 

      temp_len = len(nuc.split("_"))
      #to identified excited state nuclides (i.e., "_m1")
      if (temp_len == 2):
        pattern = re.compile("([a-zA-Z]+)([0-9]+)([_]+)([a-z0-9]+)")
        components = pattern.match(nuc)
        isotope_name = components.group(1)
        isotope_A = components.group(2) + components.group(4)
      #to indentify regular nuclides
      elif (temp_len == 1):
        pattern= re.compile("([a-zA-Z]+)([0-9]+)")
        components = pattern.match(nuc)
        isotope_name = components.group(1)
        isotope_A = components.group(2)
        
      else:
        raise ValueError("{} not in valid string format".format(nuc))

      url = 'https://t2.lanl.gov/nis/data/data/ENDFB-VII.1-neutron/{}/{}'.format(isotope_name,isotope_A)  
      #calling it 'tape20' to make automation across all nuclides easier.
      print(url)
      urllib.request.urlretrieve(url,'tape20')
      os.chdir('..')

  def gen_bps_matrix(self):
    
    """Generate matrix of breaking points which will be used to generate
       energy structure

        Returns
        ------
     
        bps_matrix : numpy.float.matrix
                     Two-column matrix of points where the first column corresponds to the points that
                     define the groups while the second one corresponds to the number of groups of that 
                     given interval.
    """
    #built-in bps_matriz depending on the number of total matrix.

    if ((self.total_groups == 501) or (self.total_groups == None)):
        groups_400KeV_3MeV = 101.
        groups_epithermal_resolved = 51.
        groups_above_10MeV = 91. 
    
    elif (self.total_groups == 301):
        groups_400KeV_3MeV = 11.
        groups_epithermal_resolved = 11.
        groups_above_10MeV = 21. 

    else:
        raise ValueError("{} not valid for total number of groups, "
                         "please choose either 301 or 501 groups".format(self.total_groups))

    bps = np.ones(len(O16Peaks)) * 11.
    bps[0] = groups_400KeV_3MeV
    bps_matrix = np.matrix([[1.E-5,0.],[4.E+4,51.]])
    last_row = [[2.E+7,91.]]
    len_bps = len(bps_matrix)
    
    for irow,row in enumerate(O16Peaks):
        new_row = [[row,bps[irow]]]
        bps_matrix = np.concatenate((bps_matrix,new_row))
    bps_matrix = np.concatenate((bps_matrix,last_row))

    return bps_matrix

  def gen_energy_struc(self,default_group_struc=False):
    """Generate energy structure based on inputs

        Parameters
        ----------
            .
 
        default_group_struc : bool
            user can either select a built-in 300 group
            structure, which includes more groups in the high energy range
            above 1 MeV, or they can input their own structure.

        Returns
        ------
        group_struc_accumulate : energy group structure.

    """
     
    bps_matrix = self.gen_bps_matrix()
    if (default_group_struc):
        self._total_groups = 300
        group_struc_accumulate = sorted(def_300group_struc)
    else: 
        #check that total groups input matches
        #the sum of all the groups provided in the
        #bps_matrix.
        sum_bps = 0
        len_matrix = len(bps_matrix)
        for i in range(1,len_matrix):
            sum_bps += bps_matrix[i,1] - 1
            if (i == 1):
                group_struc_accumulate = np.logspace(np.log10(bps_matrix[i-1,0]),np.log10(bps_matrix[i,0]),int(bps_matrix[i,1]))
      
            else:    
                current_struct = np.logspace(np.log10(bps_matrix[i-1,0]),np.log10(bps_matrix[i,0]),int(bps_matrix[i,1]))
                group_struc_accumulate = np.concatenate((group_struc_accumulate,current_struct[1:]),axis=0)
        group_struc_accumulate = sorted(group_struc_accumulate) 
        if (sum_bps != self._total_groups):
           self._total_groups = int(sum_bps)
           print('WARNING: The input for total number of groups did not match'
                 'the length of the break points matrix.')
        
    return (group_struc_accumulate)

  def generate_tallies_xml(self,filename_string,default_group_struc=False):
    """Generate tallies.xml with hybrid tallies by using the energy group structure
       as a filter together with the material filter of an original tallies.xml.

        Parameters
        ----------
        filename_string : string
            string with the path to a tallies.xml file that contains
            the reaction rate tally of a given geometry. This method will
            expand on this by adding the extra flux tally method and by 
            reducing the number of scores in the original reaction-rate tally
            from 7 to 2.
        SHEM361_only : bool
            In case the user wants to use only the SHEM361 library.

    """
    #reads a tallies.xml file with a reaction rate tally and uses the 
    #material filter to build a new tallies.xml with hybrid tallies using
    #the energy filter.  

    tree = ET.parse(filename_string)
    root = tree.getroot()
    mat_filter = root[0][0].text
    energy_filter = self.gen_energy_struc(default_group_struc)

    new_data = ET.Element('tallies')

    #attributes for hybrid tally.
    attrib_mat = {'id':'1', 'type':'material'}
    attrib_energy = {'id':'2', 'type':'energy'}
    empty_attrib = {}

    filter_mat = ET.SubElement(new_data,'filter',attrib_mat)
    filter_ene = ET.SubElement(new_data,'filter',attrib_energy)
    ET.SubElement(filter_mat,'bins',empty_attrib)
    ET.SubElement(filter_ene,'bins',empty_attrib)
    
    new_data[0][0].text = mat_filter
    new_data[1][0].text = ' '.join(str(i) for i in energy_filter)

    #writing '<tally>' object
    for key,value in enumerate(self.tallies_dict.items()):
      tally_id = value[1]['id']
      tally_name = value[1]['name']
      attrib_tallies = {'id':tally_id,'name':tally_name}
      tally_level = ET.SubElement(new_data,'tally',attrib_tallies)

      ET.SubElement(tally_level,'filters',empty_attrib)
      
      if (value[1]['nuclides'] != ''):
          ET.SubElement(tally_level,'nuclides',empty_attrib)
      ET.SubElement(tally_level,'scores',empty_attrib)
      #last_index could be either a 1 or a 2, depending on
      #whether we are doing a reaction rate tally (i.e, last_index=1) or
      #a flux tally (i.e., last_index = 2)
      last_index = 1
      filter_list = value[1]['filter']
      score_list = value[1]['scores']
      #offset by '2' indicates we are writing the tally objects below
      #the material and energy filters
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
    #renaming original 'tallies.xml' to 'tallies-original.xml'
    #to distinguish it from  the new hybrid 'tally.xml' file.
    tree = ET.ElementTree(new_data)
    #os.rename(filename_string,'tallies-original.xml')
    tree.write('./tallies-hybrid-tally.xml', xml_declaration=True,
                   encoding='utf-8', method="xml")
                                                     
  def run_njoy_arbitrary_struc(self,default_group_struc=False):

    """Run NJOY using the input structure.

        Parameters
        ----------
        SHEM361_only : bool
            In case the user wants to use only the SHEM361 library.

    """
    filename = './tape20'
    temperatures = self.temperatures
    ace='./ace'
    xsdir='./xs_info'
    moder=True
    unresr=True
    
    groupr=True
    #these are false when using GROUPR
    purr=False
    acer=False
    #this is 
    #one_over_e : bool, optional
    #Thus far, when running GROUPR we have to provide an arbitrary group structure.
    #Howerer, NJOY also needs a weight-function to compute the group cross-sections.
    #If this variable is true then NJOY will use the following weight-function: 1/E + maxwellian +  
    # fission spectrum. By default, a flat weight-function is selected. 
    #This option works well for fine-group structure like the ones that will be needed when doing flux-tallies. 
    one_over_e=False
    #'stdout' will print the NJOY output to the screen.. 
    # Always good for debugging. 'save_files' will tell the python API
    # to save a text file with the results so that it can get parsed later on. 
    stdout=True
    save_files=True
   
    #generating group structure.
    group_struc = self.gen_energy_struc(default_group_struc)

    #In this implementation we are selecting to enter a group-structure to NJOY
    #the structure is organize into 5 columns. 'energy_update' will contain the text.
    #'format_counter' helps organize the group-structure into 5 columns when 
    #writing the njoy input. This improves the aesthetics of the NJOY input and prevents 
	  #NJOY from crashing due to any problem when reading the structure. 
    #As soon as it gets to five, we reset 'format_counter'
    #to 0 and start counting again. 'total_counter' keeps going until we hit
    #the last element in the loop, which is important because 
    #The last line needs to be treated differently (there can't be white space after the upper limit,
    #which is 20 MeV in the default case).

    energy_text_input=''
    format_counter = 0
    total_counter = 0 
    for energy in group_struc:
      format_counter += 1
      total_counter += 1
      if (format_counter == 5):
        if (total_counter == (self.total_groups+1)):
          extra = ''
        else:
          extra = ' ' + '\n'
          format_counter = 0
      else:
        extra = ' '
      energy_text_input += str(energy) + extra
    #looping over all nuclides. 
    for nuc in DEPLETION_NUCLIDES:
      directory = nuc
      os.chdir(directory)
      #when running GROUPR, you must specify the total number of groups minus 1.
      openmc.data.njoy_groupr.make_xs(filename,lorde=1,temperatures=self.temperatures,background_xs=[1000.,10.],
                                    group_struct=energy_text_input,num_groups=self.total_groups,ace=ace,xsdir=xsdir,moder=moder,
                                    unresr=unresr,purr=purr,groupr=groupr,acer=acer,stdout=stdout,one_over_e=one_over_e)
   
      os.chdir('..')

  def parse_groupr(self): 
    """Parse NJOY results.
    """
    #establishing keys to parse njoy output files.
    #key1 varies from '(n,fission)','(n,p)','(n,a)','(n,g)' and '(n,xn)'
    #where x can be 2,3 and 4.
   #300k 
    single_background_xs_keys = ['for','dilution\n \n','\n \n']
    multiple_background_xs_keys = ['for','\n \n','\n \n']
    len_temps = len(self.temperatures)
    for nuc in DEPLETION_NUCLIDES:
        directory = nuc
        os.chdir(directory)
        h5name = '{}.h5'.format(nuc)
        if os.path.isfile(h5name):
           os.remove(h5name)
        njoy_file = open('njoy_output_groupr_{}'.format(self.total_groups),'r').read()
        f = h5py.File(h5name,'w')
        for itemp,temp in enumerate(self.temperatures):
                for rxn in depletion_rx_list:
                    key1=rxn
                    #Here we are matching OpenMC score strings. 
                    #Unfortunely instead of 'fission' NJOY uses '(n,fission)' and instead '(n,gamma)'
                    #NJOY uses '(n,g)'
                    if (rxn == '(n,fission)'):
                       rxn =  'fission'
                    if (rxn == '(n,g)'):
                       rxn = '(n,gamma)'
                    #not all nuclides have all 7 depletion reactions, so we
                    #first 'try' to parse it
                    try:
                            #I guess you can do the parsing in one full/long
                            #line but this is probably more readable.
                        
                            parse = njoy_file.split(key1)[1+itemp].split('\n \n')[2].split('\n \n')[0].split()
                            #since we are parsing threshold-rxns, we need to know
                            #the threshold index and the array with the actual data.
                            #'parse3' is a 1-D array, however it holds information that can be readily
                            #mapped to a 2-D array. It contains both the indexes and cross-sections (xs)
                            #adjacent to each other. For example, first element is an index, second one will 
                            #be the xs for that index, third one will be an index, fourth one will be 
                            #its xs value and so on. So really, all even elements will be indexes, while all 
                            #odd elements will correspond to xs values. 
                            start_point = int(parse[0]) 
                            #for iback,back in enumerate(self.background_xs):
                            #index_xs =
                            #index_parse = -1 * (len_background_xs + 1)
                            group_name = '{}/reactions/{}/{}K'.format(nuc,rxn,temp)
                            
                            irxs_group = f.create_group(group_name)
                            len_xs = int((len(parse)/2))
                          
                            dummy_array = np.zeros((self.total_groups))
                            for j in range(len_xs):
                                index_xs = 2*j + 1 #len_background_xs + 1
                                index_parse =  index_xs - 1#len_background_xs + 1
                                index_dummy = int(parse[index_parse])
                                #print(index_even)
                                #one of the issues with parsing the NJOY output file is that
                                #the xs values have the following format '1.5294-9' instead of '1.5294e-9'
                                #therefore, before saving the data in *h5 format we need to write it correctly. 
                                try:
                                   #print(start_point)
                                   xs_val=parse[index_xs].split('+')
                                   
                                   #print(parse3[index_even])
                                   #dummy_array[j,0] = float(parse3[index_even])
                                   dummy_array[index_dummy-1] = float(xs_val[0]+'e+'+xs_val[-1])
                                except:
                                  # print(xs_val)
                                   xs_val=parse[index_xs].split('-')
                                   #dummy_array[j,0] = float(parse3[index_even])
                                   dummy_array[index_dummy-1] = float(xs_val[0]+'e-'+xs_val[-1])
                                 
        
                            irxs_group.create_dataset('groupr',data=dummy_array)
                            irxs_group.attrs['start_point'] = start_point
 
                    except:
                            start_point = 0
                            dummy_array = np.zeros((1))
                            
                            group_name = '{}/reactions/{}/{}K/'.format(nuc,rxn,temp)
                            print("from back: "+group_name)
                            irxs_group = f.create_group(group_name)
                            irxs_group.create_dataset('groupr',data=dummy_array)
                            irxs_group.attrs['start_point'] = start_point
        
        os.chdir('..')
        f.close()

  def generating_cross_section_xml(self):
  
     """Generate 'cross_section.xml'
     """
     library = openmc.data.DataLibrary()

     for nuc in DEPLETION_NUCLIDES:
        inuc = "{}/{}.h5".format(nuc,nuc)
        library.register_file(inuc)
     library.export_to_xml()
