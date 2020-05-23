#!/usr/bin/env python3
import os
from data_nuc import DEPLETION_NUCLIDES,def_300group_struc
import openmc
import openmc.data
import openmc.data.njoy_groupr

new_pwd = os.getcwd()

DEPLETION_NUCLIDES = ['H1', 'H2','Na23']

filename = './tape20'
temperatures = [300.,1000.]
ace='./ace'
xsdir='./xs_info'
moder=True
unresr=True  
groupr=True
#these are false when using GROUPR
purr=False
acer=False
#one_over_e : bool, optional
#Thus far, when running GROUPR we have to provide an arbitrary group structure.
#Howerer, NJOY also needs a weight-function to compute the group cross-sections.
#If this variable is true then NJOY will use the following weight-function: 1/E + maxwellian +  
# fission spectrum. By default, a flat weight-function is selected. 
#which works well for fine-group structures like the ones that will be used when doing flux-tallies. 
one_over_e=True
#'stdout' will print the NJOY output to the screen.. 
# Always good for debugging. 'save_files' will tell the python API
# to save a text file with the results so that it can get parsed later on. 
stdout=True
save_files=True   
#generating group structure.
group_struc = def_300group_struc
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
 
for nuc in DEPLETION_NUCLIDES:
   directory = nuc
   os.chdir(directory)
   #when running GROUPR, you must specify the total number of groups minus 1.
   openmc.data.njoy_groupr.make_xs(filename,lorde=1,temperatures=temperatures,
                                    group_struct=None,num_groups=300,ace=ace,xsdir=xsdir,moder=moder,
                                    unresr=unresr,purr=purr,groupr=groupr,acer=acer,stdout=stdout,one_over_e=one_over_e)
   
   os.chdir('..')
