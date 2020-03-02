#!/usr/bin/env python
# coding: utf-8

# In[1]:





import mbuild as mb # all lengths in mbuild are in nm
import novice_functions
import signac
import flow
import os
from shutil import copyfile
from subprocess import call
import numpy as np
from subprocess import Popen, PIPE



####User input, ideally should read from novice.py 
functional='BLYP'
dir= '/home/siepmann/singh891/cp2k-6.1.0/data/' # your cp2k data folder with a slash('/') at the end
project_name='chlorine' # change this to whatever you want to, must be a string
temperature= 177.0 # units are in Kelvin
box_length=1.1 # box side length in nm
number_of_molecules=20 # number of molecules, int
num_atoms_one_molecule=2;
simulation_time=1 # time in ps for which you want to run the simulation
############input ends here


################
string="tail -{} molecule_opt-pos-1.xyz > opt_coor.xyz".format(num_atoms_one_molecule)
call(string,shell=True)
table=0*np.empty([1, 3])#dummy array to start with;
scaling=0*np.empty([1,3])# scaling for nm to A
with open('opt_coor.xyz') as input_data:
    for line in input_data:
        n, x, y,z = line.strip().split()
        table=np.concatenate((table, np.array([[float(x),float(y),float(z)]])), axis=0);
        #scaling=np.concatenate((scaling, np.array([[0.1,0.1,0.1]])), axis=0);
table=np.delete(table,0,0)
#scaling=np.delete(table,0,0)
table=table*0.1;
molecule=mb.load('molecule.pdb');

#write the input for md1 and md2


# In[ ]:





# In[2]:


novice_functions.run_md_pre(molecule,functional,temperature, box_length,dir,number_of_molecules,simulation_time,project_name,table)


# In[3]:





# In[4]:


novice_functions.run_md_main(molecule,functional,temperature, box_length,dir,number_of_molecules,simulation_time,project_name)


# In[16]:





# In[17]:





# In[ ]:





# In[ ]:


simulation = signac.get_project()


# In[ ]:


for job in simulation:
    base=os.getcwd();
    working_dir=job.workspace();
    copyfile('{}/md-pre.inp'.format(base),'{}/md-pre.inp'.format(working_dir))
    
    copyfile('{}/md-main.inp'.format(base),'{}/md-main.inp'.format(working_dir))
    copyfile('{base}/'.format(base=base)+project_name+'.xyz','{}/'.format(working_dir)+project_name+'.xyz')


# In[ ]:





# In[ ]:




