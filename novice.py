#!/usr/bin/env python
# coding: utf-8

# In[24]:


#import os
#import sys
#curr_dir=os.getcwd()
#sys.path.append(curr_dir)
import mbuild as mb # all lengths in mbuild are in nm
import novice_functions
import signac
import flow
import os
from shutil import copyfile
from subprocess import call


# In[3]:


class Cl2(mb.Compound): # this class builds a chlorine molecule with a bond-length given in the chlorine2 x coor (nm)
    def __init__(self):
        super(Cl2, self).__init__()
        
        chlorine1= mb.Particle(pos=[0.0, 0.0, 0.0], name='Cl')
        chlorine2= mb.Particle(pos=[0.2, 0.0, 0.0], name='Cl')
        self.add([chlorine2,chlorine1])
        self.add_bond((chlorine2,chlorine1))


# In[13]:


cl2=Cl2();
molecule=cl2;
molecule.save('molecule.pdb',overwrite='True')
# you have defined the molecule, now please define the forcefield (XC Functional), Temperature (in Kelvin) and Length (nm) of the box

functional='BLYP'
dir= '/home/siepmann/singh891/cp2k-6.1.0/data/' # your cp2k data folder with a slash('/') at the end
project_name='chlorine' # change this to whatever you want to, must be a string
temperature= 177 # units are in Kelvin
box_length=1.1 # box side length in nm
number_of_molecules=32 # number of molecules, int


# In[14]:


novice_functions.optimize_molecule(molecule,functional, box_length,dir) # now we have the file for geo_opt

# this function writes the cp2k input file and the structure of the molecule in xyz file for cp2k
# Now you have a file mol_opt.inp and mol_unopt_coord.xyz in the folder that you need to submit on the cluster


# In[17]:


simulation = signac.init_project(project_name)
sp={'Temperature':temperature,'N':number_of_molecules,'L':box_length}
simulation.open_job(sp).init()


# In[28]:


i=0;
for job in simulation:
    if i==0:
        base=os.getcwd();
        working_dir=job.workspace();
        copyfile('{}/mol_opt.inp'.format(base),'{}/mol_opt.inp'.format(working_dir))
        copyfile('{}/mol_unopt_coord.xyz'.format(base),'{}/mol_unopt_coord.xyz'.format(working_dir))
    i+=1;


# In[ ]:




