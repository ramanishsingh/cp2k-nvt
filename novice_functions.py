#!/usr/bin/env python
# coding: utf-8

# In[1]:


import mbuild as mb
import subprocess
from cssi_cp2k.classes import SIM as sim
# read this: https://www.cp2k.org/_media/events:2015_cecam_tutorial:watkins_optimization.pdf for more knowledge


# In[2]:


def info_molecule(molecule):
    num_atoms=len(molecule.atoms)
    name_atoms=[];
    mass_atoms=[];
    for i in range(num_atoms):
        name_atoms.append(molecule.atoms[i].element_name)
        mass_atoms.append(molecule.atoms[i].mass)
    molecular_mass=sum(mass_atoms)
    return name_atoms,mass_atoms,molecular_mass


# In[ ]:


def remove_duplicate(x):
    return list(dict.fromkeys(x))


# In[ ]:


def basis_set(element_symbol):
    #HERE I should define all the elements and all the basis set
    if element_symbol=='H':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='F':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='Cl':
        return "TZV2PX-MOLOPT-GTH"
    elif element_symbol=='I':
        return "DZVP-MOLOPT-SR-GTH"
        


# In[ ]:


def potential(element_symbol,functional):
    return "GTH-"+functional
        


# In[ ]:


def optimize_molecule(molecule,functional,length,dir):
    molecule.save('mol_unopt_coord.xyz',overwrite='True')
    with open('mol_unopt_coord.xyz', 'r') as fin:
        data = fin.read().splitlines(True)
    with open('mol_unopt_coord.xyz', 'w') as fout:
        fout.writelines(data[2:])
    molecule=molecule.to_parmed()
    atom_list,mass_list,total_mass=info_molecule(molecule)
    unique_atom_list=remove_duplicate(atom_list)
    num_atoms=len(atom_list)
    num_unique_atoms=len(unique_atom_list)

    mySim = sim.SIM()

    mySim.GLOBAL.RUN_TYPE = "GEO_OPT"
    mySim.GLOBAL.PROJECT  = "molecule_opt"
    mySim.GLOBAL.PRINT_LEVEL = "LOW"
    #FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD='QUICKSTEP'
    mySim.FORCE_EVAL.SUBSYS.CELL.ABC='{L} {L} {L}'.format(L=2*10*length)
    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD='mol_unopt_coord.xyz'
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms);
    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].SECTION_PARAMETERS=unique_atom_list[i]
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].BASIS_SET=basis_set(unique_atom_list[i])
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].POTENTIAL=potential(unique_atom_list[i],functional)

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME=dir+'BASIS_MOLOPT'
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME=dir+'GTH_POTENTIALS'
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT=1E-10

    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF=400
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF=40
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS=5

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE='PAIR_POTENTIAL'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE='DFTD3'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME='dftd3.dat'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF=11

    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS='ATOMIC'
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF=200
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF=1E-6

    mySim.MOTION.GEO_OPT.TYPE='MINIMIZATION'
    mySim.MOTION.GEO_OPT.OPTIMIZER='BFGS'
    mySim.MOTION.GEO_OPT.MAX_ITER=200
    mySim.MOTION.GEO_OPT.MAX_DR=1e-3

    mySim.MOTION.CONSTRAINT.FIXED_ATOMS.LIST ='1'
    mySim.write_changeLog(fn="mol_opt-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn='mol_opt.inp')


# In[7]:


def run_md_pre(molecule,functional,temperature,box_length,dir,number_of_molecules,simulation_time,project_name,table):
    current_molecule=mb.clone(molecule)
    molecule_pmd=mb.clone(current_molecule);
    molecule_pmd=molecule_pmd.to_parmed()
    atom_list,mass_list,total_mass=info_molecule(molecule_pmd)
    unique_atom_list=remove_duplicate(atom_list)
    num_atoms=len(atom_list)
    num_unique_atoms=len(unique_atom_list)
    
    box = mb.Box(lengths=[box_length,box_length,box_length])
    current_molecule.xyz=table;
    box_of_molecule= mb.fill_box(compound=current_molecule, n_compounds=number_of_molecules, box=box)
    filename=project_name+".xyz"
    box_of_molecule.save(filename)
    with open(project_name+".xyz", 'r') as fin:
        data = fin.read().splitlines(True)
    with open(project_name+".xyz", 'w') as fout:
        fout.writelines(data[2:])
        
        
    mySim = sim.SIM()

    mySim.GLOBAL.RUN_TYPE = "MD"
    mySim.GLOBAL.PROJECT  = project_name+"pre"
    mySim.GLOBAL.PRINT_LEVEL = "LOW"


    #FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD='QUICKSTEP'
    mySim.FORCE_EVAL.STRESS_TENSOR='ANALYTICAL';

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME=dir+'BASIS_MOLOPT'
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME=dir+'GTH_POTENTIALS'
    mySim.FORCE_EVAL.DFT.CHARGE=0
    mySim.FORCE_EVAL.DFT.MULTIPLICITY=1
    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF=400
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF=40
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS=4
    mySim.FORCE_EVAL.DFT.QS.METHOD='GPW'
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT=1E-8
    mySim.FORCE_EVAL.DFT.QS.EXTRAPOLATION='ASPC'
    mySim.FORCE_EVAL.DFT.POISSON.PERIODIC="XYZ"
    mySim.FORCE_EVAL.DFT.PRINT.E_DENSITY_CUBE.SECTION_PARAMETERS="OFF"
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS='ATOMIC'
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF=30
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF=1E-6

    mySim.FORCE_EVAL.DFT.SCF.OT.SECTION_PARAMETERS=".TRUE."
    mySim.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER="FULL_SINGLE_INVERSE"
    mySim.FORCE_EVAL.DFT.SCF.OT.MINIMIZER="DIIS"
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.SECTION_PARAMETERS='.TRUE.'

    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.MAX_SCF=10
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.EPS_SCF=1e-6
    mySim.FORCE_EVAL.DFT.SCF.PRINT.RESTART.SECTION_PARAMETERS='OFF'
    mySim.FORCE_EVAL.DFT.SCF.PRINT.DM_RESTART_WRITE='.TRUE.'

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE='PAIR_POTENTIAL'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE='DFTD3'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME='dftd3.dat'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF=11
    
    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD=project_name+".xyz";
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms);
    
    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].SECTION_PARAMETERS=unique_atom_list[i]
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].BASIS_SET=basis_set(unique_atom_list[i])
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].POTENTIAL=potential(unique_atom_list[i],functional)
    mySim.FORCE_EVAL.SUBSYS.CELL.ABC='{L} {L} {L}'.format(L=10*box_length)

    #MOTION SECTION
    mySim.MOTION.GEO_OPT.OPTIMIZER='BFGS'
    mySim.MOTION.GEO_OPT.MAX_ITER=100
    mySim.MOTION.GEO_OPT.MAX_DR=0.003


    mySim.MOTION.MD.ENSEMBLE = "NVT"
    mySim.MOTION.MD.STEPS  = 100
    mySim.MOTION.MD.TIMESTEP = 0.2
    
    mySim.MOTION.MD.TEMPERATURE = temperature
    mySim.MOTION.MD.THERMOSTAT.TYPE = "NOSE"
    mySim.MOTION.MD.THERMOSTAT.REGION = "MASSIVE"
    mySim.MOTION.MD.THERMOSTAT.NOSE.LENGTH = 3
    mySim.MOTION.MD.THERMOSTAT.NOSE.YOSHIDA = 3
    mySim.MOTION.MD.THERMOSTAT.NOSE.TIMECON = 1000.0
    mySim.MOTION.MD.THERMOSTAT.NOSE.MTS = 2


    #mySim.MOTION.MD.PRINT.ENERGY.EACH.MD = 20
    #mySim.MOTION.MD.PRINT.PROGRAM_RUN_INFO.EACH.MD = 20
    #mySim.MOTION.MD.AVERAGES.SECTION_PARAMETERS= ".falbmbjse."


    mySim.MOTION.PRINT.STRESS.SECTION_PARAMETERS='OFF'
    mySim.MOTION.PRINT.TRAJECTORY.EACH.MD=10
    mySim.MOTION.PRINT.VELOCITIES.SECTION_PARAMETERS='OFF'
    mySim.MOTION.PRINT.FORCES.SECTION_PARAMETERS="OFF"
    mySim.MOTION.PRINT.RESTART_HISTORY.SECTION_PARAMETERS="ON"
    mySim.MOTION.PRINT.RESTART_HISTORY.EACH.MD=500
    mySim.MOTION.PRINT.RESTART.SECTION_PARAMETERS="ON"
    mySim.MOTION.PRINT.RESTART.BACKUP_COPIES=3

    mySim.MOTION.PRINT.RESTART.EACH.MD=1

    
    mySim.write_changeLog(fn="md-pre-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn='md-pre.inp')


# In[8]:


def run_md_main(molecule,functional,temperature,box_length,dir,number_of_molecules,simulation_time,project_name):
    
    current_molecule=mb.clone(molecule)
    current_molecule=current_molecule.to_parmed()
    atom_list,mass_list,total_mass=info_molecule(current_molecule)
    unique_atom_list=remove_duplicate(atom_list)
    num_atoms=len(atom_list)
    num_unique_atoms=len(unique_atom_list)
    total_atoms=num_atoms*number_of_molecules;
    #string="tail -{} {}-pos-1.xyz > {}.xyz".format(total_atoms,project_name+"pre",project_name)
    #subprocess.call(string,shell=True)
    
    #we need to decide time_step
    lightest=min(mass_list);
    
    if lightest <1.5:
        time_step=0.5
    elif (lightest>=1.5) and (lightest<40):
        time_step=1
    if lightest>=40:
        time_step=1.5
    steps= round(simulation_time*1000/time_step)
    mySim = sim.SIM()

    mySim.GLOBAL.RUN_TYPE = "MD"
    mySim.GLOBAL.PROJECT  = project_name
    mySim.GLOBAL.PRINT_LEVEL = "LOW"


    #FORCE EVAL SECTION
    mySim.FORCE_EVAL.METHOD='QUICKSTEP'
    mySim.FORCE_EVAL.STRESS_TENSOR='ANALYTICAL';

    mySim.FORCE_EVAL.DFT.BASIS_SET_FILE_NAME=dir+'BASIS_MOLOPT'
    mySim.FORCE_EVAL.DFT.POTENTIAL_FILE_NAME=dir+'GTH_POTENTIALS'
    mySim.FORCE_EVAL.DFT.CHARGE=0
    mySim.FORCE_EVAL.DFT.MULTIPLICITY=1
    mySim.FORCE_EVAL.DFT.MGRID.CUTOFF=400
    mySim.FORCE_EVAL.DFT.MGRID.REL_CUTOFF=40
    mySim.FORCE_EVAL.DFT.MGRID.NGRIDS=4
    mySim.FORCE_EVAL.DFT.QS.METHOD='GPW'
    mySim.FORCE_EVAL.DFT.QS.EPS_DEFAULT=1E-8
    mySim.FORCE_EVAL.DFT.QS.EXTRAPOLATION='ASPC'
    mySim.FORCE_EVAL.DFT.POISSON.PERIODIC="XYZ"
    mySim.FORCE_EVAL.DFT.PRINT.E_DENSITY_CUBE.SECTION_PARAMETERS="OFF"
    mySim.FORCE_EVAL.DFT.SCF.SCF_GUESS='ATOMIC'
    mySim.FORCE_EVAL.DFT.SCF.MAX_SCF=30
    mySim.FORCE_EVAL.DFT.SCF.EPS_SCF=1E-6

    mySim.FORCE_EVAL.DFT.SCF.OT.SECTION_PARAMETERS=".TRUE."
    mySim.FORCE_EVAL.DFT.SCF.OT.PRECONDITIONER="FULL_SINGLE_INVERSE"
    mySim.FORCE_EVAL.DFT.SCF.OT.MINIMIZER="DIIS"
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.SECTION_PARAMETERS='.TRUE.'

    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.MAX_SCF=10
    mySim.FORCE_EVAL.DFT.SCF.OUTER_SCF.EPS_SCF=1e-6
    mySim.FORCE_EVAL.DFT.SCF.PRINT.RESTART.SECTION_PARAMETERS='OFF'
    mySim.FORCE_EVAL.DFT.SCF.PRINT.DM_RESTART_WRITE='.TRUE.'

    mySim.FORCE_EVAL.DFT.XC.XC_FUNCTIONAL.SECTION_PARAMETERS=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.POTENTIAL_TYPE='PAIR_POTENTIAL'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.TYPE='DFTD3'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.PARAMETER_FILE_NAME='dftd3.dat'
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.REFERENCE_FUNCTIONAL=functional
    mySim.FORCE_EVAL.DFT.XC.VDW_POTENTIAL.PAIR_POTENTIAL.R_CUTOFF=11
    mySim.FORCE_EVAL.SUBSYS.COORD.DEFAULT_KEYWORD=project_name+".xyz";
    mySim.FORCE_EVAL.SUBSYS.init_atoms(num_atoms);
    
    for i in range(num_unique_atoms):
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].SECTION_PARAMETERS=unique_atom_list[i]
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].BASIS_SET=basis_set(unique_atom_list[i])
        mySim.FORCE_EVAL.SUBSYS.KIND[i+1].POTENTIAL=potential(unique_atom_list[i],functional)
    mySim.FORCE_EVAL.SUBSYS.CELL.ABC='{L} {L} {L}'.format(L=10*box_length)

    #MOTION SECTION
    mySim.MOTION.GEO_OPT.OPTIMIZER='BFGS'
    mySim.MOTION.GEO_OPT.MAX_ITER=100
    mySim.MOTION.GEO_OPT.MAX_DR=0.003


    mySim.MOTION.MD.ENSEMBLE = "NVT"
    mySim.MOTION.MD.STEPS  = steps
    mySim.MOTION.MD.TIMESTEP = time_step
    mySim.MOTION.MD.TEMPERATURE = temperature
    mySim.MOTION.MD.THERMOSTAT.TYPE = "NOSE"
    mySim.MOTION.MD.THERMOSTAT.REGION = "MASSIVE"
    mySim.MOTION.MD.THERMOSTAT.NOSE.LENGTH = 3
    mySim.MOTION.MD.THERMOSTAT.NOSE.YOSHIDA = 3
    mySim.MOTION.MD.THERMOSTAT.NOSE.TIMECON = 1000.0
    mySim.MOTION.MD.THERMOSTAT.NOSE.MTS = 2


    #mySim.MOTION.MD.PRINT.ENERGY.EACH.MD = 20
    #mySim.MOTION.MD.PRINT.PROGRAM_RUN_INFO.EACH.MD = 20
    #mySim.MOTION.MD.AVERAGES.SECTION_PARAMETERS= ".falbmbjse."


    mySim.MOTION.PRINT.STRESS.SECTION_PARAMETERS='OFF'
    mySim.MOTION.PRINT.TRAJECTORY.EACH.MD=10
    mySim.MOTION.PRINT.VELOCITIES.SECTION_PARAMETERS='OFF'
    mySim.MOTION.PRINT.FORCES.SECTION_PARAMETERS="OFF"
    mySim.MOTION.PRINT.RESTART_HISTORY.SECTION_PARAMETERS="ON"
    mySim.MOTION.PRINT.RESTART_HISTORY.EACH.MD=500
    mySim.MOTION.PRINT.RESTART.SECTION_PARAMETERS="ON"
    mySim.MOTION.PRINT.RESTART.BACKUP_COPIES=3

    mySim.MOTION.PRINT.RESTART.EACH.MD=1

    
    mySim.write_changeLog(fn="md-main-changeLog.out")
    mySim.write_errorLog()
    mySim.write_inputFile(fn='md-main.inp')


# In[ ]:




