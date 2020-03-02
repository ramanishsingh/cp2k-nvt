import flow
import os
from shutil import copyfile
from subprocess import call
import shutil
class Project(flow.FlowProject):
    pass

#@Project.label
#def melted(job):
#    return job.isfile('melt.run1a.dat')

@Project.label
def has_input_files(job):
    return (job.isfile('mol_opt.inp') and job.isfile('mol_unopt_coord.xyz'))

@Project.label
def struc_optimized(job):
    return (job.isfile("mol_opt.out"))


# For equilibration, only run out of prod1 directory
#@Project.label
#def prod1(job):
#    return(job.statepoint()['Prod']==1)


# This project operation will run the melting stage of the VLE
@Project.operation
@Project.post(struc_optimized)
@Project.pre(has_input_files)
#@Project.pre(prod1)
def run_mol_opt(job):
    # Information for running slurm job
    home  = job.workspace()
    base=os.getcwd()

    os.chdir(home)
    print("starting cp2k")
    #call("~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i md.inp -o md.out",shell=True)
    call("mpirun -np 24 ~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i mol_opt.inp -o mol_opt.out",shell=True)
    print("Complete!")
    copyfile('{}/molecule_opt-pos-1.xyz'.format(home),'{}/molecule_opt-pos-1.xyz'.format(base))
    # Copy output files to avoid overwriting
   # copyfile("{}/config1a.dat".format(scr),job.fn("melt.config1a.dat"))
    #copyfile("{}/config1a.dat".format(scr),job.fn("fort.77"))

if __name__ == '__main__':
    Project().main()
