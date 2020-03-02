import flow
import os
from shutil import copyfile
from subprocess import call
from subprocess import Popen, PIPE

import shutil
class Project(flow.FlowProject):
    pass

#@Project.label
#def melted(job):
#    return job.isfile('melt.run1a.dat')

@Project.label
def has_input_files(job):
    return (job.isfile('md-pre.inp') and job.isfile('md-main.inp') )

@Project.label
def struc_optimized(job):
    return (job.isfile("md-main.out") and job.isfile("md-pre.out"))


# For equilibration, only run out of prod1 directory
#@Project.label
#def prod1(job):
#    return(job.statepoint()['Prod']==1)


# This project operation will run the melting stage of the VLE
@Project.operation
@Project.post(struc_optimized)
@Project.pre(has_input_files)
#@Project.pre(prod1)
def run_md(job):
    # Information for running slurm job
    home  = job.workspace()
    base=os.getcwd()

    os.chdir(home)
    print("starting cp2k")
    #call("~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i md.inp -o md.out",shell=True)
    call("mpirun -np 10 ~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i md-pre.inp -o md-pre.out",shell=True)
    print("Complete!")
    process = Popen(['grep', '@INCLUDE', 'md-pre.inp'], stdout=PIPE, stderr=PIPE)
    coor_file= process.communicate()[0].decode("utf-8").strip().split()[1]
    process = Popen(['wc', '-l', coor_file], stdout=PIPE, stderr=PIPE)
    number_of_atoms= int(process.communicate()[0].decode("utf-8").strip().split()[0])
    process = Popen(['grep', 'PROJECT', 'md-pre.inp'], stdout=PIPE, stderr=PIPE)
    project_name= process.communicate()[0].decode("utf-8").strip().split()[1]

    string="tail -{} {}-pos-1.xyz > {}".format(number_of_atoms,project_name,coor_file)
    call(string, shell=True) 
    call("mpirun -np 10 ~/test-cp2k/cp2k/exe/Linux-x86-64-intel/cp2k.popt -i md-main.inp -o md-main.out",shell=True)
    
    # Copy output files to avoid overwriting
   # copyfile("{}/config1a.dat".format(scr),job.fn("melt.config1a.dat"))
    #copyfile("{}/config1a.dat".format(scr),job.fn("fort.77"))

if __name__ == '__main__':
    Project().main()
