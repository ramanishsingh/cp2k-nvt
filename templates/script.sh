
{% extends base_script %}
{% block project_header %}
#PBS -l mem=50gb,nodes=1:ppn=10

module purge
module load mkl
module load fftw
module load intel/cluster/2018

date >> execution.log
{{ super() }}
{% endblock %}
