#!/bin/bash

rm -rf runs
mkdir runs

for rep in {1..20}
do

    # base_dir = "/global/scratch/freyman/projects/draba/analyses/divergence_times/"
    # run_num = 1
    # n_iterations = 1000
    # seed(1) 
    
    echo "base_dir = \"/global/scratch/freyman/projects/draba/analyses/divergence_times/\"" >> runs/$rep.Rev
    echo "run_num = $rep" >> runs/$rep.Rev
    echo "n_iterations = 1000" >> runs/$rep.Rev
    echo "seed($rep)" >> runs/$rep.Rev
    cat draba.Rev >> runs/$rep.Rev
    
done


echo "#!/bin/bash
# Job name:
#SBATCH --job-name=draba
#
# Partition:
#SBATCH --partition=savio
#
# Processors:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#
# Wall clock limit (max=72:00:00)
#SBATCH --time=72:00:00
#
# Account:
#SBATCH --account=fc_rothfelslab
#
# Specify Faculty Computing Allowance:
#SBATCH --qos=savio_normal
#" > runs/draba.sh
for core in {1..20}
do
    echo "/global/scratch/freyman/revbayes/projects/cmake/rb /global/scratch/freyman/projects/draba/analyses/divergence_times/runs/$core.Rev > /global/scratch/freyman/projects/draba/analyses/divergence_times/output/$core.out &" >> runs/draba.sh
done
echo "wait" >> runs/draba.sh

