#!/bin/bash

for locus in "nuclear" "chloroplast"
do
    for rep in {1..20}
    do

        # base_dir = "/global/scratch/freyman/projects/draba/analyses/divergence_times/"
        # run_num = 1
        # n_iterations = 1000
        # sample_freq = 1
        # seed(1) 
        
        echo "base_dir = \"/global/scratch/freyman/projects/draba/analyses/divergence_times/\"" >> runs/$locus$rep.Rev
        echo "run_num = \"$locus$rep\"" >> runs/$locus$rep.Rev
        echo "n_iterations = 1000" >> runs/$locus$rep.Rev
        echo "sample_freq = 1" >> runs/$locus$rep.Rev
        echo "seed($rep)" >> runs/$locus$rep.Rev
        cat draba_$locus.Rev >> runs/$locus$rep.Rev
        
    done


echo "#!/bin/bash
# Job name:
#SBATCH --job-name=draba_$locus
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
#" > runs/draba_$locus.sh
for core in {1..20}
do
    echo "/global/scratch/freyman/revbayes/projects/cmake/rb /global/scratch/freyman/projects/draba/analyses/divergence_times/runs/$locus$core.Rev > /global/scratch/freyman/projects/draba/analyses/divergence_times/output/$locus$core.out &" >> runs/draba_$locus.sh
done
echo "wait" >> runs/draba_$locus.sh

done
