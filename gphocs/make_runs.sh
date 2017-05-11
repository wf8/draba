#!/bin/bash

rm -rf runs
mkdir runs

for rep in {1..20}
do

    #GENERAL-INFO-START
    #seq-file            /global/scratch/freyman/projects/draba/gphocs/data/ramos_filtered4.gphocs
    #trace-file          /global/scratch/freyman/projects/draba/gphocs/output/replicate.log
    #locus-mut-rate      VAR 1.0
    #random-seed         1234

    echo "GENERAL-INFO-START" >> runs/${rep}.ctl
    echo "seq-file            /global/scratch/freyman/projects/draba/gphocs/data/ramos_filtered4.gphocs" >> runs/${rep}.ctl
    echo "trace-file          /global/scratch/freyman/projects/draba/gphocs/output/${rep}.log" >> runs/${rep}.ctl
    echo "locus-mut-rate      VAR 1.0" >> runs/${rep}.ctl
    echo "random-seed         ${rep}" >> runs/${rep}.ctl
    cat template.ctl >> runs/${rep}.ctl
    
done


echo "#!/bin/bash
# Job name:
#SBATCH --job-name=draba-gph
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
    echo "/global/scratch/freyman/G-PhoCS/bin/G-PhoCS-1-2-3 /global/scratch/freyman/projects/draba/gphocs/runs/${core}.ctl > /global/scratch/freyman/projects/draba/gphocs/output/${core}.out &" >> runs/draba.sh
done
echo "wait" >> runs/draba.sh

