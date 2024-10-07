#!/bin/bash

# -----
# Script for postprocessing all .inference.trees files for one experimental setup
# -----

# make sure that you have java installed/ a java module loaded!
# beside the beasts2.7.jar file, a local installation of BEAST2 with Tree annotator is necessary

# example job submission on cluster
# Slurm:
#OMP_NUM_THREADS=1 sbatch --time=120 --mem-per-cpu=10240 --job-name='postprocessing' --array=1-20 --wrap='bash run_tree_postprocessing.sh -d dir $SCRATCH/celldev_data/Typewriter/baseline'
# LSF:
#OMP_NUM_THREADS=1 bsub -W 48:00 -R 'rusage[mem=10240]' -J 'postprocessing_Typewriter_baseline[1-20]' "bash run_tree_postprocessing.sh -d $SCRATCH/celldev_data/Typewriter/baseline"

# if run locally, loop over the seeds (for seed in `seq 1 20` do ... done)

# seed of alignment
seed=$SLURM_ARRAY_TASK_ID # Slurm
# seed=$LSB_JOBINDEX # LSF

# tree and dirs
while getopts "d:" flag
do
    case ${flag} in
        d) dir=${OPTARG};;
    esac
done

# for each tree
for tree in tree_s tree_ss tree_sd tree_sds tree_bd
do
  # for each .trees file
  if [ -e ${dir}/inferenceOutput/${tree}_${seed}.inference.trees ]
  then

    # construct 95% credible set of trees
    java -Xmx15g -cp $HOME/beasts2.7.jar beastfx.app.tools.TreeTraceAnalysis -trees ${dir}/inferenceOutput/${tree}_${seed}.inference.trees -out ${dir}/inferenceOutput/${tree}_${seed}.hpd.trees

    # construct MCC tree
    $HOME/beast/bin/treeannotator -heights median ${dir}/inferenceOutput/${tree}_${seed}.inference.trees ${dir}/inferenceOutput/${tree}_${seed}.mcc.tree

    echo "${dir}/${tree}_${seed}: done" >> postprocessing.out
  fi
done

# check file numbers (should be the same as *.inference.trees)
#ls /path/to/inferenceOutput/*.hpd.trees | wc -l
#ls /path/to/inferenceOutput/*.mcc.tree | wc -l

# afterwards: compress .inference.trees files using 
# Slurm:
#sbatch --job-name=compress --wrap='gzip --best /path/to/inferenceOutput/*.inference.trees'
# LSF:
#bsub -J 'compress' "gzip --best /path/to/inferenceOutput/*.inference.trees"
