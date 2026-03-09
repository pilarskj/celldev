#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=processing_Typewriter_baseline_large
#SBATCH --array=1-20
#SBATCH --output=/cluster/scratch/jpilarski/celldev_data/homog/outs/%x_%a.out

# postprocessing all .inference.trees files for one experimental setup
# make sure that you have java installed/ a java module loaded!
# beside the .jar file, a local installation of BEAST2 with Tree annotator is necessary
# if run locally, loop over the seeds (for seed in `seq 1 20` do ... done)

# seed of alignment
seed=$SLURM_ARRAY_TASK_ID 

# tree and dirs
# run script from celldev project directory
dir="/cluster/scratch/jpilarski/celldev_data/homog/Typewriter/baseline_large" # here and in header, adapt path

# for each tree
for tree in tree_s tree_ss tree_sd tree_sds tree_bd
do
  # for each .trees file
  if [ -e ${dir}/inferenceOutput/${tree}_${seed}.inference.trees ]
  then
  
    # thin .trees file (otherwise exceeds heap space) - in case of large trees, otherwise used original .trees file
    $HOME/beast/bin/logcombiner -b 0 -resample 25000 -log ${dir}/inferenceOutput/${tree}_${seed}.inference.trees -o ${dir}/inferenceOutput/${tree}_${seed}.thinned.trees

    # construct 95% credible set of trees
    java -Xmx8g -cp software/simbundle.jar beastfx.app.tools.TreeTraceAnalysis -trees ${dir}/inferenceOutput/${tree}_${seed}.thinned.trees -out ${dir}/inferenceOutput/${tree}_${seed}.hpd.trees

    # construct MCC tree
    $HOME/beast/bin/treeannotator -lowMem true -height median ${dir}/inferenceOutput/${tree}_${seed}.thinned.trees ${dir}/inferenceOutput/${tree}_${seed}.mcc.tree

    echo "${tree}_${seed}: done" 
  fi
done

# check file numbers (should be the same as *.inference.trees)
#ls *.hpd.trees | wc -l
#ls *.mcc.tree | wc -l

# afterwards: compress .inference.trees files using 
#sbatch --job-name=compress --wrap='gzip --best *.inference.trees'
