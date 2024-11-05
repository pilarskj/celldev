#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=5120
#SBATCH --job-name=inference_TiDe_hierarchical
#SBATCH --array=1-20

# use (TiDe|Typewriter)_(distinct|hierarchical)
# for distinct: seeds 1,2,5,6,7,9,11,13,16,19,21,23,24,27,31,32,33,34,35,36

# set paths
cd ~/Projects/celldev_data/heterog # make sure that output directories exist, if not, make them
xmlFile='~/Projects/celldev/heterog/inference/inference_TiDe.xml'

# inference id
seed=$SLURM_ARRAY_TASK_ID

# recorder
method="TiDe"

# prototype of multitype model
prototype="hierarchical"

# run inference
java -jar $HOME/beasts2.7.jar -statefile ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.xml.state -overwrite \
-seed ${seed} -D "prototype=${prototype}" ${xmlFile} >> ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.out

# job submission (Slurm): sbatch < run_inference.sh
# Note: make sure that Java is installed/ module is loaded, and the location of beasts2.7.jar is correct