#!/bin/bash
#SBATCH --time 8:00:00
#SBATCH --mem-per-cpu=10240
#SBATCH --job-name=postprocessing_TiDe_hierarchical
#SBATCH --array=1-20

# use (TiDe|Typewriter)_(distinct|hierarchical)
# for distinct: seeds 1,2,5,6,7,9,11,13,16,19,21,23,24,27,31,32,33,34,35,36

# set paths
cd ~/Projects/celldev_data/heterog 
xmlFile='~/Projects/celldev/heterog/inference/mapping_TiDe.xml'

# inference id
seed=$SLURM_ARRAY_TASK_ID

#recorder
method="TiDe"

# prototype of multitype model
prototype="hierarchical"

# for each .trees file
if [ -e ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.trees ]
then

  # map types on trees
  java -Xmx15g -jar $HOME/beasts2.7.jar -overwrite -seed ${seed} -D "prototype=${prototype}" ${xmlFile}

  # construct 95% credible set of trees
  java -Xmx15g -cp $HOME/beasts2.7.jar beastfx.app.tools.TreeTraceAnalysis -trees ${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees -out ${method}/${prototype}/inferenceOutput/tree_${seed}.hpd.trees

  # construct MCC tree
  $HOME/beast/bin/treeannotator -heights median ${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees ${method}/${prototype}/inferenceOutput/tree_${seed}.mcc.tree

fi

# job submission (Slurm): sbatch < run_tree_postprocessing.sh
# Note: make sure that Java is installed/ module is loaded, the location of beasts2.7.jar is correct and an installation of BEAST2 with Tree annotator is available
