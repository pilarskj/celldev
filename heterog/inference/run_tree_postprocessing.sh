#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=10240
#SBATCH --job-name=postprocessing_TiDe_hierarchical
#SBATCH --array=1-20

# use (TiDe|Typewriter)_(distinct|hierarchical)
# for distinct: seeds 1,2,5,6,7,9,11,13,16,19,21,23,24,27,31,32,33,34,35,36

# set paths
cd ~/Projects/celldev
dir="/Users/jpilarski/Projects/celldev_data/heterog" # make sure that output directories exist, if not, make them

# recorder
method="TiDe"
xmlFile="heterog/inference/mapping_${method}.xml"

# prototype of multitype model
prototype="hierarchical"

# inference id
seed=$SLURM_ARRAY_TASK_ID


# for each .trees file
if [ -e ${dir}/${method}/${prototype}/inferenceOutput/tree_${seed}.inference.trees ]
then

  # map types on trees
  java -Xmx15g -jar software/simbundle.jar \
  -version_file software/beast2_version.xml \
  -version_file software/feast_version.xml \
  -version_file software/tidetree_version.xml \
  -version_file software/sciphy_version.xml \
  -version_file software/bdmmprime_version.xml \
  -overwrite -seed ${seed} -D "dir=${dir},prototype=${prototype}" \
  ${xmlFile}

  # construct 95% credible set of trees
  java -Xmx15g -cp software/simbundle.jar beastfx.app.tools.TreeTraceAnalysis -trees ${dir}/${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees -out ${dir}/${method}/${prototype}/inferenceOutput/tree_${seed}.hpd.trees

  # construct MCC tree
  $HOME/beast/bin/treeannotator -lowMem true -height median ${dir}/${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees ${dir}/${method}/${prototype}/inferenceOutput/tree_${seed}.mcc.tree

fi

# job submission (Slurm): sbatch < run_tree_postprocessing.sh
# Note: make sure that Java is installed/ module is loaded, the location of .jar is correct and an installation of BEAST2 with Tree annotator is available
