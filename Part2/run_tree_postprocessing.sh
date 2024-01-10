#!/bin/sh
#BSUB -W 48:00
#BSUB -R "rusage[mem=10240]"
#BSUB -J "pp_TiDe_hierarchical[10]"

# inference id
seed=$LSB_JOBINDEX

#recorder
method="TiDe"

# prototype of multitype model
prototype="hierarchical"

module load openjdk/17.0.0_35

# for each .trees file
if [ -e ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.trees ]
then

  # map types on trees
  java -Xmx15g -jar $HOME/beasts2.7.jar -overwrite -seed ${seed} -D "prototype=${prototype}" mapping_${method}.xml

  # construct 95% credible set of trees
  java -Xmx15g -cp $HOME/beasts2.7.jar beastfx.app.tools.TreeTraceAnalysis -trees ${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees -out ${method}/${prototype}/inferenceOutput/tree_${seed}.hpd.trees

  # construct MCC tree
  $HOME/beast/bin/treeannotator -heights median ${method}/${prototype}/inferenceOutput/tree_${seed}.typed.node.trees ${method}/${prototype}/inferenceOutput/tree_${seed}.mcc.tree

fi
