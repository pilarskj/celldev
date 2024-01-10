#!/bin/sh
#BSUB -W 120:00
#BSUB -R "rusage[mem=5120]"
#BSUB -J "inference_TiDe_hierarchical[1-25]"
# or BSUB -J "inference_Typewriter_hierarchical[1-25]"
# or BSUB -J "inference_TiDe_distinct[1,2,5,6,7,9,11,12,13,16,19,21,22,23,24,27,31,32,33,34,35,36,37,38,39]"
# or BSUB -J "inference_Typewriter_distinct[1,2,5,6,7,9,11,12,13,16,19,21,22,23,24,27,31,32,33,34,35,36,37,38,39]"

# inference id
seed=$LSB_JOBINDEX

#recorder
method="TiDe"

# prototype of multitype model
prototype="hierarchical"

module load openjdk/17.0.0_35

# run inference
java -jar $HOME/beasts2.7.jar -statefile ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.xml.state -overwrite \
-seed ${seed} -D "prototype=${prototype}" inference_${method}.xml >> ${method}/${prototype}/inferenceOutput/tree_${seed}.inference.out

# command: bsub < run_inference.sh
