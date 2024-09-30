#!/bin/bash

# seed of alignment
seed=$SLURM_ARRAY_TASK_ID # Slurm
# seed=$LSB_JOBINDEX # LSF

# tree and dirs
while getopts "t:i:o:h:d:" flag
do
    case ${flag} in
        t) tree=${OPTARG};;
        i) inDir=${OPTARG};;
        o) outDir=${OPTARG};;
        h) scarringHeight=${OPTARG};;
        d) scarringDuration=${OPTARG};;
    esac
done

if [[ "tree_s tree_sd" =~ ${tree} ]]
then
  rho=1
else
  rho=0.1
fi

# run inference
echo "seed=${seed}, tree=${tree}, rho=${rho}" > ${outDir}/${tree}_${seed}.inference.out

java -jar $HOME/beasts2.7.jar -statefile ${outDir}/${tree}_${seed}.inference.xml.state -overwrite \
-seed ${seed} -D "tree=${tree},rho=${rho},inDir=${inDir},outDir=${outDir},scarringHeight=${scarringHeight},scarringDuration=${scarringDuration}" \
inference_TiDe.xml >> ${outDir}/${tree}_${seed}.inference.out
