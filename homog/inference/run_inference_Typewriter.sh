#!/bin/bash

# seed of alignment
seed=$SLURM_ARRAY_TASK_ID # Slurm
# seed=$LSB_JOBINDEX # LSF

# tree and dirs
while getopts "t:x:i:o:n:l:" flag
do
    case ${flag} in
        t) tree=${OPTARG};;
        x) xmlFile=${OPTARG};;
        i) inDir=${OPTARG};;
        o) outDir=${OPTARG};;
        n) nTapes=${OPTARG};;
        l) tapeLength=${OPTARG};;
    esac
done

if [[ "tree_s tree_sd" =~ ${tree} ]]
then
  rho=1
else
  rho=0.1 # or 0.01 in case of lower sampling
fi

# run inference
echo "seed=${seed}, tree=${tree}, rho=${rho}" > ${outDir}/${tree}_${seed}.inference.out

java -Dglass.platform=Monocle -Dmonocle.platform=Headless --module-path=$HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
-jar software/simbundle.jar \
-version_file software/beast2_version.xml \
-version_file software/feast_version.xml \
-version_file software/tidetree_version.xml \
-version_file software/sciphy_version.xml \
-version_file software/bdsky_version.xml \
-statefile ${outDir}/${tree}_${seed}.inference.state \
-overwrite \
-seed ${seed} \
-D "tree=${tree},rho=${rho},inDir=${inDir},outDir=${outDir},nTapes=${nTapes},tapeLength=${tapeLength}" \
${xmlFile} >> ${outDir}/${tree}_${seed}.inference.out

# Note: make sure that Java is installed/ module is loaded, and the location of .jar is correct
