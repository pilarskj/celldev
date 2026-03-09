#!/bin/bash

# set paths
code_dir="/Users/jpilarski/Projects/celldev"
data_dir="/Users/jpilarski/Projects/celldev_data/heterog"

cd ${code_dir}

# run simulation of trees (adapt the number)
## make sure that output directories exists and modify the prototype (distinct vs. hierarchical) as needed!
for seed in `seq 1 20` 
do
  java -jar software/simbundle.jar \
  -version_file software/beast2_version.xml \
  -version_file software/feast_version.xml \
  -version_file software/bdmmprime_version.xml \
  -overwrite \
  -seed ${seed} \
  -D "dir=${data_dir}" \
  heterog/simulation_trees/simulation_trees.xml
done