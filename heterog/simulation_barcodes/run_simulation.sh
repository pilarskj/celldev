#!/bin/bash

# set paths
cd ~/Projects/celldev
dataDir="/Users/jpilarski/Projects/celldev_data/heterog"
# make sure that output directories exist, if not, make them

xmlFile="heterog/simulation_barcodes/simulation_TiDe.xml" # or Typewriter
trans="hierarchical"

# run simulation of alignments
for seed in `seq 1 20` # or 1 2 5 6 7 9 11 13 16 19 21 23 24 27 31 32 33 34 35 36 for "distinct" alignments
do
  java -jar software/simbundle.jar \
  -version_file software/beast2_version.xml \
  -version_file software/feast_version.xml \
  -version_file software/tidetree_version.xml \
  -version_file software/sciphy_version.xml \
  -overwrite \
  -seed ${seed} \
  -D "dir=${dataDir},trans=${trans}" \
  ${xmlFile}
done

# command: bash run_simulation.sh
