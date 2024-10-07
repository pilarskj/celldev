#!/bin/bash

# set paths
cd ~/Projects/celldev_data/heterog # make sure that output directories exist, if not, make them
xmlFile='~/Projects/celldev/heterog/simulation_barcodes/simulation_TiDe.xml' # or Typewriter, adapt protopype in xml's

# run simulation of alignments
for seed in `seq 1 20` # or 1 2 5 6 7 9 11 13 16 19 21 23 24 27 31 32 33 34 35 36 for TiDeTree alignments
do
  java -jar $HOME/beasts2.7.jar -overwrite -seed ${seed} ${xmlFile}
done

# command: bash run_simulation.sh
