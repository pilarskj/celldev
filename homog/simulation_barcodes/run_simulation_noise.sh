#!/bin/bash

# run from ~/Projects/celldev/homog/simulation_barcodes

treeDir="/Users/jpilarski/Projects/celldev_data/homog/LargeTrees"
outDir="/Users/jpilarski/Projects/celldev_data/homog/Typewriter/baseline_noise/simulationOutputRaw" # or TiDe

# create directory if not available
[ -d "$outDir" ] || mkdir -p "$outDir"

for tree in tree_{s,ss,sd,sds,bd}
do
  for seed in `seq 1 20`
  do
    #treeFile="$treeDir/${tree}_${seed}.newick"
    #outFile="$outDir/${tree}_${seed}.nexus"
    #echo "treeFile=$treeFile,outFile=$outFile"
    echo "treeDir=$treeDir,outDir=$outDir,tree=$tree,seed=$seed"
    java -jar /Users/jpilarski/simbundle.jar \
    -version_file /Users/jpilarski/intellij/SciPhy/version.xml \
    -version_file /Users/jpilarski/intellij/tidetree-dropout/version.xml \
    -version_file /Users/jpilarski/intellij/feast/version.xml \
    -version_file /Users/jpilarski/intellij/Beast2/version.xml \
    -overwrite \
    -seed $seed \
    -D "treeDir=$treeDir,outDir=$outDir,tree=$tree,seed=$seed" \
    simulation_Typewriter_noise.xml # or TiDe
  done
done

# for TiDe: -D "treeFile=$treeFile,outFile=$outFile" 
