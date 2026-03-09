#!/bin/bash

# run from ~/Projects/celldev
dataDir="/Users/jpilarski/Projects/celldev_data/homog"

treeDir="$dataDir/LargeTrees"
outDir="$dataDir/Typewriter/baseline_noise/simulationOutputRaw" # or TiDe

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
    java -jar software/simbundle.jar \
    -version_file software/beast2_version.xml \
    -version_file software/feast_version.xml \
    -version_file software/tidetree_version.xml \
    -version_file software/sciphy_version.xml \
    -overwrite \
    -seed $seed \
    -D "treeDir=$treeDir,outDir=$outDir,tree=$tree,seed=$seed" \
    homog/simulation_barcodes/simulation_Typewriter_noise.xml # or TiDe
  done
done

# for TiDe: -D "treeFile=$treeFile,outFile=$outFile" 
