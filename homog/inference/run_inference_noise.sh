#!/bin/bash
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=1024
#SBATCH --job-name=inference_TiDe_noise
#SBATCH --array=1-100
#SBATCH --output=/cluster/scratch/jpilarski/celldev_data/homog/outs/%x_%a.out

method="TiDe"
inDir="/cluster/scratch/jpilarski/celldev_data/homog/$method/baseline_noise/simulationOutput"
#inDir="~/Projects/celldev_data/homog/TiDe/baseline_noise/simulationOutput"
statsFile="$inDir/filtering_stats.csv"
xmlFile="/cluster/home/jpilarski/celldev/homog/inference/inference_$method.xml"
#xmlFile="~/Projects/celldev/homog/inference/inference_$method.xml"
outDir="/cluster/scratch/jpilarski/celldev_data/homog/$method/baseline_noise/inferenceOutput"
#outDir="$HOME/Projects/celldev_data/homog/$method/baseline_noise/inferenceOutput"

# create directory if not available
[ -d "$outDir" ] || mkdir -p "$outDir"

# find setting
i=$SLURM_ARRAY_TASK_ID

# read all information from file
read tree seed nTargets taxLabel samplingProp <<< "$(
awk -F',' -v i="$i" '
NR==1 { for (j=1;j<=NF;j++) col[$j]=j; next }
NR==i+1 { print $(col["tree_model"]), $(col["seed"]), $(col["num_targets"]), $(col["tax_label"]), $(col["sampling_prop"]) }
' $statsFile
)"

#echo "run=$i, tree=$tree, seed=$seed, nTargets=$nTargets, taxLabel=$taxLabel, rho=$samplingProp" > ${outDir}/${tree}_${seed}.inference.out

if [ $method = "TiDe" ]; then
    D="tree=$tree,taxLabel=$taxLabel,rho=$samplingProp,inDir=$inDir,outDir=$outDir,scarringHeight=40,scarringDuration=40"
elif [ $method = "Typewriter" ]; then
    D="tree=$tree,taxLabel=$taxLabel,rho=$samplingProp,inDir=$inDir,outDir=$outDir,nTapes=$nTargets,tapeLength=5"
fi

echo $D

# run inference
java -Dglass.platform=Monocle -Dmonocle.platform=Headless --module-path=$HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
-jar software/simbundle.jar \
-version_file software/beast2_version.xml \
-version_file software/feast_version.xml \
-version_file software/tidetree_version.xml \
-version_file software/sciphy_version.xml \
-version_file software/bdsky_version.xml \
-statefile ${outDir}/${tree}_${seed}.inference.state \
-overwrite \
-seed $seed \
-D $D \
${xmlFile} >> ${outDir}/${tree}_${seed}.inference.out
