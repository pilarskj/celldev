# -----
# Script for postprocessing all .inference.trees files for one experimental setup
# -----

# example command: OMP_NUM_THREADS=1 bsub -W 48:00 -R 'rusage[mem=10240]' -J 'postprocessing_Typewriter_baseline[1-20]' "sh run_tree_postprocessing.sh -s Typewriter/baseline"

# seed of alignment
seed=$LSB_JOBINDEX

module load openjdk/17.0.0_35

# tree and dirs
while getopts "s:" flag
do
    case ${flag} in
        s) setup=${OPTARG};;
    esac
done

# for each tree
for tree in tree_s tree_ss tree_sd tree_sds tree_bd
do
  # for each .trees file
  if [ -e ${setup}/inferenceOutput/${tree}_${seed}.inference.trees ]
  then

    # construct 95% credible set of trees
    java -Xmx15g -cp $HOME/beasts2.7.jar beastfx.app.tools.TreeTraceAnalysis -trees ${setup}/inferenceOutput/${tree}_${seed}.inference.trees -out ${setup}/inferenceOutput/${tree}_${seed}.hpd.trees

    # construct MCC tree
    $HOME/beast/bin/treeannotator -heights median ${setup}/inferenceOutput/${tree}_${seed}.inference.trees ${setup}/inferenceOutput/${tree}_${seed}.mcc.tree

    echo "${setup}/${tree}_${seed}: done" >> postprocessing.out
  fi
done
