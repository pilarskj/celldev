# run simulation of trees or alignments
for seed in `seq 1 25` # or 1 2 5 6 7 9 11 12 13 16 19 21 22 23 24 27 31 32 33 34 35 36 37 38 39 for TiDeTree alignments
do
  java -jar ../beasts2.7.jar -overwrite -seed ${seed} simulation_trees.xml
done

# command: sh run_simulation.sh
