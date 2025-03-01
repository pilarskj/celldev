# -----
# Snakemake pipeline for alignment simulations (for homogeneous cell populations) and inference
# -----

# SETUP
import os

# paths (adapt!)
code_dir = '~/Projects/celldev/homog'
data_dir = '~/Projects/celldev_data/homog'
os.chdir(data_dir)

# trees
trees = ['tree_s', 'tree_ss', 'tree_sd', 'tree_sds', 'tree_bd']

# simulation paramaters
nSim = 20
seeds = list(range(1, nSim + 1))

# experimental scenario
method = 'Typewriter'
setup = 'baseline'

# experimental parameters
editRate = 0.05 # [0.01, 0.05, 0.1, 0.45]
if method == 'TiDe':
    nTargets = 20 # [5, 20, 40]
    scarringHeight = 40 # [40, 30, 20]
    scarringDuration = 40 # [40, 20]
elif method == 'Typewriter':
    nTapes = 20 # [1, 5, 20, 40]
    tapeLength = 5 # [2, 5, 10, 20]

# output files
if method == 'TiDe':
    simFiles = expand('TiDe/' + setup + '/simulationOutput/' + '{tree}_{seed}.alignment.nexus', tree = trees, seed = seeds)
    infFiles = expand('TiDe/' + setup + '/inferenceOutput/' + '{tree}_1.inference.out', tree = trees)
elif method == 'Typewriter':
    simFiles = expand('Typewriter/' + setup + '/simulationOutput/' + '{tree}_{seed}.alignment_1.nexus', tree = trees, seed = seeds)
    infFiles = expand('Typewriter/' + setup + '/inferenceOutput/' + '{tree}_1.inference.out', tree = trees)

rule all:
    input:
        simFiles,
        infFiles


# SIMULATION

rule simulate_TiDe:
    input:
        code_dir + '/simulation_barcodes/simulation_TiDe.xml'
    output:
        'TiDe/' + setup + '/simulationOutput/{tree}_{seed}.alignment.nexus'
    shell:
        "java -jar $HOME/beasts2.7.jar -overwrite -seed {wildcards.seed}\
        -D 'tree={wildcards.tree},editRate={editRate},nTargets={nTargets},scarringHeight={scarringHeight},scarringDuration={scarringDuration},outFile={output}'\
        {input}"


rule simulate_Typewriter:
    input:
        code_dir + '/simulation_barcodes/simulation_Typewriter.xml'
    params:
        outDir = 'Typewriter/' + setup + '/simulationOutput'
    output:
        'Typewriter/' + setup + '/simulationOutput/{tree}_{seed}.alignment_1.nexus'
    shell:
        "java -jar $HOME/beasts2.7.jar -overwrite -seed {wildcards.seed}\
        -D 'tree={wildcards.tree},editRate={editRate},nTapes={nTapes},tapeLength={tapeLength},outDir={params.outDir}'\
        {input}"


# INFERENCE

rule infer_TiDe:
    input:
        code_dir + '/inference/inference_TiDe.xml',
        ['TiDe/' + setup + '/simulationOutput/{tree}_' + str(seed) + '.alignment.nexus' for seed in seeds]
    params:
        # Slurm:
        jobName = 'inference_TiDe_' + setup + '_{tree}',
        nArray = str(nSim),
        # LSF:
        #jobName = 'inference_TiDe_' + setup + '_{tree}[1-' + str(nSim) + ']',
        inDir = 'TiDe/' + setup + '/simulationOutput',
        outDir = 'TiDe/' + setup + '/inferenceOutput',
        xmlFile = code_dir + '/inference/inference_TiDe.xml'
    output:
        'TiDe/' + setup + '/inferenceOutput/{tree}_1.inference.out'
    shell:
        #"bash run_inference_TiDe.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}"
        # Slurm: 
        "sbatch --time=120:00:00 --mem-per-cpu=1024 --job-name={params.jobName} --array=1-{params.nArray} --wrap='bash run_inference_TiDe.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}'"
        # LSF: 
        #"bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'bash run_inference_TiDe.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}'"


rule infer_Typewriter:
    input:
        code_dir + '/inference/inference_Typewriter.xml',
        ['Typewriter/' + setup + '/simulationOutput/{tree}_' + str(seed) + '.alignment_1.nexus' for seed in seeds]
    params:
        # Slurm:
        jobName = 'inference_Typewriter_' + setup + '_{tree}',
        nArray = str(nSim),
        # LSF:
        #jobName = 'inference_Typewriter_' + setup + '_{tree}[1-' + str(nSim) + ']',
        inDir = 'Typewriter/' + setup + '/simulationOutput',
        outDir = 'Typewriter/' + setup + '/inferenceOutput',
        xmlFile = code_dir + '/inference/inference_Typewriter.xml'
    output:
        'Typewriter/' + setup + '/inferenceOutput/{tree}_1.inference.out'
    shell:
        #"bash run_inference_Typewriter.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}"
        # Slurm: 
        "sbatch --time=120:00:00 --mem-per-cpu=1024 --job-name={params.jobName} --array=1-{params.nArray} --wrap='bash run_inference_Typewriter.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
        # LSF: 
        #"bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'bash run_inference_Typewriter.sh -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
