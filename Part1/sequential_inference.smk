# -----
# Snakemake pipeline for alignment simulations (for homogeneous cell populations) and inference
# -----

# SETUP

# trees
trees = ['tree_s', 'tree_ss', 'tree_sd', 'tree_sds', 'tree_bd']

# simulation paramaters
nSim = 20
seeds = list(range(1, nSim + 1))

# experimental scenario
method = 'Typewriter'
setup = 'sequential'

# experimental parameters
editRate = 0.45 # [0.01, 0.05, 0.1, 0.45]
if method == 'TiDe':
    nTargets = 20 # [5, 20, 40]
    scarringHeight = 40 # [40, 30, 20]
    scarringDuration = 40 # [40, 20]
elif method == 'Typewriter':
    nTapes = 1 # [1, 5, 20, 40]
    tapeLength = 20 # [2, 5, 10, 20]

# output files
if method == 'TiDe':
    simFiles = expand('TiDe/' + setup + '/simulationOutput/' + '{tree}_{seed}.alignment.nexus', tree = trees, seed = seeds)
    infFiles = expand('TiDe/' + setup + '/inferenceOutput/' + '{tree}_1.inference.out', tree = trees)
elif method == 'Typewriter':
    simFiles = expand('simulationOutput/' + '{tree}_{seed}.alignment_1.nexus', tree = trees, seed = seeds)
    infFiles = expand('inferenceOutput/' + '{tree}_1.inference.out', tree = trees)

rule all:
    input:
        simFiles,
        infFiles


# SIMULATION

rule simulate_TiDe:
    input:
        'simulation_TiDe.xml'
    output:
        'TiDe/' + setup + '/simulationOutput/{tree}_{seed}.alignment.nexus'
    shell:
        "java -jar $HOME/beasts2.7.jar -overwrite -seed {wildcards.seed}\
        -D 'tree={wildcards.tree},editRate={editRate},nTargets={nTargets},scarringHeight={scarringHeight},scarringDuration={scarringDuration},outFile={output}'\
        {input}"


rule simulate_Typewriter:
    input:
        'simulation_Typewriter.xml'
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
        'inference_TiDe.xml',
        ['TiDe/' + setup + '/simulationOutput/{tree}_' + str(seed) + '.alignment.nexus' for seed in seeds]
    params:
        jobName = 'inference_TiDe_' + setup + '_{tree}[1-' + str(nSim) + ']',
        inDir = 'TiDe/' + setup + '/simulationOutput',
        outDir = 'TiDe/' + setup + '/inferenceOutput'
    output:
        'TiDe/' + setup + '/inferenceOutput/{tree}_1.inference.out'
    shell:
        #"sh run_inference_TiDe.sh -t {wildcards.tree} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}"
        "bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'sh run_inference_TiDe.sh -t {wildcards.tree} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}'"


rule infer_Typewriter:
    input:
        'inference_Typewriter_sequential.xml',
        ['simulationOutput/{tree}_' + str(seed) + '.alignment_1.nexus' for seed in seeds]
    params:
        jobName = 'inference_Typewriter_' + setup + '_{tree}',
        nArray = str(nSim),
        inDir = 'simulationOutput',
        outDir = 'inferenceOutput',
    output:
        'inferenceOutput/{tree}_1.inference.out'
    shell:
        #"sh run_inference_Typewriter.sh -t {wildcards.tree} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}"
        "sbatch --time=120:00:00 --mem-per-cpu=1024 --job-name={params.jobName} --array=1-{params.nArray} --wrap='sh run_inference_Typewriter.sh -t {wildcards.tree} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
        #"bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'sh run_inference_Typewriter.sh -t {wildcards.tree} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
