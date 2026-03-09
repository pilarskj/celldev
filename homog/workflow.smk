# -----
# Snakemake pipeline for alignment simulations (for homogeneous cell populations) and inference
# -----

# SETUP
import os

# paths (adapt!)
codeDir = '/cluster/home/jpilarski/celldev'
dataDir = '/cluster/scratch/jpilarski/celldev_data'
os.chdir(codeDir)
treeDir = dataDir + '/Trees'

# trees
trees = ['tree_s', 'tree_ss', 'tree_sd', 'tree_sds', 'tree_bd']

# simulation paramaters
nSim = 20
seeds = list(range(1, nSim + 1))

# experimental scenario
method = 'TiDe'
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
    simFiles = expand(dataDir + '/TiDe/' + setup + '/simulationOutput/' + '{tree}_{seed}.alignment.nexus', tree = trees, seed = seeds)
    infFiles = expand(dataDir + '/TiDe/' + setup + '/inferenceOutput/' + '{tree}_1.inference.out', tree = trees)
elif method == 'Typewriter':
    simFiles = expand(dataDir + '/Typewriter/' + setup + '/simulationOutput/' + '{tree}_{seed}.alignment_1.nexus', tree = trees, seed = seeds)
    infFiles = expand(dataDir + '/Typewriter/' + setup + '/inferenceOutput/' + '{tree}_1.inference.out', tree = trees)

rule all:
    input:
        simFiles,
        infFiles


# SIMULATION

rule simulate_TiDe:
    input:
        codeDir + '/homog/simulation_barcodes/simulation_TiDe.xml'
    output:
        dataDir + '/TiDe/' + setup + '/simulationOutput/{tree}_{seed}.alignment.nexus'
    shell:
        # adapt to running environment
        "java -Dglass.platform=Monocle -Dmonocle.platform=Headless --module-path=$HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
        -jar software/simbundle.jar \
        -version_file software/beast2_version.xml \
        -version_file software/feast_version.xml \
        -version_file software/tidetree_version.xml \
        -version_file software/sciphy_version.xml \
        -overwrite \
        -seed {wildcards.seed} \
        -D 'treeDir={treeDir},tree={wildcards.tree},editRate={editRate},nTargets={nTargets},scarringHeight={scarringHeight},scarringDuration={scarringDuration},outFile={output}' \
        {input} >> {dataDir}/outs/simulation_TiDe_{setup}.{wildcards.tree}_{wildcards.seed}.out"


rule simulate_Typewriter:
    input:
        codeDir + '/homog/simulation_barcodes/simulation_Typewriter.xml'
    params:
        outDir = dataDir + '/Typewriter/' + setup + '/simulationOutput'
    output:
        dataDir + '/Typewriter/' + setup + '/simulationOutput/{tree}_{seed}.alignment_1.nexus'
    shell:
        "java -Dglass.platform=Monocle -Dmonocle.platform=Headless --module-path=$HOME/javafx-sdk-17.0.6-linux-monocle/lib --add-modules=javafx.base,javafx.fxml \
        -jar software/simbundle.jar \
        -version_file software/beast2_version.xml \
        -version_file software/feast_version.xml \
        -version_file software/tidetree_version.xml \
        -version_file software/sciphy_version.xml \
        -overwrite \
        -seed {wildcards.seed} \
        -D 'treeDir={treeDir},tree={wildcards.tree},editRate={editRate},nTapes={nTapes},tapeLength={tapeLength},outDir={params.outDir}' \
        {input} >> {dataDir}/outs/simulation_Typewriter_{setup}.{wildcards.tree}_{wildcards.seed}.out"


# INFERENCE

rule infer_TiDe:
    input:
        codeDir + '/homog/inference/inference_TiDe.xml',
        [dataDir + '/TiDe/' + setup + '/simulationOutput/{tree}_' + str(seed) + '.alignment.nexus' for seed in seeds]
    params:
        # Slurm:
        jobName = 'inference_TiDe_' + setup + '_{tree}',
        nArray = str(nSim),
        # LSF:
        #jobName = 'inference_TiDe_' + setup + '_{tree}[1-' + str(nSim) + ']',
        inDir = dataDir + '/TiDe/' + setup + '/simulationOutput',
        outDir = dataDir + '/TiDe/' + setup + '/inferenceOutput',
        xmlFile = codeDir + '/homog/inference/inference_TiDe.xml',
        shFile = codeDir + '/homog/inference/run_inference_TiDe.sh'
    output:
        dataDir + '/TiDe/' + setup + '/inferenceOutput/{tree}_1.inference.out'
    shell:
        #"bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}"
        # Slurm: 
        "sbatch --time=120:00:00 --mem-per-cpu=2G --job-name={params.jobName} --array=1-{params.nArray} --output={dataDir}/outs/%x_%a.out --wrap='bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}'"
        # LSF: 
        #"bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -h {scarringHeight} -d {scarringDuration}'"


rule infer_Typewriter:
    input:
        codeDir + '/homog/inference/inference_Typewriter.xml',
        [dataDir + '/Typewriter/' + setup + '/simulationOutput/{tree}_' + str(seed) + '.alignment_1.nexus' for seed in seeds]
    params:
        # Slurm:
        jobName = 'inference_Typewriter_' + setup + '_{tree}',
        nArray = str(nSim),
        # LSF:
        #jobName = 'inference_Typewriter_' + setup + '_{tree}[1-' + str(nSim) + ']',
        inDir = dataDir + '/Typewriter/' + setup + '/simulationOutput',
        outDir = dataDir + '/Typewriter/' + setup + '/inferenceOutput',
        xmlFile = codeDir + '/homog/inference/inference_Typewriter.xml',
        shFile = codeDir + '/homog/inference/run_inference_Typewriter.sh'
    output:
        dataDir + '/Typewriter/' + setup + '/inferenceOutput/{tree}_1.inference.out'
    shell:
        #"bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}"
        # Slurm: 
        "sbatch --time=120:00:00 --mem-per-cpu=2G --job-name={params.jobName} --array=1-{params.nArray} --output={dataDir}/outs/%x_%a.out --wrap='bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
        # LSF: 
        #"bsub -W 120:00 -R 'rusage[mem=1024]' -J {params.jobName} 'bash {params.shFile} -t {wildcards.tree} -x {params.xmlFile} -i {params.inDir} -o {params.outDir} -n {nTapes} -l {tapeLength}'"
