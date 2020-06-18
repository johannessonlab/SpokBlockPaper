# TEManualCuration: A pipeline for recovering alignments for manual curation of TE libraries of Podospora spp

The genomes to be annotated must be in a folder `data/genomes/`. So from the working path, a given sample looks like `data/genomes/{sample}.fa`, where `sample` comes from the configuration file.

It also needs a folder `scripts` with the script `repeatModelerPipeline4.pl` available [here](https://github.com/genomicrocosm/physaliaTEcourse/blob/master/Practical3_Manual_curation/repeatModelerPipeline4.pl).

So the necessary files are (that can be modified in the configuration file):

    data:
    genomes/

    scripts:
    repeatModelerPipeline4.pl

## Building the environment

For this pipeline I relied partially on locally installed software in [Uppmax](https://www.uppmax.uu.se/) and partially in conda environments.

First, one can start by updating conda.

    $ conda update -n base conda

I made a small yaml file for an environment containing t-coffee created on the fly for the corresponding rule.

    $ cat envs/tcoffee.yaml
    channels:
      - bioconda
      - defaults
      - conda-forge
    dependencies:
      - t-coffee=12.00.7fb08c2=hfc679d8_1

## Run pipeline in Uppmax

Get into the folder:

    $ cd TEManualCuration

For this I'm relying on Uppmax libraries:

    $ module load bioinfo-tools snakemake/5.2.3 

First, to get an idea of how the pipeline we can make a rulegraph:

    $ snakemake --snakefile TEManualCuration.smk --configfile TEManualCuration_config.yml --rulegraph | dot -Tpng > rulegraph.png

Test the pipeline:

    $ snakemake --snakefile TEManualCuration.smk --configfile TEManualCuration_config.yml -pn

Run it:
    
    $ screen -R TEManualCuration
    $ module load bioinfo-tools snakemake/5.2.3 perl_modules/5.24.1 blast/2.7.1+ MAFFT/7.407
    $ unset $MAFFT_BINARIES
    $ snakemake --snakefile TEManualCuration.smk --configfile TEManualCuration_config.yml -p --cluster "sbatch -A snic2017-1-567 -p core -n {params.threads} -t {params.time} --mail-user sandra.ament@evobio.eu --mail-type=ALL" -j 15 --keep-going --use-conda &> TEManualCuration.log &
    [1] 3390

Notice the `--use-conda` argument!!

A few samples (PaWa100p, PaWa58m, and CBS112042p) run out of memory sometimes. They don't make Snakemake fail, but some of the final alignments will be empty. Check the end of the slurm files for a warning. For those, I asked for the "fat" nodes that have 256 GB of memory, instead of the "thin" nodes that have 128 GB. However the queue is MUCH slower.

    $ rm repeatModelerPipeline/PaWa100p/PaWa100p.dummy repeatModelerPipeline/PaWa58m/PaWa58m.dummy repeatModelerPipeline/CBS112042p/CBS112042p.dummy
    $ snakemake --snakefile TEManualCuration.smk --configfile TEManualCuration_config.yml -p --cluster "sbatch -A snic2017-1-567 -p core -n {params.threads} -t {params.time} -C fat --mail-user sandra.ament@evobio.eu --mail-type=ALL" -j 15 --keep-going --use-conda &> TEManualCuration.log2 &
    [1] 19899

So all the files needed for the manual curation are in the path `TEManualCuration/repeatModelerPipeline`.


