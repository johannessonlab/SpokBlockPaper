# -*- snakemake -*-

from Bio import SeqIO
from Bio.Alphabet import generic_dna

### A pipeline to annotate Podospora complex genomes
#############################################################################
# Pipeline to produce de novo repeat libraries with RepeatModeler
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018/11/22 - 2019/11/14
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
renamer = config["renamer"]

MINSIZE = config["MINSIZE"] # Minimum size of a contig to get annotated
MINSIZEkb = int(MINSIZE/1000)

AllSamples = config["samplesID"]

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: fastareduce, builddb, renameconsensi #, referencegenome
# ----------

rule all:
	input:
		# RM results
		expand("results/{sample}/{sample}_RM.lib", sample = AllSamples),

# ---------------------------------

rule fastareduce:
	""" Remove small contigs from the input assemblies """
	input:
		genome = "data/genomes/{sample}.fa"
	output:
		simplifiedfasta = "data/reduced/{sample}_{MINSIZEkb}kp.fa"
	run:
		filteredseqs = [] # Setup an empty list
		for seq_record in SeqIO.parse(input.genome, "fasta"):
			if len(seq_record) >= MINSIZE:
				filteredseqs.append(seq_record)
		# Write to a new file
		SeqIO.write(filteredseqs, output.simplifiedfasta, "fasta")

rule builddb:
	""" Prepare BLAST database """
	input:
		genome = "data/reduced/{sample}_" + str(MINSIZEkb) + "kp.fa"
	output:
		"RepeatModeler/{sample}/{sample}_db/{sample}.nhr",
		"RepeatModeler/{sample}/{sample}_db/{sample}.nin",
		"RepeatModeler/{sample}/{sample}_db/{sample}.nnd",
		"RepeatModeler/{sample}/{sample}_db/{sample}.nni",
		"RepeatModeler/{sample}/{sample}_db/{sample}.nog",
		"RepeatModeler/{sample}/{sample}_db/{sample}.nsq",
		"RepeatModeler/{sample}/{sample}_db/{sample}.translation",	
	shell:
		"""
		# Format FASTA files for use with RepeatModeler
		BuildDatabase -name RepeatModeler/{wildcards.sample}/{wildcards.sample}_db/{wildcards.sample} -engine ncbi {input.genome}
		"""

rule repeatmodeler:
	""" Create new TE libraries per genome """
	input:
		genome = "data/reduced/{sample}_" + str(MINSIZEkb) + "kp.fa",
		database = "RepeatModeler/{sample}/{sample}_db/{sample}.nhr",
	output:
		"RepeatModeler/{sample}/{sample}_consensi.fa",
		"RepeatModeler/{sample}/{sample}_consensi.fa.classified",
		"RepeatModeler/{sample}/{sample}_consensi.fa.masked",
	params:
		time = "5-00:00:00",
		threads = 8,
	shell:
		"""
		# RepeatModeler starts the folder in the working directory so move there
		cd RepeatModeler/{wildcards.sample}

		# Run RepeatModeler
		RepeatModeler -pa {params.threads} -database {wildcards.sample}_db/{wildcards.sample}
		# RepeatModeler -pa {params.threads} -database {wildcards.sample}_db/{wildcards.sample} -srand 1234 # a higher version of RepeatModeler that doesn't work for me so far

		# Rename output
		cp RM_*/consensi.fa {wildcards.sample}_consensi.fa
		cp RM_*/consensi.fa.classified {wildcards.sample}_consensi.fa.classified
		cp RM_*/consensi.fa.masked {wildcards.sample}_consensi.fa.masked 
		"""

rule renameconsensi:
	input:
		rawlib = "RepeatModeler/{sample}/{sample}_consensi.fa.classified",
	output:
		"results/{sample}/{sample}_RM.lib"
	params:
		renamer = renamer
	shell:
		""" 
		perl {params.renamer} {input.rawlib} {wildcards.sample} {output}
		"""


