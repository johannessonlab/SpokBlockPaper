# -*- snakemake -*-

from Bio.Seq import Seq
from glob import glob

### PoJellyfish: Producing the distribution of kmers of a given size in *Podospora* genomes
#############################################################################
# Pipeline to calculate the distribution of k-mers (substring) of length six
# in the genome(s) of Podospora species using Jellyfish
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/09/21
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1: GitHub version

# -------------------------------------------------
samples = config["SampleIDs"]
path2assemblies = config["path2assemblies"]
kmers = config["kmers"]
KmerDist = config["KmerDist"]
# -------------------------------------------------

## Make a dictionary of the samples (key) and their assemblies (values)
assemblies = [(glob(path2assemblies + "/{sample}*.fa*".format(sample=sample))) for sample in samples]
ASSEMBLIESDIC = dict(zip(samples, assemblies))

# print(ASSEMBLIESDIC)

lenkmers = list(set([len(kmer) for kmer in kmers])) # Get lengths of relevant kmers

rule all:
	input:
		expand("results/{sample}_kmer_hist_k{klen}.pdf", sample = samples, klen = set([len(kmer) for kmer in kmers]))
		
rule jellyfish:
	""" Get the kmers of size relevant to the kmers """
	input:
		assembly = lambda wildcards: ASSEMBLIESDIC[wildcards.sample]
	output:
		counts = "jellyfish/{sample}/{sample}_mer_counts_k{klen}.jf"
	threads: 4
	shell:
		"jellyfish count -m {wildcards.klen} -s 100M -t {threads} -C {input.assembly} -o {output.counts}"
		# "echo jellyfish count -m {wildcards.klen} -s 100M -t {threads} -C {input.assembly} -o {output.counts};"

rule dump_kmers:
	""" Get table of kmers """
	input:
		"jellyfish/{sample}/{sample}_mer_counts_k{klen}.jf"
	output:
		temp("jellyfish/{sample}/{sample}_mer_counts_k{klen}.tbl_raw")
	shell:
		"jellyfish dump {input} -c > {output}"
		# -c table output instead of fasta


def revcom_kmers(kmers_list):
	""" Get Reverse complement of kmers """
	allkmers = []

	for kmer in kmers_list:
		kseq = Seq(kmer) # Make it a sequence
		rvkmer = kseq.reverse_complement()
		allkmers.append(str(rvkmer))

	return(allkmers)


rule fixkmerlist:
	""" Keep the kmers as requested and not their reverse complement name in the list """
	input:
		tbl = "jellyfish/{sample}/{sample}_mer_counts_k{klen}.tbl_raw"
	output:
		tbl = "jellyfish/{sample}/{sample}_mer_counts_k{klen}.tbl"
	run:
		ofile = open(output.tbl, 'w')
		tabopen = open(input.tbl, 'r')

		rckmers = revcom_kmers(kmers) # The reverse complement of the kmers of interest

		for line in tabopen:
			tab = line.rstrip("\n").split(" ")
			kmer = tab[0]
		
			if kmer in rckmers:
				dereverse = revcom_kmers([kmer])[0] # reverse complement it to be back to the original
				newline = f"{dereverse}\t{tab[1]}\n" # Replace the Jellyfish kmer for the specified one
				ofile.write(newline)	
			else:
				ofile.write(line)	

rule getTSDcount:
	""" Get the TSD count for each TSD """
	input:
		# "jellyfish/{sample}/{sample}_mer_counts_k{klen}.tbl",
		expand("jellyfish/{sample}/{sample}_mer_counts_k{{klen}}.jf", sample = samples)
	output:
		"counts/Counts_k{klen}_{kmer}.txt"	
	run:
		for n in range(0, len(input)):
			sample = samples[n]
			file = input[n]

			shell("printf '{sample} '>> {output}")
			shell("jellyfish query {file} {wildcards.kmer} >> {output}")

rule mergecounts:
	""" Put all counts in a long format table for R """
	input:
		expand("counts/Counts_k{klen}_{kmer}.txt", zip, kmer = kmers, klen = [len(kmer) for kmer in kmers]),
	output:
		tbl = "counts/Counts_per_sample.txt"
	run: 
		ofile = open(output.tbl, 'w')
		## Jellyfish reports the kmer or its reverse complement, whatever comes first, but I do want the specified kmer
		for n in range(0, len(input)):
			kmer = kmers[n] # Get the kmer name
			file = input[n]

			tabopen = open(file, 'r')
			for line in tabopen:
				tab = line.rstrip("\n").split(" ")
				newline = f"{tab[0]}\t{kmer}\t{tab[2]}\n" # Replace the Jellyfish kmer for the specified one
				ofile.write(newline)

rule kmerdistribution:
	""" How rare is this kmer compared to others? """
	input:
		counts = "jellyfish/{sample}/{sample}_mer_counts_k{klen}.tbl",
		script = KmerDist
	output:
		hist = "results/{sample}_kmer_hist_k{klen}.pdf"
	params:
		kmers = kmers
	script:
		KmerDist
