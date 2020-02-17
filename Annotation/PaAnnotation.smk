# -*- snakemake -*-

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os.path

### A pipeline to annotate Podospora complex genomes: Spok block paper
#############################################################################

# Second version of the annotation. The pipeline is changed slightly for the
# naming, but now hopefully it's a bit more efficient. And I attempt to add a
# few gene products. 
# I decided to drop the annotation transfer of Podan with MAKER alone because
# it's really bad and not useful at all.

# https://bioinformaticsworkbook.org/dataAnalysis/GenomeAnnotation/Intro_To_Maker.html

# TODO:
# - change the evm.ctl to give more weight to genemark
# - Check that MAKERcpus works (not tested yet)

#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2018/11/22
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 2

# -------------------------------------------------

annotationversion = "{0:.2f}".format(config["annotationversion"])

## Samples
sampleslongreads = config["sampleslongreads"]
samplesillu = config["samplesillu"]

## RNAseq data
samplesRNA = config["samplesRNA"]
parentals = config["parentals"]

## TEs
TElib = config["PodoTElib"]

## Scripts
gff2gff3 = config["gff2gff3"]
runTransDecoder = config["runTransDecoder"]
gff3tofasta = config["gff3tofasta"]
GffgenesIDFix = config["GffgenesIDFix"]
GFFnumerator = config["GFFnumerator"]
gff3addproduct = config["gff3addproduct"]

## Samples called
AllSamples = sampleslongreads + samplesillu

## Training files
snapHMM = config["snapHMM"]
GeneMarkMod = config["GeneMarkMod"]

## Podan2
refgenesfas = config["Podan2genesfas"]

## Known gene products
geneproducts = config["geneproducts"]

## Values sample specific
periden = config["periden"]
sppcode = config["sppcode"]

def getperind(wildcards):
	# Get index of this sample 
	sampindex = AllSamples.index(wildcards.sample)
	thissampleperiden = periden[sampindex]
	return(thissampleperiden)

def getsppcode(wildcards):
	# Get index of this sample 
	sampindex = AllSamples.index(wildcards.sample)
	thissamplesppcode = sppcode[sampindex]
	return(thissamplesppcode)

# The same thing but different names
def typesample(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"temp/genomes/{wildcards.sample}_tempnames.fa"
	else:
		genome = f"data/genomes/{wildcards.sample}.fa"
	return(genome)

def typesample_rm(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"RepeatMasker/{wildcards.sample}/{wildcards.sample}.repeatmasker.illu.gff"
	else:
		genome = f"RepeatMasker/{wildcards.sample}/{wildcards.sample}.fa.out.gff"
	return(genome)

def typesample_maker(wildcards):
	""" Function to make a distinction between short and long-read assemblies """
	if wildcards.sample in samplesillu:
		genome = f"MAKER/{wildcards.sample}/MAKER_output/{wildcards.sample}.all.illu.gff"
	else:
		genome = f"MAKER/{wildcards.sample}/MAKER_output/{wildcards.sample}.all.gff"
	return(genome)

# ------
MINSIZE = config["MINSIZE"] # Minimum size of a contig to get annotated
MAKERcpus = 2 # Jobs are not very efficient, so it's not worth it to give that many threads.

# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: gtf2gff, gff2gff3, makerconfig, collectmakerresults, collectRMresults, getfastamaker, addgeneproducts, sppcodenicegff, newnames, renameMAKERresults, renameRMresults
# ----------

rule all:
	input:
		### Results
		# #### RepeatMasker	
		expand("results/{sample}/{sample}.repeatmasker.gff", sample = AllSamples ), # RepeatMasker

		### RNAseq
		# Stats
		expand("STAR/{genome}/{rnasample}-to-{genome}_flagstat.txt", zip, rnasample = samplesRNA, genome = parentals),
		# Index for IGV
		expand("STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam.bai", zip, rnasample = samplesRNA, genome = parentals),
		# Transcripts
		expand("Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff3", zip, rnasample = samplesRNA, genome = parentals),
		# TransDecoder
		expand("TransDecoder/{genome}/{rnasample}/{rnasample}.genome.transdecoder.gff3", zip, rnasample = samplesRNA, genome = parentals),

		#### MAKER
		expand("results/{sample}/{sample}.all.gff3", sample = AllSamples), # MAKER
		expand("results/{sample}/{sample}.nice-{version}.gff3", sample = AllSamples, version = annotationversion), # Renamed, final models

# ---------------------------
### RNAseq
# ---------------------------
rule starindex:
	""" Make an index of the genome for STAR """
	input:
		genome = "data/genomes/{sample}.fa",
	output:
		"STAR/GenomeStarIndex/{sample}_GenomeDir/SAindex"
	params:
		indexdir = "STAR/GenomeStarIndex/{sample}_GenomeDir",
		time = "30:00",
		threads = 1,
		ram = int(3 * 6.8 * 1000000000) # A Rackham node contains 128 GB of RAM and 20 compute cores (each core gets at most 6.8 GB).
	shell:
		"""
		mkdir -p temp/STAR/

		STAR --runMode genomeGenerate --genomeDir {params.indexdir} --genomeFastaFiles {input.genome} --runThreadN {params.threads} --limitGenomeGenerateRAM {params.ram} --genomeLoad NoSharedMemory --outTmpDir temp/STAR/{wildcards.sample} # --genomeSAindexNbases 3
		
		# genomeLoad=NoSharedMemory, shared memory is not used. This option is recommended if the shared memory is not configured properly on your server.
		# --genomeSAIndexNbases 4 or 5 for small genomes (formula min(14, log2(GenomeLength)/2 - 1)) # log2(37000000)/2 -1 = 11.57051
		"""

rule star:
	""" Map the RNAseq reads to a genome using STAR """
	input:
		genome = "data/genomes/{genome}.fa",
		index = "STAR/GenomeStarIndex/{genome}_GenomeDir/SAindex",
		read1 = "data/rnaseq/{rnasample}_postQC.1.fq.gz",
		read2 = "data/rnaseq/{rnasample}_postQC.2.fq.gz",
	output:
		output = "STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	params:
		indexdir = "STAR/GenomeStarIndex/{genome}_GenomeDir",
		time = "10:00:00",
		threads = 3,
	shell:
		"""
		STAR --genomeDir {params.indexdir} --readFilesIn {input.read1} {input.read2} \
		--runThreadN {params.threads} \
		--alignIntronMax 1000 \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR/{wildcards.genome}/{wildcards.rnasample}-to-{wildcards.genome}_ \
		--outSAMstrandField intronMotif \
		--outSAMtype BAM SortedByCoordinate

		# --outSAMattributes NH HI AS nM XS 

		# outSAMattributes To make it compatible with SAMtools
		# outSAMstrandField intronMotif to make it compatible with Cufflinks/Cuffdiff

		# For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments
		# with XS strand attribute, which STAR will generate with --outSAMstrandField
		# intronMotif option. As required, the XS strand attribute will be generated for
		# all alignments that contain splice junctions. The spliced alignments that have
		# undefined strand (i.e. containing only non-canonical unannotated junctions)
		# will be suppressed.

		# If you have stranded RNA-seq data, you do not need to use any specific STAR
		# options. Instead, you need to run Cufflinks with the library option --library-
		# type options. For example, cufflinks ... --library-type fr-firststrand should
		# be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq.
		# This option has to be used only for Cufflinks runs and not for STAR runs.

		# --outSAMtype BAM SortedByCoordinate output sorted by coordinate Aligned.sortedByCoord.out.bam le, similar to samtools sort command.
		"""

rule bamstats:
	""" Get some statistics of the STAR BAM file """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{genome}/{rnasample}-to-{genome}_flagstat.txt"
	params:
		time = "30:00",
		threads = 1,
	shell:
		"samtools flagstat {input} > {output}"

rule bamindex:
	""" Make an index for the BAM file to read it in IGV """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam.bai"
	params:
		time = "30:00",
		threads = 1,
	shell:
		"samtools index {input}"

rule cufflinks:
	""" Produce transcript models with Cufflinks """
	input:
		"STAR/{genome}/{rnasample}-to-{genome}_Aligned.sortedByCoord.out.bam",
	output:
		"Cufflinks/{genome}/{rnasample}/transcripts.gtf"
	params:
		time = "05:00:00",
		threads = 10,
	shell:
		"""
		# Stranded RNAseq 
		cufflinks {input} --library-type fr-firststrand -L {wildcards.rnasample}_CUFF --output-dir Cufflinks/{wildcards.genome}/{wildcards.rnasample} --num-threads {params.threads} 2>&1  # Redirect the stderr to stdout
		# fr-firststrand --> dUTP libraries

		# mv Cufflinks/{wildcards.genome}/transcripts.gtf Cufflinks/{wildcards.genome}/{wildcards.rnasample}-to-{wildcards.genome}_transcripts.gtf 
		"""

rule gtf2gff:
	""" Transform the gtf Cufflinks file into gff3 """
	input:
		"Cufflinks/{genome}/{rnasample}/transcripts.gtf"
	output:
		temp("Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff")
	shell:
		"""
		gffread -E {input} -o- > {output} # add | tail -n +3 to take out the headers
		"""

rule gff2gff3:
	""" Transform the gtf Cufflinks file into gff3 """
	input:
		"Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff"
	output:
		"Cufflinks/{genome}/{rnasample}/{rnasample}_transcripts.gff3"
	params:
		gff2gff3 = gff2gff3
	shell:
		"""
		{params.gff2gff3} {input} --prefix {wildcards.rnasample}_CUFF > {output}
		"""

rule transdecoder:
	""" Run TransDecoder to get ORFs from the transcripts """
	input:
		transcripts = "Cufflinks/{genome}/{rnasample}/transcripts.gtf",
		genome = "data/genomes/{genome}.fa",
	output:
		"TransDecoder/{genome}/{rnasample}/{rnasample}.genome.transdecoder.gff3",
		"TransDecoder/{genome}/{rnasample}/{rnasample}.fa" # Transcripts
	params:
		runTransDecoder = runTransDecoder,
		time = "03:00:00",
		threads = 1,
	shell:
		"""
		cd TransDecoder/{wildcards.genome}/{wildcards.rnasample}
		{params.runTransDecoder} ../../../{input.transcripts} ../../../{input.genome} {wildcards.rnasample}
		"""

# ---------------------------
### Repeatmasking
# ---------------------------

# Illumina assemblies from SPAdes have names that are too long.
rule newnames:
	""" Make a new fasta temporary names for the scaffolds for Illumina assemblies """
	input:
		genome = "data/genomes/{sample}.fa"
	output:
		tempnames = "temp/genomes/{sample}_tempnames.txt", # To keep the original names somewhere
		newgenome = "temp/genomes/{sample}_tempnames.fa"
	# wildcard_constraints: # to avoid ambiguity with tempnames
	# 	sample="\d+"
	run:
		count = 1
		dictnames = {} # Setup a dictionary
		output_handle = open(output.newgenome, "w")

		for seq_record in SeqIO.parse(input.genome, "fasta"):
			if seq_record.id in dictnames.keys(): # Just in case
				print(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
				raise ValueError(f"ERROR: The sequence {seq_record.id} is more than once in the fasta file {input.genome}")
			else:
				newname = "seq{0:06d}".format(count)
				dictnames[seq_record.id] = newname # record the original name somewhere
				count += 1 # Increase the count for the next sequence
				seq_record.id = newname
				seq_record.description = '' #Annoying extra info
				SeqIO.write(seq_record, output_handle, "fasta")

		# Keep the names somewhere
		with open(output.tempnames, 'w') as result:
			for key, value in dictnames.items():
				result.write(f"{key}\t{value}\n")

rule repeatmasker:
	""" Run RepeatMasker on the input genomes """
	input:
		genome = typesample
	output:
		"RepeatMasker/{sample}/{sample}.fa.out.gff"
	params:
		TElib = TElib, # Custom library
		time = "10:00:00",
		threads = 10,
	run:
		shell("RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker/{wildcards.sample} {input.genome}")
		if os.path.exists("RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff"):
			shell("mv RepeatMasker/{wildcards.sample}/{wildcards.sample}_tempnames.fa.out.gff {output}") # Rename to deal with tempname

	# shell:
	# 	"RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {params.TElib} -dir RepeatMasker/{wildcards.sample} {input.genome}; "
	# 	"mv RepeatMasker/{wildcards.sample}/{wildcards.sample}*.fa.out.gff {output}" # Rename to deal with tempname

# ---------------------------
### MAKER
# ---------------------------

rule makerconfig:
	""" Produce MAKER configuration files for the sample """
	# This depends on Uppmax module
	# module load bioinfo-tools maker/3.01.2-beta
	# See http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
	input:
		genome = typesample, # If it's illumina, you need the simple names

		# EST evidence
		transcWa63 = "TransDecoder/PaWa63p/PaWa63m_RNA/PaWa63m_RNA.fa",
		transcWa58BC = "TransDecoder/PaWa58m/PaWa58BCvsS/PaWa58BCvsS.fa",
		transcWa131 = "TransDecoder/PcWa139m/PcWa131m_RNA/PcWa131m_RNA.fa",

		# Other evidence
		podan2aa = "data/Podan2/Podan2_AssemblyScaffoldsmtGenesEd_aa.fas",
		comata = "data/OtherSpp/Pcomata_aa.fas",
		curated = "data/OtherSpp/curated_aa.fas",

	output:
		"MAKER/{sample}/{sample}_bopts.ctl",
		"MAKER/{sample}/{sample}_exe.ctl",
		"MAKER/{sample}/{sample}_opts.ctl",
	params:
		TElib = TElib,
		snapHMM = snapHMM,
		GeneMarkMod = GeneMarkMod,
		repeatmaskergff = "", # I wanted this to be used, but I suspect the gff won't be good
		cpus = MAKERcpus,
	shell:
		"""	
		cd MAKER/{wildcards.sample}
	
		# Produce an empty configuration file
		maker -CTL

		# Rename them to match the sample
		mv maker_bopts.ctl {wildcards.sample}_bopts.ctl
		mv maker_exe.ctl {wildcards.sample}_exe.ctl
		mv maker_opts.ctl {wildcards.sample}_opts.ctl
		mv maker_evm.ctl {wildcards.sample}_evm.ctl

		### Fill the config files with the sample info
		# Set the current genome
		sed -i "s|^genome=|genome=../../{input.genome}|g" {wildcards.sample}_opts.ctl # Add that to the opts file
		
		# #-----Set evidence
		sed -i "s|^est=|est=../../{input.transcWa63},../../{input.transcWa58BC},../../{input.transcWa131}|g" {wildcards.sample}_opts.ctl 	# The RNAseq STAR+Cufflinks transcripts
		sed -i "s|^protein=|protein=../../{input.podan2aa},../../{input.comata},../../{input.curated}|g" {wildcards.sample}_opts.ctl # Add that to the opts file

		# #-----Repeat Masking
		sed -i "s|^rmlib=|rmlib={params.TElib}|g" {wildcards.sample}_opts.ctl  #pre-identified repeat elements in a fasta file	
		# sed -i "s|^rm_gff=|rm_gff=../../{params.repeatmaskergff}|g" {wildcards.sample}_opts.ctl  #pre-identified repeat elements from an external GFF3 file

		# Notice I'm leaving it with softmask=1

		# #-----Gene Prediction
		# SNAP HMM file
		sed -i "s|^snaphmm=|snaphmm={params.snapHMM}|g" {wildcards.sample}_opts.ctl  
		# GeneMark HMM file
		sed -i "s|^gmhmm=|gmhmm={params.GeneMarkMod}|g" {wildcards.sample}_opts.ctl  

		# #infer predictions from protein homology, 1 = yes, 0 = no # Otherwise many small genes from Podan2 are lost, even tho some of them might not be real in fact (but some are for sure base on RNAseq)
		sed -i "s|^protein2genome=0|protein2genome=1|g" {wildcards.sample}_opts.ctl  
		# find tRNAs with tRNAscan, 1 = yes, 0 = no
		sed -i "s|^trna=0|trna=1|g" {wildcards.sample}_opts.ctl  
		# also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
		sed -i "s|^unmask=0|unmask=1|g" {wildcards.sample}_opts.ctl  
				
		# #-----External Application Behavior Options
		sed -i "s|^cpus=1|cpus={params.cpus}|g" {wildcards.sample}_opts.ctl  #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)
		
		# #-----MAKER Behavior Options
		# skip genome contigs below this length (under 10kb are often useless)
		sed -i "s|^min_contig=1|min_contig=10000|g" {wildcards.sample}_opts.ctl 
		# extra steps to force start and stop codons, 1 = yes, 0 = no
		sed -i "s|^always_complete=0|always_complete=1|g" {wildcards.sample}_opts.ctl # Jesper set this to 1, but I'm actually not sure.
		# length for the splitting of hits (expected max intron size for evidence alignments)
		sed -i "s|^split_hit=10000|split_hit=1000|g" {wildcards.sample}_opts.ctl # Jesper set this to 1000
		# #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
		sed -i "s|^single_exon=0|single_exon=1|g" {wildcards.sample}_opts.ctl # It's dangerous but I hope to recover genes/things expressed that are not in Podan2
		# #min length required for single exon ESTs if 'single_exon is enabled'
		sed -i "s|^single_length=250|single_length=500|g" {wildcards.sample}_opts.ctl # to compensate a bit
		
		# #limits use of ESTs in annotation to avoid fusion genes (This results in the loss of three prime UTR on the five prime gene and 
		# # loss of five prime UTR on the three prime gene, but it is better than a merged gene.)
		sed -i "s|^correct_est_fusion=0|correct_est_fusion=1|g" {wildcards.sample}_opts.ctl
		
		## Danger
		# #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
		sed -i "s|^clean_up=0|clean_up=1|g" {wildcards.sample}_opts.ctl
		
		# This option will help save disk space by deleting individual results
		# files (such as blast, exonerate, and gene predictor outputs) once
		# they are no longer needed. If you have the disk space it is usually
		# best to keep this set to 0. Having those files around will make
		# rerunning MAKER much faster if necessary.

		"""

rule maker:
	""" Run MAKER """
	# This depends on Uppmax module
	# module load bioinfo-tools maker/3.01.2-beta
	# See http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
	input:
		"MAKER/{sample}/{sample}_bopts.ctl",
		"MAKER/{sample}/{sample}_exe.ctl",
		"MAKER/{sample}/{sample}_opts.ctl",
	output:
		"MAKER/{sample}/{sample}.dummy",
	params:
		time = "10-00:00:00",
		threads = MAKERcpus,
	shell:
		"""	
		cd MAKER/{wildcards.sample}
		maker -fix_nucleotides {wildcards.sample}_opts.ctl {wildcards.sample}_bopts.ctl {wildcards.sample}_exe.ctl {wildcards.sample}_opts.ctl && touch {wildcards.sample}.dummy
		"""

rule makermerge:
	""" Recover output of MAKER """
	input:
		"MAKER/{sample}/{sample}.dummy",
	output:
		"MAKER/{sample}/MAKER_output/{sample}.all.gff",
	params:
		time = "15:00",
		threads = 1,
	run:
		sample = wildcards.sample
		if wildcards.sample in samplesillu: # To deal with illumina samples
			sample += "_tempnames"

		shell("""
			cd MAKER/{wildcards.sample}/MAKER_output
			fasta_merge -d ../{sample}.maker.output/{sample}_master_datastore_index.log
			gff3_merge -d ../{sample}.maker.output/{sample}_master_datastore_index.log
			""")

		if wildcards.sample in samplesillu: # Rename the output so it works with the rest of the rulegraph
			shell("mv MAKER/{wildcards.sample}/MAKER_output/{sample}.all.gff {output}")
	# shell:
	# 	"""
	# 	cd MAKER/{wildcards.sample}/MAKER_output

	# 	fasta_merge -d ../{wildcards.sample}.maker.output/{wildcards.sample}_master_datastore_index.log
	# 	gff3_merge -d ../{wildcards.sample}.maker.output/{wildcards.sample}_master_datastore_index.log
	# 	"""

# ---------------------------
### Renaming
# ---------------------------

rule renameMAKERresults:
	input:
		gff = "MAKER/{sample}/MAKER_output/{sample}.all.gff",
		tempnames = "temp/genomes/{sample}_tempnames.txt"
	output:
		realnames = "MAKER/{sample}/MAKER_output/{sample}.all.illu.gff",
	run:
		# Make a dictionary
		names = {}

		with open(input.tempnames, 'r') as temps:
			for line in temps:
				cleanline = line.rstrip("\n").split("\t")
				names[cleanline[1]] = cleanline[0]

		# Keep the names somewhere
		ofile = open(output.realnames, 'w')

		with open(input.gff, 'r') as gff:
			for line in gff:
				if "#" in line:
					ofile.write(line)
				else:
					cleanline = line.split("\t")
					newname = cleanline[0].rstrip("\n").strip('>')

					if newname in names.keys():
						newline = line.replace(newname, names[newname])
						ofile.write(newline)
					else:
						ofile.write(line)

rule collectmakerresults:
	input:
		maker = typesample_maker,
		# maker = "MAKER/{sample}/MAKER_output/{sample}.all.gff",
	output:
		maker = "results/{sample}/{sample}.all.gff3",
		models = "renaming/{sample}/{sample}.models.gff3",
	shell:
		"""
		# Maker output
		cat {input.maker} | grep -v -P ".\\tcontig\\t" > {output.maker} # I dislike having the full contig when visualizing

		grep -P ".\\tmaker\\t" {output.maker} > {output.models}
		"""

rule newcodesgff:
	""" Sort and give new codes to the MAKER gene models to make them cleaner """
	input:
		namedgff = "renaming/{sample}/{sample}.models.gff3"
	output:
		sortednamedgff = "renaming/{sample}/{sample}.models.sort.gff3",
		panomgff = temp("renaming/{sample}/{sample}.models.sort.PaNom.gff3")
	params:
		time = "10:00",
		threads = 1,
		GFFnumerator = GFFnumerator,
	shell:
		"""
		## Sort
		igvtools sort {input.namedgff} {output.sortednamedgff}
	
		## Give new codes to the sorted models; notice the period after the sample name
		{params.GFFnumerator} {output.sortednamedgff} -s {wildcards.sample}. > {output.panomgff}  
		"""

rule getfastamaker:
	""" Get the CDS sequences of the MAKER models """
	input:
		gff = "renaming/{sample}/{sample}.models.sort.PaNom.gff3",
		genome = "data/genomes/{sample}.fa",
		# podan2gff = "data/Podan2/Podan2_AssemblyScaffoldsmtGenesEd.gff",
		# podan2genome= "data/Podan2/Podan2.AssemblyScaffoldsmt.fa", 
	output:
		cds = temp("renaming/{sample}/{sample}.models.sort.PaNom.fas"),
	params:
		gff3tofasta = gff3tofasta,
	shell:
		"""
		{params.gff3tofasta} {input.genome} {input.gff} --type noutrs --onlyids --output renaming/{wildcards.sample}/{wildcards.sample}.models.sort.PaNom > {input.gff}_weirdgenes.log
		"""

rule migratenames:
	""" Migrate the Podan2 names to the new annotation """
	input:
		makergenes = "renaming/{sample}/{sample}.models.sort.PaNom.fas",
		models = "renaming/{sample}/{sample}.models.sort.PaNom.gff3", # Use the unsorted GFF3 because IGV freaks out with the sorted one, and my script can handle it
	output:
		log = "renaming/{sample}/{sample}.models.sort.PaNom.gff3_ID.log",
		namedgff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.gff3"
	params:
		podan2 = refgenesfas, # The podan2 genes in the right order (fasta file)
		time = "40:00",
		GffgenesIDFix = GffgenesIDFix,
		periden = getperind,
		threads = 2,
		tempdir = "renaming/{sample}/" # This is to prevent colliding with other samples
	shell:
		"""
		{params.GffgenesIDFix} {params.podan2} {input.makergenes} {input.models} --identity {params.periden} --superverbose --clean --threads {params.threads} --output {output.namedgff} --temp {params.tempdir} > {output.log}
		"""
		# {params.GffgenesIDFix} {params.podan2} {input.makergenes} {input.models} --dbtype prot --identity {params.periden} --superverbose --clean --threads {params.threads} --output {output.namedgff} --temp {params.tempdir} > {output.log}

rule addgeneproducts:
	""" Add gene products to genes that were identified """
	input:
		gff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.gff3",
		geneproducts = geneproducts
	output:
		gff = temp("renaming/{sample}/{sample}.models.sort.PaNom.ID.prods.gff3"),
	params:
		gff3addproduct = gff3addproduct
	shell:
		"{params.gff3addproduct} {input.gff} {input.geneproducts} > {output.gff}"


rule sppcodenicegff:
	""" Change the species codes of the gene names """
	input:
		panomgff = "renaming/{sample}/{sample}.models.sort.PaNom.ID.prods.gff3"
	output:
		nicegff = "results/{sample}/{sample}.nice-" + annotationversion + ".gff3"
	params:
		sppcode = getsppcode,
	shell:
		"""
		sed "s;Pa_;{params.sppcode}_;" {input.panomgff} > {output.nicegff} 
		"""


# ---------------------------
### Collect results
# ---------------------------

rule renameRMresults:
	input:
		repeatmasker = "RepeatMasker/{sample}/{sample}.fa.out.gff",
		tempnames = "temp/genomes/{sample}_tempnames.txt"
	output:
		realnames = "RepeatMasker/{sample}/{sample}.repeatmasker.illu.gff",
	run:
		# Make a dictionary
		names = {}

		with open(input.tempnames, 'r') as temps:
			for line in temps:
				cleanline = line.rstrip("\n").split("\t")
				names[cleanline[1]] = cleanline[0]
		
		# Keep the names somewhere
		ofile = open(output.realnames, 'w')

		with open(input.repeatmasker, 'r') as gff:
			for line in gff:
				if "#" in line:
					ofile.write(line)
				else:
					cleanline = line.rstrip("\n").split("\t")
					cleanline[0] = names[cleanline[0]] # Replace the temporary name with the original
					newline = '\t'.join(cleanline) + '\n'
					ofile.write(newline)


rule collectRMresults:
	input:
		repeatmasker = typesample_rm,
	output:
		repeatmasker = "results/{sample}/{sample}.repeatmasker.gff",
	shell:
		"cp {input.repeatmasker} {output.repeatmasker}"

