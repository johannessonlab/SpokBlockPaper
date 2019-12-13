# -*- snakemake -*-

### Recovering alignments for manual curation of TE libraries of Podospora spp
#############################################################################
# A Snakemake pipeline to produce alignments of a TE library ready for manual curation
# This is an adaptation of the Alex Suh Lab pipeline.
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/03/31
# ++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
### Variables
rmlibpath = config["rmlibpath"]

repeatModelerPipeline = config["repeatModelerPipeline"]

AllSamples = config["samplesID"]
# AllSamples = "Podan2"
# -------------------------------------------------

# ----------
# Rules not submitted to a job
localrules: blastdbs, cleanlib, makealignments, lessgapsalignment # If I ask for the modules, then it wants parameters too
# ----------

rule all:
	input:
		# RM results
		expand("repeatModelerPipeline/{sample}/{sample}.dummy", sample = AllSamples),
		expand("repeatModelerPipeline/{sample}/final/{sample}.dummy", sample = AllSamples),
		expand("repeatModelerPipeline/{sample}/{sample}_gap.dummy", sample = AllSamples),
# ---------------------------------

# ---------------
rule blastdbs:
	""" Make a database of the Nanopore and SPAdes* assemblies """
	input:
		assembly = "data/genomes/{sample}.fa",
	output:
		temp("databases/{sample}_db/{sample}_db.nhr"),
		temp("databases/{sample}_db/{sample}_db.nin"),
		temp("databases/{sample}_db/{sample}_db.nog"),
		temp("databases/{sample}_db/{sample}_db.nsd"),
		temp("databases/{sample}_db/{sample}_db.nsi"),
		temp("databases/{sample}_db/{sample}_db.nsq"),

	shell:
		"""
		echo "Making database..."
		makeblastdb -in {input.assembly} -out databases/{wildcards.sample}_db/{wildcards.sample}_db -dbtype nucl -parse_seqids
		"""

rule cleanlib:
	""" Make a copy of consensi.fa.classified.clean without "/" to be compatible with repeatModelerPipeline (sequence names are used as alignment file names) """
	input:
		rmlib = rmlibpath + "/{sample}/{sample}_RM.lib",
	output:
		cleanrmlib = "data/libraries/{sample}_RM.lib.clean",
		# cleanrmlibnoslash = "data/libraries/{sample}/{sample}_RM.lib.clean_noSlash"
	shell:
		"""
		# First remove long header from lib
		awk '{{print$1;}}' {input.rmlib} | sed 's|/|_|g' > {output.cleanrmlib} 
		"""

rule repeatModelerPipeline:
	""" Run Alex Suh Lab's pipeline """
	input:
		genome = "data/genomes/{sample}.fa",
		database = "databases/{sample}_db/{sample}_db.nhr",
		rmlib = "data/libraries/{sample}_RM.lib.clean",
		# rmlib = rmlibpath + "/{sample}/{sample}_RM.lib",

		# Other files, so they get used and then removed by the temp file
		databasenhr = "databases/{sample}_db/{sample}_db.nhr",
		databasenin = "databases/{sample}_db/{sample}_db.nin",
		databasenog = "databases/{sample}_db/{sample}_db.nog",
		databasensd = "databases/{sample}_db/{sample}_db.nsd",
		databasensi = "databases/{sample}_db/{sample}_db.nsi",	
		databasensq = "databases/{sample}_db/{sample}_db.nsq",
	output:
		"repeatModelerPipeline/{sample}/{sample}.dummy"
	params:
		repeatModelerPipeline = repeatModelerPipeline,
		time = "10:00:00", # Fast for Podospora (<20 TEs and 36Mb genome)
		threads = 16,
	shell:
		"""
		cd repeatModelerPipeline/{wildcards.sample}

		# Run the pipeline
		perl {params.repeatModelerPipeline} ../../{input.genome} ../../databases/{wildcards.sample}_db/{wildcards.sample}_db ../../{input.rmlib} && touch {wildcards.sample}.dummy
		
		# Cleanup of temporary files:
		rm *emp.out
		"""

rule makealignments:
	""" Remove redundant sequences from alignments and save in folder "final" """ 
	input:
		"repeatModelerPipeline/{sample}/{sample}.dummy"
	output:
		"repeatModelerPipeline/{sample}/final/{sample}.dummy"
	shell:
		""" 
		cd repeatModelerPipeline/{wildcards.sample}/aligned

		for sample in `ls *.fa`; do 
			name=${{sample%.fa}}
			cat $sample | perl -ne 'chomp;s/>\s+/>/;if(/>(\S+)/){{$id{{$1}}++;$id2=$1;}}if($id{{$id2}}==1){{print "$_\\n"}}' > ../final/$name.fa 
		done && touch ../final/{wildcards.sample}.dummy
		"""

rule lessgapsalignment:
	""" Remove alignment positions with gaps in more than 90% of the sequences """ 
	input:
		"repeatModelerPipeline/{sample}/final/{sample}.dummy"
	output:
		"repeatModelerPipeline/{sample}/{sample}_gap.dummy",
	conda: 
		"envs/tcoffee.yaml"
	shell:
		""" 
		cd repeatModelerPipeline/{wildcards.sample}/final
		export CACHE_4_TCOFFEE=../../../temp/{wildcards.sample} # this is where T-Coffee stores any data expensive to obtain: PDB files, structural alignments ...

		for file in `ls *.fa`; do 
			name=${{file%.fa}}
			t_coffee -other_pg seq_reformat -in $file -action +rm_gap 90 > $name.nogaps90.fa
		done && touch ../{wildcards.sample}_gap.dummy
		"""






