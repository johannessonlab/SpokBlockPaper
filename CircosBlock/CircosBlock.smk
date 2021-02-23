# -*- snakemake -*-

### Circos plot of the Spok block vs the host genome
#############################################################################
#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/11/28 - 2019/12/13 - 2021/02/23
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
block = config["block"]
host = config["host"]
otherblocks = config["otherblocks"]
blockgff = config["blockgff"]

hostTEs = config["hostTEs"]
hostgenes = config["hostgenes"]
# -------------------------------------------------

rule all:
	input:
		"circos/circos.png",
		"results/TEcontent.txt" # A table of repeat content in the reference genome and block, for funsies

# ------------ Karyotype ------------

rule makekaryotype:
	""" Produce a karyotype file for circos """
	input:
		host = host,
		block = block
	output:
		"circos/karyotype.txt"
	params:
		threads = 1,
	shell:
		"""
		samtools faidx {input.host}
		samtools faidx {input.block}

		cat {input.host}.fai | awk '{{seq=$1; label=gensub(/(.+)chromosome_([\.0-9]*)/, "Chr\\\\2", "g", seq); print "chr","-",$1,label,1,$2,"grey"'}} > {output}
		cat {input.block}.fai | awk {{'print "chr","-",$1,$1,1,$2,"vdgrey"'}} >> {output}

		rm {input.host}.fai
		rm {input.block}.fai
		"""	

# ------------ BLAST ------------
rule blastdbs:
	""" Make a database of the Nanopore and SPAdes* assemblies """
	input:
		assembly = ancient(host), # Ignoring timestamps
	output:
		"databases/host_db/host_db.nhr",
		"databases/host_db/host_db.nin",
		"databases/host_db/host_db.nog",
		"databases/host_db/host_db.nsd",
		"databases/host_db/host_db.nsi",
		"databases/host_db/host_db.nsq",
	params:
		# time = "30:00",
		threads = 1,
	shell:
		"""
		echo "Making database..."
		makeblastdb -in {input.assembly} -out databases/host_db/host_db -dbtype nucl -parse_seqids
		"""

rule blastn:
	""" Use BLASTn to get nucleotide alignments """
	input:
		query = ancient(block),
		database = "databases/host_db/host_db.nhr",
	output:
		blast = "BLAST/block-to-host_tabularn.txt",
	params:
		# time = "30:00",
		threads = 6,
	version: "2"
	shell:
		"""
		outputblast=$(echo {output.blast}| sed 's/.txt//')

		database=$(echo {input.database}| sed 's/.nhr//')

		echo "Blasting ..."
		# The entire tabular file
		tblastx -query {input.query} -out {output.blast} -db $database -outfmt 6 -num_threads {params.threads}
		"""

rule makelinks_blast:
	""" Prepare a links with the tabular BLAST hits """
	input:
		"BLAST/block-to-host_tabularn.txt",
	output:
		"links/blast.txt",
	shell:
		"""
		cat {input} | awk '{{i=$3; if ( i >= 70 ) {{print $1,$7,$8,$2,$9,$10,"color=blue_a5"}} else {{print $1,$7,$8,$2,$9,$10,"color=vlpurple_a5"}} }}' > {output}
		# cat {input} | awk '{{print $1,$7,$8,$2,$9,$10,"value="$3}}' > {output}
		"""

# ------------ MUMmer ------------

rule mummer_block_host:
	""" Make nucmer alignments """
	input:
		query = block,
		reference = host,
	output:
		delta = "mummer/block-to-host.delta",
		deltafilter = "mummer/block-to-host.filter",
		coords = "mummer/block-to-host.coords",
		coordsfilter = "mummer/block-to-host.filter.coords",	
	params:
		time = "1:00:00",
		threads = 8,
	shell:
		"""
		basename=$(echo {output.delta} | cut -f1 -d".")
		echo $basename

		echo
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 40 --maxmatch --nosimplify -p $basename {input.reference} {input.query} -t {params.threads}

		# Filter the delta
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# --mum  Use anchor matches that are unique in both the reference and query
		# --mumreference  Use anchor matches that are unique in in the reference
        #           but not necessarily unique in the query (default behavior)
        # --maxmatch      Use all anchor matches regardless of their uniqueness (false)
        # --nosimplify    Don't simplify alignments by removing shadowed clusters. 
        #           Use this option when aligning a sequence to itself to look for repeats (false)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)


rule makelinks_mummer:
	""" Prepare a links with the MUMmer alignments """
	input:
		"mummer/block-to-host.coords",
	output:
		allsites = "links/mummer.txt",
	params:
		threads = 1,
	shell:
		"""
		cat {input} | awk '{{ print $10,$1,$2,$11,$3,$4,"color=green_a5"}}' > {output}
		"""
		# http://circos.ca/documentation/tutorials/recipes/complex_histograms/configuration
		# http://circos.ca/documentation/tutorials/recipes/heatmap_links/

# ------------ Annotation ------------

rule blockTEtrack:
	""" Prepare data for TE track """
	input:
		block = block,
		blockgff = blockgff,
	output:
		"tracks/BlockTEs.txt"
	params:
		threads = 1,
	shell:
		"""
		# Get name of the block sequence
		seqname=$(grep '>' {input.block} | sed 's/>//')
		
		# Get TEs into a track format
		grep 'RepeatMasker' {input.blockgff} | awk -v seq="$seqname" '{{ print seq,$4,$5,"id=TE" }}' > {output} 
		"""

rule blockGenesTrack:
	""" Prepare data for the genes track """
	input:
		block = block,
		blockgff = blockgff,
	output:
		"tracks/BlockGenes.txt"
	params:
		threads = 1,
	shell:
		"""
		# Get name of the block sequence
		seqname=$(grep '>' {input.block} | sed 's/>//')
		
		# Get genes into a track format
		grep -P '\\tgene\\t' {input.blockgff} | awk -v seq="$seqname" '{{ann=$9; label=gensub(/ID=(.+)=(Pa_[0-9]_[0-9]+)?(_.+)/, "\\\\2", "g", ann); if( label ~ "Pa_") {{print seq,$4,$5,"fill_color=orange"}} else {{print seq,$4,$5,"fill_color=dpurple"}} }}' > {output} 
		# x ~ y   True if the string x matches the regexp denoted by y 
		"""
# ------------ Distribution of TEs in the host ------------

rule hostbed:
	""" Create and index-like bed file for the host """
	input:
		host
	output:
		temp("tracks/host.bed")
	params:
		threads = 1,
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print $1"\\t"1"\\t"$2'}} > {output} # It has to be tabs or BEDtools won't like it.
		rm {input}.fai 
		"""

rule makewindows:
	""" Use the BEDtools makewindows to get windows of the host genome """
	input:
		"tracks/host.bed"
	output:
		temp("tracks/hostwins.bed")
	params:
		threads = 1,
	shell:
		"""
		bedtools makewindows -b {input} -w 50000 -s 10000 > {output}
		"""	

rule TEhostbed: 
	""" Prepare data for TE bed from the host """
	input:
		hostTEs
	output:
		temp("tracks/hostTEs.bed")
	params:
		threads = 1,
	shell:
		"""
		grep 'RepeatMasker' {input} | awk {{'print $1"\\t"$4"\\t"$5'}} > {output}
		"""

rule BEDtools_TEs:
	""" Use BEDtools coverage to produce a distribution of TEs along the host genome """
	# https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
	input:
		chrs = "tracks/hostwins.bed",
		tes = "tracks/hostTEs.bed",
	output:
		"tracks/hostTEdistribution.txt"
	params:
		threads = 1,
	shell:
		"""
		bedtools coverage -a {input.chrs} -b {input.tes} | awk '{{ print $1,$2,$3,$7 }}' > {output}
		"""

# ------------ Distribution of genes in the host ------------

rule TEgenesbed: 
	""" Prepare data for TE bed from the host """
	input:
		hostgenes
	output:
		"tracks/hostgenes.bed"
	params:
		threads = 1,
	shell:
		"""
		grep -P '\\tgene\\t' {input} | grep -vP '_mt\\t' | awk {{'print $1"\\t"$4"\\t"$5'}} > {output}
		"""

rule BEDtools_genes:
	""" Use BEDtools coverage to produce a distribution of TEs along the host genome """
	input:
		chrs = "tracks/hostwins.bed",
		tes = "tracks/hostgenes.bed",
	output:
		"tracks/hostgenesdistribution.txt"
	params:
		threads = 1,
	shell:
		"""
		bedtools coverage -a {input.chrs} -b {input.tes} | awk '{{ print $1,$2,$3,$7 }}' > {output}
		"""

# ------------ Conservation of the Spok block ------------

rule mummer_blocks:
	""" Make nucmer alignments of blocks """
	input:
		query = otherblocks,
		reference = block,
	output:
		delta = "mummer/blocks-to-TheBlock.delta",
		deltafilter = "mummer/block-to-TheBlock.filter",
		coords = "mummer/block-to-TheBlock.coords",
		coordsfilter = "mummer/block-to-TheBlock.filter.coords",	
	params:
		time = "1:00:00",
		threads = 8,
	shell:
		"""
		basename=$(echo {output.delta} | cut -f1 -d".")
		echo $basename

		echo
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 40 --mum -p $basename {input.reference} {input.query} -t {params.threads}

		# Filter the delta
		delta-filter -q {output.delta} > {output.deltafilter}

		# To view a summary of all the alignments produced by NUCmer
		echo "Running show-coords"
		# For Ribbon http://genomeribbon.com/
		echo "...for Ribbon"
		show-coords -r -lTH {output.delta} > {output.coords}
		show-coords -r -lTH {output.deltafilter} > {output.coordsfilter}

		"""
		# Default: Use anchor matches that are unique in in the reference
	    # --mum           Use anchor matches that are unique in both the reference and query (false)
	    # --maxmatch      Use all anchor matches regardless of their uniqueness (false)
        # --nosimplify    Don't simplify alignments by removing shadowed clusters. 
        #           Use this option when aligning a sequence to itself to look for repeats (false)
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)

rule mummer2bed:
	""" Process the coords file into a bed """
	input:
		"mummer/block-to-TheBlock.coords"
	output:
		temp("tracks/Blocks.bed")
	params:
		threads = 1,
	shell:
		"""
		cat {input} | awk '{{ print $10"\\t"$1"\\t"$2"\\t"$11 }}' > {output} # It has to be tabs or BEDtools won't like it.
		"""

# ----
def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.
	return(result)
# ----

rule removeoverlap:
	""" Fuse overlapping ranges in a BED file """
	# Basically the same as the script totalcovergff.py but for bed files
	# This is necessary because there are multiple alignments of repetitive regions in each block
	input:
		"tracks/Blocks.bed"
	output:
		"tracks/Blocks_collapsed.bed"
	run:
		# Make a dictionary of ranges
		beddic = {} # key: contig, value: (start, end)

		with open(input[0], 'r') as file: 
			for line in file:
				if '#' in line:
					pass
				elif line not in ['\n', '\r\n']: # Ignore empty lines
					cols = line.rstrip("\n").split("\t")

					contig = cols[3] # The contig that matched the reference block
					start = int(cols[1])
					end = int(cols[2])

					if contig in list(beddic.keys()):
						beddic[contig].append([start, end])
					else: # contig is new
						beddic[contig] = [[start, end]]

					# Get name of reference block (this will happen many times, so it's inefficient)
					ref = cols[0]

			# Reduce the overlaps
			for ctg in beddic.keys():
				beddic[ctg] = remove_overlap(beddic[ctg])

			# Print them in a new file
			with open(output[0], 'w') as result:
				result.write('#Ref\tStart\tEnd\tQuery\n') # header
				for ctg in beddic.keys():
					for interval in beddic[ctg]:
						result.write('{}\t{}\t{}\t{}\n'.format(ref, interval[0], interval[1], ctg))


rule blockbed:
	""" Create and index-like bed file for the host """
	input:
		block
	output:
		temp("tracks/block.genome")
		# temp("tracks/block.bed")
	params:
		threads = 1,
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print $1"\\t"$2'}} > {output} # It has to be tabs or BEDtools won't like it.
		rm {input}.fai
		"""

rule genomecov_blocks:
	""" Use BEDtools genomecov to produce a distribution of conservation along the reference block """
	input:
		block = "tracks/block.genome",
		# block = "tracks/blockwins.bed",
		otherblocks = "tracks/Blocks_collapsed.bed",
	output:
		"tracks/blocksconservation.txt"
	params:
		threads = 1,
	shell:
		"""
		bedtools genomecov -bga -g {input.block} -i {input.otherblocks} | awk '{{ print $1,$2,$3,$4 }}' > {output}
		"""


# ------------ Make a table to report the content ------------

rule tableTEcontent:
	input:
		chrs = "tracks/host.bed",
		tes = "tracks/hostTEs.bed",
		block = "tracks/block.bed",
		tesblock = "tracks/BlockTEs.txt",
	output:
		tesblock = temp("temp/BlockTEs.bed"),
		tecontent = "results/TEcontent.txt"
	params:
		threads = 1,
	shell:
		"""
		# Replace white spaces with tabs
		cat {input.tesblock} | awk '{{ print $1"\\t"$2"\\t"$3 }}' > {output.tesblock}
	
		# Genome
		bedtools coverage -a {input.chrs} -b {input.tes} | awk '{{ print $1,$2,$3,$7 }}' > {output.tecontent}
		
		# Block
		bedtools coverage -a {input.block} -b {output.tesblock} | awk '{{ print $1,$2,$3,$7 }}' >> {output.tecontent}
		"""	

# ------------ Run Circos ------------

rule circos:
	""" Run Circos to plot the alignment """
	input:
		"circos/karyotype.txt",
		"circos/etc/circos.conf",
		"links/mummer.txt",
		"tracks/BlockTEs.txt",
		"tracks/BlockGenes.txt",
		"tracks/hostTEdistribution.txt",
		"tracks/hostgenesdistribution.txt",
		"tracks/blocksconservation.txt",
		"links/blast.txt",
	output:
		"circos/circos.png"
	shell:
		"cd circos; circos"