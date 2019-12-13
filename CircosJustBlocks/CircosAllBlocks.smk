# -*- snakemake -*-

### Circos plot of the Spok blocks
#############################################################################
#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/12/05
# +++++++++++++++++++++++++++++++++++++++++++++++++
# Version 1

# -------------------------------------------------
allblocks = config["allblocks"]
samples = config["samples"]
gffs = config["gffs"]

# -------------------------------------------------

rule all:
	input:
		"circos/circos.png",
		

rule makekaryotype:
	""" Produce a karyotype file for circos """
	input:
		allblocks
	output:
		"circos/karyotype.txt"
	params:
		threads = 1,
	shell:
		"""
		samtools faidx {input}

		cat {input}.fai | awk {{'print "chr","-",$1,$1,1,$2,"vdgrey"'}} > {output}
		rm {input}.fai
		"""	

rule allvsall:
	input:
		query = allblocks,
		reference = allblocks,
	output:
		delta = "mummer/allblocks.delta",
		deltafilter = "mummer/allblocks.filter",
		coords = "mummer/allblocks.coords",
		coordsfilter = "mummer/allblocks.filter.coords",	
	params:
		time = "1:00:00",
		threads = 8,
	shell:
		"""
		echo
		echo "MUMmer alignment ..."
		nucmer -b 200 -c 40 --maxmatch -p mummer/allblocks {input.reference} {input.query} -t {params.threads}

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
        # -c|mincluster   Sets the minimum length of a cluster of matches (default 65)
        # -b|breaklen     Set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default 200)


rule makelinks:
	""" Prepare a links with the MUMmer alignments """
	input:
		coord = "mummer/allblocks.coords",
		# karyotype = "circos/karyotype.txt"
	output:
		allsites = "links/mummer.txt",
	shell:
		"""
		cat {input.coord} | awk '{{i=$10; if ( i == "PcWa139m" ) {{print $10,$1,$2,$11,$3,$4,"color=vdgreen_a2"}} else {{print $10,$1,$2,$11,$3,$4,"color=vlpurple_a5"}} }}' > {output}
		# cat {input.coord} | awk '{{i=$10; if ( i == "PcWa139m_SpokBlock" ) {{print $10,$1,$2,$11,$3,$4,"color=vdgreen_a2"}} else {{print $10,$1,$2,$11,$3,$4,"color=vlpurple_a5"}} }}' > {output}
		"""
		# http://circos.ca/documentation/tutorials/recipes/complex_histograms/configuration
		# http://circos.ca/documentation/tutorials/recipes/heatmap_links/
# ------------ Annotation ------------

def getgff(wildcards):
	for file in gffs:
		if wildcards.sample in file:
			return(file)
			break

## TEs
rule TEtrack:
	""" Prepare data for TE track """
	input:
		getgff
	output:
		temp("tracks/{sample}_TEs.txt")
	params:
		threads = 1,
	shell:
		"""
		# Get TEs into a track format
		grep 'RepeatMasker' {input} | awk -v seq={wildcards.sample} '{{ print seq,$4,$5,"id=TE" }}' > {output} 
		"""

rule mergetracks_TE:
	input:
		expand("tracks/{sample}_TEs.txt", sample = samples)
	output:
		"tracks/TEs.txt"
	params:
		threads = 1,
	shell:
		"cat {input} > {output}"

## Genes
rule genetrack:
	""" Prepare data for genes track """
	input:
		getgff
	output:
		temp("tracks/{sample}_genes.txt")
	params:
		threads = 1,
	shell:
		"""
		# Get genes into a track format
		grep -P '\\tgene\\t' {input} | awk -v seq={wildcards.sample} '{{ann=$9; label=gensub(/ID=(.+)=(Pa_)?(.+)/, "\\\\2", "g", ann); if( label == "Pa_") {{print seq,$4,$5,"fill_color=orange"}} else {{print seq,$4,$5,"fill_color=dpurple"}} }}' > {output} 
		"""
		# We used CDS instead of genes because Kirc is very broken by TEs in PcWa139m

rule mergetracks_genes:
	input:
		expand("tracks/{sample}_genes.txt", sample = samples)
	output:
		"tracks/genes.txt"
	params:
		threads = 1,
	shell:
		"cat {input} > {output}"

# ------------ Run Circos ------------

rule circos:
	""" Run Circos to plot the alignment """
	input:
		"circos/etc/circos.conf",
		"circos/karyotype.txt",
		"links/mummer.txt",
		"tracks/TEs.txt",
		"tracks/genes.txt",
	output:
		"circos/circos.png"
	shell:
		"cd circos; circos"

