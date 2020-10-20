#!/bin/bash

#SBATCH -A XXXXXXX
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH --mail-user xxxxxxx@xxxx
#SBATCH --mail-type=ALL

# ================ runTransDecoder =================
# Script to run TransDecoder on the output of STAR + Cufflinks in order to get
# candidate coding regions within transcript sequences.

# http://transdecoder.github.io/
# For now I won't include homology searches
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2017/01/01
# +++++++++++++++++++++++++++++++++++++++++++++++++


# --- Variables ---
VERSION="2.0"
# -----------------

TRANSCRIPTGTF=$1
GENOME=$2

if [[ -z "$TRANSCRIPTGTF" || -z "$GENOME" ]] ; then 
        echo 'Must specify a GTF file (e.g. from Cufflinks), the genome used for mapping and a name for output (optional)'
        echo "Usage: $0 transcripts.gtf genome.fas [ sampleID ]"
        echo "Version $VERSION"
        exit 1
fi

# **** The number of cores used ****
# Set number of cores based on default sbatch command or by user
if [[ $SLURM_CPUS_ON_NODE ]] ; then # If the variable exists
	CORES=$SLURM_CPUS_ON_NODE
elif [[ condition ]]; then
	CORES=1
fi
# **********************************

# A name for the sample
basegenome=$(basename ${TRANSCRIPTGTF%.gtf*})
sampleID=${3:-$basegenome}

echo "Transcripts in gtf: $TRANSCRIPTGTF"
echo "Output basename: $sampleID"
echo "$0 Version $VERSION"

echo " ... Construct the transcript fasta file using the genome and the transcripts.gtf file ... "
gtf_genome_to_cdna_fasta.pl $TRANSCRIPTGTF $GENOME > $sampleID.fa
# $TransDecoderPATH/util/cufflinks_gtf_genome_to_cdna_fasta.pl $TRANSCRIPTGTF $GENOME > $sampleID.fa

echo " ... Convert the transcript structure GTF file to an alignment-GFF3 formatted file ... "
gtf_to_alignment_gff3.pl $TRANSCRIPTGTF > $sampleID.gff3
# $TransDecoderPATH/util/cufflinks_gtf_to_alignment_gff3.pl $TRANSCRIPTGTF > $sampleID.gff3

echo " ... Generate the best candidate ORF predictions ... "
TransDecoder.LongOrfs -S -t $sampleID.fa
# $TransDecoderPATH/TransDecoder.LongOrfs -S -t $sampleID.fa

# By default, TransDecoder.LongOrfs will identify ORFs that are at least 100
# amino acids long. You can lower this via the -m parameter, but know that the
# rate of false positive ORF predictions increases drastically with shorter
# minimum length criteria.

# If the transcripts are oriented according to the sense strand, then include
# the -S flag to examine only the top strand.

echo " ... Predict the likely coding regions ... " # This takes the longest time!
TransDecoder.Predict -t $sampleID.fa --cpu $CORES

echo " ... Generate a genome-based coding region annotation file ... "
cdna_alignment_orf_to_genome_orf.pl $sampleID.fa.transdecoder.gff3 $sampleID.gff3 $sampleID.fa > $sampleID.genome.transdecoder.gff3




