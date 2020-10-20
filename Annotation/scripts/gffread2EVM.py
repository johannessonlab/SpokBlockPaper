#!/usr/bin/env python

# ================== EDIT Cofflinks/gffread GFF3 file for EVM compatibility =================
# Script to make the GFF3 output of Cofflinks + gffread compatible with
# EVidenceModeler by appending the missing ID to the exons. The header can be
# turn off. It's largely based on the script SNAPgff2_gff3.py.

# NOTE: Actually this script only makes the output of gffread a bit closer to
# GFF3 specifications. However, EVM requires a modified version of the
# original GFF3. To make it completely compatible with EVM, do:
# $ grep 'exon\t' output_of_gffread2EVM.py.gff3 | sed 's:ID=.*Parent:ID:' > gffreallycompatiblewithEVM.gff3

# 2018.12.07: Version 1.1 - update to python 3
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2014/07/31
# +++++++++++++++++++++++++++++++++++++++++++++++++
# python gffread2EVM.py InputFiles/Thamnolia_RNAcuff_MAT.gff3 -n > new.gff3

# ------------------------------------------------------
import sys
import re
import datetime
import argparse # For the fancy options
# ------------------------------------------------------
version = 1.1
versiondisplay = "{0:.2f}".format(version)

# Make a nice menu for the user
parser = argparse.ArgumentParser(description="* Make gffread (Cufflinks) output compatible with EVidenceModeler *", epilog="The output goes to stdout. Use ONLY for SNAP GFF2 output!") # Create the object using class argparse

# Add options
parser.add_argument('cuffgff', help="GFF3 obtained after using gffread on Cufflinks output")
parser.add_argument('--noheader', '-n', help="Turn off header", default=False, action='store_true')
parser.add_argument('--prefix', '-p', help="Prefix of the model names (default = 'CUFF')", default='CUFF')
parser.add_argument('--version', "-v", action='version', version='%(prog)s '+ versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	cuffgffopen = open(args.cuffgff, 'rU')
except IOError as msg:  # Check that the file exists
    parser.error(str(msg)) 
    parser.print_help()

# ------------------------------------------------------

# ---- Functions taken from Joshua Orvis' script convert_aat_btab_to_gff3.py 
next_ids = {'exon':1}
def get_next_id(type, prefix):
    if prefix == None:
        id = type + ".{0:06d}".format( next_ids[type] )
    else:
        id = prefix + "." + type + ".{0:06d}".format( next_ids[type] )
    
    next_ids[type] += 1

    return id
# ---- I got these functions from http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside by unutbu user
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
# ------------------------------------------------------

# -------- Print first lines to stdout --------
if not args.noheader: 				# Depending if the user wants the header or not
	# Get the date
	today = datetime.date.today()
	reportdate = "# Date " + str(today) + "\n"
	about = "# Modified GFF3 of Cufflinks plus gffread with the script " + sys.argv[0] + " " + versiondisplay + "\n"

	# Print to stdout
	sys.stdout.write("##gff-version 3\n")
	sys.stdout.write(reportdate)
	sys.stdout.write(about)

# -------- Read the GFF3 file --------
# Append an ID to the exons
try:
	for line in cuffgffopen:
		cols = line.rstrip("\n").split("\t")
		linestring = ''
		if '#' in cols[0]:
			pass  # Ignore comments in the file
		else:
			for tab in cols[0:8]: 						# Keep the all the fields except attribute as they are in tab format
				linestring += ''.join(str(tab)) + "\t"
			if cols[2] == 'exon' and 'ID' not in cols[8]:  # For exons only and double checking they don't have the ID already
				exon_id = get_next_id('exon', args.prefix)  # Get new exon ID
				attribute = "ID={0};{1}".format(exon_id, cols[8])  # Build new ID
				finalline = linestring + attribute + "\n" 
				sys.stdout.write(finalline)
			else:
				sys.stdout.write(line)
except IOError:  	# Catch some errors if piping in the command line
	try:
	    sys.stdout.close()
	except:
	    pass
	try:
	    sys.stderr.close()
	except:
	    pass