#!/usr/bin/env python
# encoding: utf-8

# ================== gff3addproduct =================
# Script to add the product name of a gene (from Podan2) to a gff file
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/11/13
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import datetime
# ------------------------------------------------------

version = 1
versiondisplay = "{0:.2f}".format(version)

# Input from console
try:
	GFFopen = open(sys.argv[1], 'r')
	products = open(sys.argv[2], 'r')
except:
	print("Usage: python " + sys.argv[0] + " annotation.gff products.txt")
	print("Version " + versiondisplay)
	sys.exit(1)

# --------------

# Make a dictionary with the products file
proddic = {}
for line in products:
	cols = line.rstrip("\n").split("\t")
	proddic[cols[0]] = cols[1]

# Process the GFF
for line in GFFopen:
	if '##gff-version 3' in line: 
		sys.stdout.write(line) # The very fist line
		# Add a line to mark the file with this script
		now = datetime.datetime.now()
		sys.stdout.write("# Gene products added with " + os.path.basename(sys.argv[0]) + " v. " + str(versiondisplay) + ' on ' + str(now) + '\n')

	elif '#' in line: # Print comment as it is
		sys.stdout.write(line)
	elif line not in ['\n', '\r\n']: # Ignore empty lines
		if ('\tgene\t' in line): # Look for the genes definition
			cols = line.rstrip("\n").split("\t")		# break the line into columns defined by the tab
			attributes = cols[8].split(";") # Get individual attributes

			genename = [element for element in attributes if "Name=" in element][0] # The name of this gene
			
			for gene in proddic.keys(): # Is it in the known products dictionary?
				fullgenestr = "Name=" + gene
				if fullgenestr == genename: # to avoid naming genes that overlap, like Pa_5_10 and Pa_5_101
					line = line.rstrip("\n") + ";product_name=" + proddic[gene] + "\n" # Replace the current line
		sys.stdout.write(line)
