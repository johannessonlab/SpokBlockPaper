#!/usr/bin/env python
# encoding: utf-8

# ================== GFFnumerator =================

# Script to re-name the IDs of the genes and features of a GFF3 file. It
# prints to the standard output, and it's designed for python3. It depends on
# the library gffutils.

# Sources
# https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# http://daler.github.io/gffutils/database-ids.html#merge-strategy
# https://github.com/daler/gffutils/blob/master/gffutils/test/test.py
# https://github.com/daler/gffutils/issues/61
# https://pythonhosted.org/gffutils/database-schema.html
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2019/03/26
# +++++++++++++++++++++++++++++++++++++++++++++++++

# ------------------------------------------------------
import sys # To exit the script 
import os # For the input name
import argparse  # For the fancy options
# import time
import datetime
import gffutils
import gffutils.inspect as inspect
# ------------------------------------------------------

version = 1
versiondisplay = "{0:.2f}".format(version)

# ============================
# Make a nice menu for the user
# ============================
parser = argparse.ArgumentParser(description="* Script to re-name the IDs of the genes and features of a GFF3 file *")  # Create the object using class argparse

# Add options
parser.add_argument('GFF', help="GFF3 file sorted beforehand")
parser.add_argument('--sample', '-s', help="String representing sample that gets appended into the gene IDs")
parser.add_argument("--printdb", "-p", help="Print the database into a file instead of saving it in memory", default=False, action='store_true')
parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + versiondisplay)

try:
	# ArgumentParser parses arguments through the parse_args() method You can
	# parse the command line by passing a sequence of argument strings to
	# parse_args(). By default, the arguments are taken from sys.argv[1:]
	args = parser.parse_args()
	GFFopen = open(args.GFF, 'r')
except IOError as msg:  # Check that the file exists
	parser.error(str(msg)) 
	parser.print_help()
# ============================

# ---------------------------------
# Make database
# ---------------------------------
# t0 = time.time()
# This will parse the file, infer the relationships among the features in the file, and store the features and relationships
# See https://pythonhosted.org/gffutils/autodocs/gffutils.create_db.html
# I expect a MAKER gff, where CDS have no unique ID

if args.printdb:
	input_base = os.path.splitext(args.GFF)[0] # Taking out the prefix of the file
	input_name = os.path.basename(input_base) # Remove the path
	dbfnchoice = input_name + '.db'
else:
	dbfnchoice = ':memory:'

id_spec={"gene": ["ID", "Name"], "mRNA": ["ID", "transcript_id"]} # http://daler.github.io/gffutils/database-ids.html

db = gffutils.create_db(data = args.GFF, 
	dbfn = dbfnchoice,
	force = True, # force=True overwrite any existing databases.
	id_spec = id_spec, 
	merge_strategy = "create_unique") # Add an underscore an integer at the end for each consecutive occurrence of the same ID 
	# verbose = True,) 

# t1 = time.time()
# db_results = inspect.inspect(db) # Report
# print("\n\nIt took {0:.1f}s to create database".format(t1 - t0))
# ---------------------------------

## Get iterator of all genes
genes = [gene for gene in db.features_of_type("gene")]
nogenes = db.count_features_of_type("gene")

# Get ids of genes in database
allIDsGenes = [gene.id for gene in db.features_of_type("gene")]

## Make new IDs for the genes
if args.sample is not None:
	newIDgene = [args.sample + "{0:05d}".format(n) for n in range(1, nogenes + 1)]
else:
	newIDgene = ["gene.{0:05d}".format(n) for n in range(1, nogenes + 1)]

# ---------------------------
### Change IDs
# ---------------------------

def getnewID(id):
	indexgene = allIDsGenes.index(id)
	newid = newIDgene[indexgene]
	# print(newid)
	return(newid)

def gen():
	sys.stdout.write("##gff-version 3\n")
	# Add a line to mark the file with this script
	now = datetime.datetime.now()
	newhead = '# Original file ' + os.path.basename(args.GFF) + ' modified with GFFnummerator.py v. ' + str(versiondisplay) + ' on ' + str(now) +  '\n'
	sys.stdout.write(newhead)

	# Actual GFF
	for gene in db.features_of_type('gene'):
	# for gene in db.features_of_type('gene', order_by='start'):
		newidgene = getnewID(gene.id)
		gene['ID'] = newidgene

		print(gene)

		mRNA = 1
		for child in list(db.children(gene)):	
			if child.featuretype == "mRNA":
				child['Parent'] = newidgene

				# make a new ID for the mRNA
				newidrna = newidgene + "-mRNA" + str(mRNA)
				mRNA += 1
				child['ID'] = newidrna

				print(child)

				# Update all the features that depend on the mRNA, assuming each feature has a SINGLE parent
				typeids = {'gene':1, 'mRNA':1, 'exon':1, 'CDS':1, 'five_prime_UTR':1, 'three_prime_UTR':1}
				for grandchild in list(db.children(gene, level = 2)):
					grandchild['Parent'] = newidrna # TODO Add something here, a conditional, in case there are multiple parents

					typefea = grandchild.featuretype # What type are we dealing with?

					newidchild = newidrna + '-' + typefea + "{0:01d}".format(typeids[typefea])
					grandchild['ID'] = newidchild

					typeids[typefea] += 1 # Increase the count of the corresponding type for this gene
					print(grandchild)
gen()
