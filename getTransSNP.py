import numpy as np 
import json 
import csv
import os 
import sys 
import gzip
import argparse
from collections import OrderedDict, defaultdict

def argsParser():

	# Define argument parser
	parser = argparse.ArgumentParser(description='A script that cleans the  oxforf gen file')

	# Add arguments
	parser.add_argument('-ChrIndex', help= 'int: number of the chromosome to be parsed. Only for slurm', required=False)

	# Initialize variables
	args = parser.parse_args()

	# If no CHR specified, analyze as specified in options
	if not args.ChrIndex:
		# Load options
		with open("options.json", "r") as jsonFile:
			options = json.load(jsonFile)
		if isinstance(options['chrArray']['getSNP'], int):
			chrArray = ["CHR_" + str(options)]
		elif 'all' in options['chrArray']['getSNP']: 
			chrArray = ["CHR_" + str(s) for s in range(1,23)]

	# Else, analyze the inputed chromosome (slurm)
	elif len(args.ChrIndex) > 1:
		chrArray = ["CHR_" + str(args.ChrIndex)]
		print(chrArray)
	elif len(args.ChrIndex) == 1:
		chrArray = ["CHR_" + str(args.ChrIndex)]
		print(chrArray)

	return chrArray

def extractPos(options, chr):

	# Open file 
	print("Opening " + options['folder']['Data'] + chr +'.csv.gz')
	OpenFile = gzip.open(options['folder']['Data'] + chr +'.csv.gz')
	next(OpenFile)

	# Create dict where key is a position
	SNPS = defaultdict(list)
	line = 1
	print("Processing file")
	for snp in OpenFile:

		try:
			if line % 10000 == 0:
				print('Current line: ' +  str(line))

			# Filter by pvalue (continue if bigger than pvalue threshold) and pair of genes
			if float(str(snp).split(',')[3]) > options['getTrans']['pval_thres'] or '&' in str(snp).split(',')[1]:
				line += 1
				continue

			# Take position 
			if 'rs' in str(snp).split(',')[0]: 
				pos = str(snp).split(',')[0].split(':')[-1]
			else:
				pos = str(snp).split(',')[0][2:]

			# Using dicts is too slow, use
			SNPS[pos].append([str(snp).split(',')[1], str(snp).split(',')[3]])

		except:
			continue

		line += 1

	return SNPS 


def filterSNPS(SNPS):

	print('Filtering file')

	# For each position
	line = 1
	print("Filtering")
	for pos in list(SNPS.keys()):

		# Verbose
		if line % 10000 == 0:
			print(str(line))
		line += 1

		# Get pvals 
		pvals = np.asarray([pval[-1] for pval in SNPS[pos]])

		# Get idx
		minpval = np.argmin(pvals)

		# Keep 
		SNPS[pos] = SNPS[pos][minpval]

	return SNPS

def main(chrArray):

	for chr in chrArray:

		# Read options
		with open('options.json','r') as inFile:
			options = json.load(inFile)

		# Extract positions 
		SNPS = extractPos(options, chr)

		# Filter: take minimum pvalue for each position
		SNPS_filt = filterSNPS(SNPS)

		# Store 
		with open(options['folder']['filtChrPlot'] + chr + '.csv', 'w') as outFile:
			writer = csv.writer(outFile)
			for pos in list(SNPS_filt.keys()):
				item = [pos, SNPS_filt[pos]]
				writer.writerow(item)

		
if __name__ == "__main__":
	chrArray = argsParser()
	main(chrArray)

