import numpy as np 
import json 
import csv
import os 
import sys 
import gzip
import argparse

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
		if isinstance(options['chrArray']['CHRpreprocess'], int):
			chrArray = ["CHR_" + str(options)]
		else: 
			chrArray = ["CHR_" + str(s) for s in options['chrArray']['CHRpreprocess']]

	# Else, analyze the inputed chromosome (slurm)
	elif len(args.ChrIndex) > 1:
		chrArray = ["CHR_" + str(args.ChrIndex)]
	elif len(args.ChrIndex) == 1:
		chrArray = ["CHR_" + str(args.ChrIndex)]

	return chrArray

def extractPos(options, chr):

    # Open file 
    OpenFile = gzip.open(options['folder']['Data'] + chr +'.csv.gz')
    next(OpenFile)

    # Create dict where key is a position
    SNPS = dict()
    line = 1
    print("Processing file")
    for snp in OpenFile:

        if line % 100000 == 0:
            print('Current line: ' +  str(line))

        # Filter by pvalue (continue if bigger than pvalue threshold)
        if float(str(snp).split(',')[3]) > options['getTrans']['pval_thres'] or '&' in str(snp).split(',')[1]:
            continue

        # Take position 
        if 'rs' in str(snp).split(',')[0]: 
            pos = str(snp).split(',')[0].split(':')[-1]
        else:
            pos = str(snp).split(',')[0][2:]

        # Append to dict
        if pos in list(SNPS.keys()):
            SNPS[pos].append([str(snp).split(',')[1], str(snp).split(',')[3]])
        else:
            SNPS[pos] = [[str(snp).split(',')[1], str(snp).split(',')[3]]]
        line += 1

    return SNPS 

def filterSNPS(SNPS):

    # For each position 
    for pos in list(SNPS.keys()):

        # Get pvalues
        pvals = np.asarray([float(gene[1]) for gene in SNPS[pos]])
        
        # Get index of the minimum pvalue
        pval_idx = np.where(pvals == min(pvals))

        # Keep only the gene with minimum pvalue
        SNPS[pos] = SNPS[pos][pval_idx[0][0]]

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

        # Store the array 
        SNPS_lists = list()
        for pos in list(SNPS_filt.keys()):
            item_app = [pos]
            item_app.extend(SNPS_filt[pos])
            SNPS_lists.append(item_app)

        with open(options['folder']['chrPlot'] + chr + '.csv', 'w') as outFile:
            writer = csv.writer(outFile)
            writer.writerows(SNPS_lists)

        
if __name__ == "__main__":
    chrArray = argsParser()
    main(chrArray)

