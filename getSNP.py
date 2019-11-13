import numpy as np 
import json 
import csv
import os 
import sys 
import gzip

def importGeneSymbols(options):
    
    with open(options['file']['geneData'], 'r') as inFile:
        reader = csv.reader(inFile)
        geneSymbols = next(reader)[2:]
    
    return geneSymbols

def findSNPs(options, SNP_range):

    # Open CHR file 
    CHR_file = gzip.open(options['file']['CHR14Data'])
    next(CHR_file)

    snp_pos = list()
    for snp in CHR_file:
        if 'rs' in str(snp).split(',')[0]:
            if int(str(snp).split(',')[0].split(':')[-1]) > SNP_range[0] and int(str(snp).split(',')[0].split(':')[-1]) < SNP_range[1]:
                app_item = str(snp).split(',')
                app_item[0] = app_item[0][2:]
                app_item[-1] = app_item[-1][:-3]
                snp_pos.append(app_item)
        else:
            if int(str(snp).split(',')[0][2:]) > SNP_range[0] and int(str(snp).split(',')[0][2:]) < SNP_range[1]:
                app_item = str(snp).split(',')
                app_item[0] = app_item[0][2:]
                app_item[-1] = app_item[-1][:-3]
                snp_pos.append(app_item)
    return snp_pos



def main():

    # Read options 
    with open('options.json','r') as jsonFile:
        options = json.load(jsonFile)

    # Import geneData
    geneSymbols = importGeneSymbols(options)

    # Get postions (NEED DATASET)
    SNP_range = np.array([105863196, 105863260])

    # Add 500kb
    SNP_range = SNP_range + np.array([-500e3, 500e3])

    # Find SNPS
    snp_pos = findSNPs(options, SNP_range)

    # Save
    with open('tst.csv', 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerows(snp_pos)
    


if __name__ == "__main__":
    pass

main()