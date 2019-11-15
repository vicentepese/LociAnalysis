import numpy as np 
import json 
import csv
import os 
import sys 
import gzip
import argparse

def extractPos(options, chr):

    # Open file 
    OpenFile = gzip.open(options['folder']['Data'] + chr +'.csv.gz')

    # Create dict where key is a position
    SNPS = dict()
    line = 1
    for snp in OpenFile:
        print('Current line ' + str(line) + '\r')

        # Take position 
        pos = str(snp)

        line += 1

    return  

def main():

    # Read options
    with open('options.json','r') as inFile:
        options = json.load(inFile)

    # Extract positions 
    extractPos(options, 'CHR_14')
    
if __name__ == "__main__":
    main()

main()