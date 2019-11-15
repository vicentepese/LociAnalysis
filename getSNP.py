import numpy as np 
import json 
import csv
import os 
import sys 
import gzip
import argparse

def gene2PosRange(options):

    # Load protein data
    with open(options['file']['geneData'], 'r') as inFile:
        reader = csv.reader(inFile)
        geneSymbols = next(reader)[2:]

    # Open file 
    openFile = gzip.open(options['file']['GRCh37'], 'rt')

    genePos = list()
    i = 1
    for row in openFile:
        if '\t' not in row:
            continue
        if row.split('\t')[2] == 'gene':
            att = row.split('\t')[-1]
            gene = att.split(';')[5].split('=')[-1]
            print(gene)
            if any(gene in symbol for symbol in geneSymbols):
                init = row.split('\t')[3]
                end = row.split('\t')[4]
                if int(init) > options['IGHlow'] and int(end) < options['IGHhigh']:
                    genePos.append([gene, int(init), int(end)])

        i += 1

    # Save 
    genePos.insert(0, ['Gene', 'Init', 'End'])
    with open(options['file']['genePos'], 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerows(genePos)

    return genePos

def getCisSNP(options, genePos, chrArray):

    AA_OpenFile = gzip.open(options['folder']['Data'] + 'CHR_14.csv.gz')
    next(AA_OpenFile)

    # Genes in AA
    genes = [gene[0] for gene in genePos]

    # SNP in AA results
    cisSNPreg = {gene: [] for gene in genes}
    line = 0
    for snp in AA_OpenFile:
        print('Current line ' +str(line) + '\r')
    
        # If the SNP is not modulating a pair:
        if '&' not in str(snp).split(',')[1]:
                
            # Get pos
            if 'rs' in str(snp):
                pos = int(str(snp).split(',')[0].split(':')[-1])
                pos_range = [pos-options['posRange'], pos+options['posRange']]
            else:
                pos = int(str(snp).split(',')[0][2:])
                pos_range = [pos-options['posRange'], pos+options['posRange']]

            # If the position range is overlapping (at some extent) with the any gene range (check cis)
            if any((pos_range[1] > int(gene[1])) and (pos_range[0] < int(gene[2])) for gene in genePos):
            
                # Get index of genes it overlaps with
                idx = [idx for idx, gene in enumerate(genePos) if
                (pos_range[1] > int(gene[1])) and (pos_range[0] < int(gene[2]))]

                # Check if any the genes the SNP range is overlapping with is the gene it is regulating
                if any(genePos[i][0] in str(snp).split(',')[1] for i in idx):

                    # Get idx
                    idx = np.where(np.asarray([genePos[i][0] in str(snp).split(',')[1] for i in idx]) == True)

                    # If TRUE, the SNP is regulating the gene that is encoding 
                    item_app = str(snp).split(',')
                    item_app[0] = item_app[0][2:]
                    item_app[-1] = item_app[-1][:-3]
                    if 'rs' in item_app[0]:
                        pos = item_app[0].split(':')[-1]
                    else:
                        pos = item_app[0]
                    item_app.insert(1, pos)
                    item_app.insert(0, int(genePos[idx[0][0]][2]))
                    item_app.insert(0, int(genePos[idx[0][0]][1]))

                    cisSNPreg[genePos[idx[0][0]][0]].append(item_app)
        
        # Verbose
        line += 1

    saveTotalCisChr14 = cisSNPreg
    for gene in list(cisSNPreg.keys()):
        if len(cisSNPreg[gene]) > 0:
            cisSNPreg[gene].insert(0,['LowerPos','HighPos','snpsRS','SNP','geneSNP','statistic','pvalue','FDR','beta'])


    with open(options['folder']['Outputs'] + 'CHR_14_cis.csv', 'w') as outFile:
        writer = csv.writer(outFile)
        writer.writerow(['gene','LowerPos','HighPos','snpsRS','SNP','geneSNP','statistic','pvalue','FDR','beta'])
        for gene in list(saveTotalCisChr14.keys()):
            for snp in saveTotalCisChr14[gene][1:]:
                snp.insert(0, gene)
                writer.writerow(snp)

def main(chrArray):

    # Read options 
    with open('options.json','r') as inFile:
        options = json.load(inFile)

    # Get position range for each gene
    if 'genePos.csv' not in os.listdir(options['folder']['Resources']):
        genePos = gene2PosRange(options)
    else:
        with open(options['file']['genePos'],'r') as inFile:
            reader = csv.reader(inFile)
            next(reader)
            genePos = list()
            for row in reader:
                genePos.append(row)

    # Find SNP Cis-Regulator
    getCisSNP(options, genePos)



if __name__ == "__main__":
    main()