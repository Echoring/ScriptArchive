#!/usr/bin/env python3
# Usage infile > file
import sys
infile = sys.argv[1]

genescoreAB = {}
genescoreBA = {}
with open(infile, 'r') as r:
    for line in r:
        if not line.startswith('#'):
            geneA = line.split()[0]
            geneB = line.split()[1]
            score = int(line.split()[2].strip('L'))
            if geneA not in genescoreAB or genescoreAB[geneA][1] < score:
                genescoreAB[geneA] = [geneB, score]
            if geneB not in genescoreBA or genescoreBA[geneB][1] < score:
                genescoreBA[geneB] = [geneA, score]
                
for geneA in genescoreAB:
    geneB = genescoreAB[geneA][0]
    if genescoreBA[geneB][0] == geneA:
        print(f'{geneA}\t{genescoreAB[geneA][0]}')