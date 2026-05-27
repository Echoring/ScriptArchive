#!/usr/bin/env python3
# Usage: ... <ingfa> (output as ingfa(-last suffix).fasta)
# Convert gfa format to fasta format.
import sys
inputgfa = sys.argv[1]

def changeSuffix(filename, newsuffix):
    namelist = filename.split('.')
    namelist[-1] = newsuffix
    return '.'.join(namelist)


with open(inputgfa, 'r') as gfa, open(changeSuffix(inputgfa, 'fasta'), 'w') as fasta:
    for line in gfa:
        flag = line.split()[0].strip()
        seqid = line.split()[1].strip()
        seq = line.split()[2].strip()
        if flag == 'S':
            fasta.write(f'>{seqid}\n{seq}\n')
