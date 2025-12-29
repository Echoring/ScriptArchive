#!/usr/bin/env python3
# convert count matrix to TPM matrix.
import argparse
import collections
parser = argparse.ArgumentParser()
parser.add_argument('-c', dest='input_count', required=True, help='(*Required) Input count matrix.')
parser.add_argument('-r', dest='ref_fasta', required=True, help='(*Required) all counted CDS sequence, FASTA format')
incount = parser.parse_args().input_count
inref = parser.parse_args().ref_fasta

def readFastaAsDict(fastafile):
    fastaDict = {}
    with open(fastafile, 'r') as inputfile:
        allline = inputfile.read()
        eachidseq = allline.split('>')
        for idseq in eachidseq:
            if idseq != '':
                sidraw, seqraw = idseq.split('\n', 1)
                sid = sidraw.split()[0].strip()
                seq = seqraw.replace('\n', '').upper()
                fastaDict[sid] = seq
    return fastaDict
cdsfasta = readFastaAsDict(inref)
cdslendict = {}
for sid, seq in cdsfasta.items():
    cdslendict[sid] = len(seq)

countdict = {}
with open(incount, 'r') as c:
    for line in c:
        if line.startswith('\t'):
            title = line.strip('\n')
        else:
            geneid = line.strip().split()[0]
            countdict[geneid] = line.strip().split()[1:]

rpkdict = collections.defaultdict(list)
totalrpklist = [0] * len(countdict[list(countdict.keys())[0]])
for geneid, countlist in countdict.items():
    for i in range(len(countdict[list(countdict.keys())[0]])):
        rpk = float(countlist[i]) / (cdslendict[geneid] / 1000)
        rpkdict[geneid].append(rpk)
        totalrpklist[i] += rpk
    
rpkmlist = [x / 1000000 for x in totalrpklist]

tpmdict = collections.defaultdict(list)
for geneid, rpklist in rpkdict.items():
    for i in range(len(rpklist)):
        tpmdict[geneid].append(str(round(rpklist[i] / rpkmlist[i], 2)))

print(title)
for geneid, tpmlist in tpmdict.items():
    tpm = '\t'.join(tpmlist)
    print(f'{geneid}\t{tpm}')