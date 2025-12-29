#!/usr/bin/env python3
import sys
infile = sys.argv[1]

def readFastaAsDict(fastafile):
    fastaDict = {}
    allline = open(fastafile, 'r').read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

with open(infile, 'r') as r:
    fastadict = readFastaAsDict(infile)
    outfastadict = {}
    for sid, seq in fastadict.items():
        geneid = sid.split('.')[0]
        if geneid not in outfastadict or len(seq) > len(outfastadict[geneid]):
            outfastadict[geneid] = seq

for sid, seq in outfastadict.items():
    print(f'>{sid}\n{seq}')
