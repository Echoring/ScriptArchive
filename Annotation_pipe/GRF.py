#!/usr/bin/env python3
# Usage: python GRF.py <genome file> <gff3 file> <min pep length> <outfile prefix> [--rename-expr RENAME_EXPR] [--skip-filters]
# rename Expr: a python expression to rename gene id, available variables: genecount, geneobject:
# e.g. "f'Avadeh1c{geneobject.chrom[3:5]}g{str(genecount).zfill(5)}'"
# if not provided, gene and transcript id will be kept as original.
import sys
import collections
import argparse
import contextlib

def readFastaAsDict(fastafile):
    fastaDict = {}
    with open(fastafile, 'r') as f:
        allline = f.read()
    eachidseq = allline.split('>')
    for idseq in eachidseq:
        if idseq != '':
            sidraw, seqraw = idseq.split('\n', 1)
            sid = sidraw.split()[0].strip()
            seq = seqraw.replace('\n', '').upper()
            fastaDict[sid] = seq
    return fastaDict

def reversedseq(seq: str):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    seq = seq[::-1]
    seq = ''.join(complement.get(base, base) for base in seq)
    return seq

def processgff3(gff3file, genomedict, minlen, skip_filters=False):
    class Gene:
        
        def __init__(self, record):
            self.id = None
            chrom, source, feature, start, end, score, strand, phase, attributes = record.strip().split('\t')
            self.attributes = []
            for attribute in attributes.split(';'):
                if attribute.startswith('ID='):
                    self.id = attribute.split('=')[1].strip()
                else:
                    self.attributes.append(attribute.strip())
            if not self.id:
                print(f'[Error] No ID attribute at line: {record}', file=sys.stderr)
                sys.exit(1)
            
            self.chrom = chrom
            self.source = source
            self.feature = feature
            self.start = int(start)
            self.end = int(end)
            self.score = score
            self.strand = strand
            self.phase = phase
            
            self.mrna = []
            self.seq = ''
            self.cdslong = '' 
            self.peplong = '' # end in empty if not a valid gene
            
        def getlongestCDS(self, skip_filters):
            if not self.mrna:
                if not skip_filters:
                    print(f'[Warning] Remove gene {self.id} for it has no mRNA feature', file=sys.stderr)
                    infodict[f'Removed gene: no mRNA feature'] += 1
                return ''
            if not skip_filters:
                self.mrna = [mrna for mrna in self.mrna if mrna.pep != '']
                if not self.mrna:
                    print(f'[Warning] Remove gene {self.id} for it has no valid mRNA', file=sys.stderr)
                    infodict[f'Removed gene: no valid mRNA'] += 1
                    return ''
            return sorted(self.mrna, key=lambda x: len(x.cds), reverse=True)[0]
        
        def updateseq(self, skip_filters):
            for mrna in self.mrna:
                mrna.updateseq(skip_filters)
            if self.strand == '+':
                self.seq = genomedict[self.chrom][self.start-1:self.end]
            else:
                self.seq = reversedseq(genomedict[self.chrom][self.start-1:self.end])
            mrnaobject = self.getlongestCDS(skip_filters)
            if mrnaobject == '':
                return
            self.cdslong = mrnaobject.cds
            self.peplong = mrnaobject.pep
            self.mrna.sort(key=lambda x: x.start)
        
        def updategeneid(self, newid):
            i = 1
            self.id = newid
            for mrna in self.mrna:
                mrna.transcriptid = f't{i}'
                mrna.updategeneid(self.id)
                i += 1
        
        def shortrecord(self):
            record = f'{self.chrom}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t' + \
            f'ID={self.id};{";".join(self.attributes)}\n'
            for mrna in self.mrna:
                record += mrna.shortrecord()
            return record

    class mRNA:
        
        def __init__(self, record):
            self.id = None
            self.parent = None
            chrom, source, feature, start, end, score, strand, phase, attributes = record.strip().split('\t')
            self.attributes = []
            for attribute in attributes.split(';'):
                if attribute.startswith('ID='):
                    self.id = attribute.split('=')[1].strip()
                elif attribute.startswith('Parent='):
                    self.parent = attribute.split('=')[1].strip()
                else:
                    self.attributes.append(attribute.strip())
            if not self.id:
                print(f'[Error] No ID attribute at line: {record}', file=sys.stderr)
                sys.exit(1)
            if not self.parent:
                print(f'[Error] No parent attribute at line: {record}', file=sys.stderr)
                sys.exit(1)   

            self.chrom = chrom
            self.source = source
            self.feature = feature
            self.start = int(start)
            self.end = int(end)
            self.score = score
            self.strand = strand
            self.phase = phase
            
            self.transcriptid = None
            self.child = []
            self.cds = ''
            self.pep = '' # end in empty if not a valid mRNA
        
        def joinCDS(self, skip_filters):
            cdschild = [child for child in self.child if child.feature == 'CDS']
            if not cdschild:
                if not skip_filters:
                    print(f'[Warning] Remove mRNA {self.id} for it has no valid CDS', file=sys.stderr)
                    infodict[f'Removed mRNA: no valid CDS'] += 1
                return ''
            if self.strand == '+':
                cdslist = sorted(cdschild, key=lambda x: x.start)
            else:
                cdslist = sorted(cdschild, key=lambda x: x.end, reverse=True)
            seq = ''
            chrom_seq = genomedict[self.chrom]
            for cds in cdslist:
                if self.strand == '+':
                    seq += chrom_seq[cds.start-1:cds.end]
                else:
                    seq += reversedseq(chrom_seq[cds.start-1:cds.end])
            if not skip_filters:
                if 'N' in seq:
                    print(f'[Warning] Remove mRNA {self.id} for N in CDS sequence: {seq}', file=sys.stderr)
                    infodict[f'Removed mRNA: N in CDS sequence'] += 1
                    return ''
                if len(seq) % 3 != 0:
                    print(f'[Warning] Remove mRNA {self.id} for CDS not divisible by 3: {seq}', file=sys.stderr)
                    infodict[f'Removed mRNA: CDS not divisible by 3'] += 1
                    return ''
                if seq.startswith('ATG') == False:
                    print(f'[Warning] Remove mRNA {self.id} for no start codon: {seq}', file=sys.stderr)
                    infodict[f'Removed mRNA: no start codon'] += 1
                    return ''
                if seq.endswith('TAA') == False and seq.endswith('TAG') == False and seq.endswith('TGA') == False:
                    print(f'[Warning] Remove mRNA {self.id} for no stop codon: {seq}', file=sys.stderr)
                    infodict[f'Removed mRNA: no stop codon'] += 1
                    return ''
            return seq
        
        def translate(self, seq, minlen, skip_filters):
            codon_table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
            'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
            }
            if seq == '':
                return ''
            protein = ''
            for i in range(0, len(seq), 3):
                codon = seq[i:i+3]
                if not skip_filters:
                    if codon not in codon_table:
                        print(f'[Warning] Remove mRNA {self.id} for invalid codon {codon} at position {i+1}: {seq}', file=sys.stderr)
                        infodict[f'Removed mRNA: invalid codon'] += 1
                        return ''
                    if codon == '' and i+3 < len(seq):
                        print(f'[Warning] Remove mRNA {self.id} for in-frame stop codon at position {i+1}: {seq}', file=sys.stderr)
                        infodict[f'Removed mRNA: in-frame stop codon'] += 1
                        return ''
                protein += codon_table.get(codon, 'X')  # if invalid, put X or something, but since skip, but to be safe
            if not skip_filters and len(protein) < minlen:
                print(f'[Warning] Remove mRNA {self.id} for product length < {minlen}: {protein}', file=sys.stderr)
                infodict[f'Removed mRNA: product length < {minlen}'] += 1
                return ''
            return protein
        
        def updateseq(self, skip_filters):
            self.child.sort(key=lambda x: x.start)
            self.cds = self.joinCDS(skip_filters)
            self.pep = self.translate(self.cds, minlen, skip_filters)
            
        def updategeneid(self, parentnewid):
            self.parent = parentnewid
            self.id = f'{self.parent}.{self.transcriptid}'
            for child in self.child:
                child.updateparentid(self.id)
            
        def shortrecord(self):
            record = f'{self.chrom}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t' + \
            f'ID={self.id};Parent={self.parent};{";".join(self.attributes)}\n'
            for child in self.child:
                record += child.shortrecord()
            return record
            
    class child:
        
        def __init__(self, record):
            self.id = None
            self.parent = None
            chrom, source, feature, start, end, score, strand, phase, attributes = record.strip().split('\t')
            self.attributes = []
            for attribute in attributes.split(';'):
                if attribute.startswith('ID='):
                    self.id = attribute.split('=')[1].strip()
                elif attribute.startswith('Parent='):
                    self.parent = attribute.split('=')[1].strip()
                else:
                    self.attributes.append(attribute.strip())
            if not self.id:
                print(f'[Error] No ID attribute at line: {record}', file=sys.stderr)
                sys.exit(1)
            if not self.parent:
                print(f'[Error] No parent attribute at line: {record}', file=sys.stderr)
                sys.exit(1)   
                
            self.chrom = chrom
            self.source = source
            self.feature = feature
            self.start = int(start)
            self.end = int(end)
            self.score = score
            self.strand = strand
            self.phase = phase
            
        def updateparentid(self, parentnewid):
            idx = [c for c in mRNAdict[self.parent].child if c.feature == self.feature].index(self)+1
            self.parent = parentnewid
            self.id = f'{self.parent}.{self.feature}{idx}'
            
        def shortrecord(self):
            record = f'{self.chrom}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t' + \
            f'ID={self.id};Parent={self.parent};{";".join(self.attributes)}\n'
            return record
    
    def priority(line):
        if line.startswith('#'):
            return 4
        if len(line.split('\t')) != 9:
            return 4
        feature = line.split('\t')[2]
        if feature == 'gene':
            return 1
        elif feature == 'mRNA':
            return 2
        else:
            return 3    

    genedict = {}
    mRNAdict = {}
    infodict = collections.defaultdict(int)
    with open(gff3file, 'r') as gff3:
        gff3lines = gff3.readlines()
    gff3lines.sort(key=priority)
    for line in gff3lines:
        linepriority = priority(line)
        if linepriority == 1:
            geneobject = Gene(line)
            genedict[geneobject.id] = geneobject
        elif linepriority == 2:
            mrnaobject = mRNA(line)
            genedict[mrnaobject.parent].mrna.append(mrnaobject)
            mRNAdict[mrnaobject.id] = mrnaobject
        elif linepriority == 3:
            childobject = child(line)
            mRNAdict[childobject.parent].child.append(childobject)
    infodict['Imported gene'] = len(genedict)
    infodict['Imported mRNA'] = len(mRNAdict)
    for geneobject in genedict.values():
        geneobject.updateseq(skip_filters)
    if skip_filters:
        geneobjects = list(genedict.values())
    else:
        geneobjects = [geneobject for geneobject in genedict.values() if geneobject.peplong != '']
    geneobjects.sort(key=lambda x: (x.chrom, x.start))
    return geneobjects, infodict

def main():
    parser = argparse.ArgumentParser(description='Filter and process GFF3 files with genome sequences.')
    parser.add_argument('genome_file', help='Path to the genome FASTA file')
    parser.add_argument('gff3_file', help='Path to the GFF3 file')
    parser.add_argument('outfile_prefix', help='Prefix for output files')
    parser.add_argument('--min_pep_length', type=int, default=0, help='Minimum peptide length')
    parser.add_argument('--rename-expr', default='', help='Python expression to rename gene IDs')
    parser.add_argument('--skip-filters', action='store_true', help='Skip all filters and export all genes')
    
    args = parser.parse_args()
    
    genomedict = readFastaAsDict(args.genome_file)
    geneobjects, infodict = processgff3(args.gff3_file, genomedict, args.min_pep_length, args.skip_filters)
    print(f'''[Info] {args.outfile_prefix}:
[Info] Imported gene: {infodict["Imported gene"]}
[Info] Imported mRNA: {infodict["Imported mRNA"]}
[Info] Removed gene: no mRNA feature: {infodict["Removed gene: no mRNA feature"]}
[Info] Removed gene: no valid mRNA: {infodict["Removed gene: no valid mRNA"]}
[Info] Removed mRNA: no valid CDS: {infodict["Removed mRNA: no valid CDS"]}
[Info] Removed mRNA: N in CDS sequence: {infodict["Removed mRNA: N in CDS sequence"]}
[Info] Removed mRNA: CDS not divisible by 3: {infodict["Removed mRNA: CDS not divisible by 3"]}
[Info] Removed mRNA: no start codon: {infodict["Removed mRNA: no start codon"]}
[Info] Removed mRNA: no stop codon: {infodict["Removed mRNA: no stop codon"]}
[Info] Removed mRNA: in-frame stop codon: {infodict["Removed mRNA: in-frame stop codon"]}
[Info] Removed mRNA: invalid codon: {infodict["Removed mRNA: invalid codon"]}
[Info] Removed mRNA: product length < {args.min_pep_length}: {infodict[f"Removed mRNA: product length < {args.min_pep_length}"]}
[Info] Retained gene: {len(geneobjects)}
[Info] Retained mRNA: {sum(len(geneobject.mrna) for geneobject in geneobjects)}''', file=sys.stdout)
    
    with contextlib.ExitStack() as stack:
        gff3out = stack.enter_context(open(f'{args.outfile_prefix}.gff3', 'w')) if not args.skip_filters else None
        cdsout = stack.enter_context(open(f'{args.outfile_prefix}.cds.fasta', 'w'))
        pepout = stack.enter_context(open(f'{args.outfile_prefix}.pep.fasta', 'w'))
        genomicout = stack.enter_context(open(f'{args.outfile_prefix}.genomic.fasta', 'w'))
        cdslongout = stack.enter_context(open(f'{args.outfile_prefix}.cdslong.fasta', 'w'))
        peplongout = stack.enter_context(open(f'{args.outfile_prefix}.peplong.fasta', 'w'))
        genecount = 1
        for geneobject in geneobjects:
            if args.rename_expr:
                geneobject.updategeneid(eval(args.rename_expr))
            if gff3out:
                gff3out.write(geneobject.shortrecord())
            for mrnaobject in geneobject.mrna:
                cdsout.write(f'>{mrnaobject.id}\n{mrnaobject.cds}\n')
                pepout.write(f'>{mrnaobject.id}\n{mrnaobject.pep}\n')
            genomicout.write(f'>{geneobject.id}\n{geneobject.seq}\n')
            cdslongout.write(f'>{geneobject.id}\n{geneobject.cdslong}\n')
            peplongout.write(f'>{geneobject.id}\n{geneobject.peplong}\n')
            genecount += 1
        
if __name__ == '__main__':
    main()
