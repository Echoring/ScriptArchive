#!/usr/bin/env sh
# Software: jcvi, seqkit, ParaAT, mafft, ParaFly, KaKs_Calculator
# Script: orthologfilter.py
# Usage: <this script> <gffA> <cdsA> <pepA> <gffB> <cdsB> <pepB> <threads> <outprefix> <mode>
# Calculate Ks between all ortholog of two genomes
# cds and pep file: each gene gives only one transcript.
# mode: single=ignore gene with multiple ortholog; RBH=only calculate Reciprocal Best Hit ortholog pair; all=calculate all pairs, allow 1-to-many
# mamba activate jcvi
set -euxo pipefail
Agff3=$(readlink -f $1)
Acds=$(readlink -f $2)
Apep=$(readlink -f $3)
Bgff3=$(readlink -f $4)
Bcds=$(readlink -f $5)
Bpep=$(readlink -f $6)
mkdir $8.tmp
cd $8.tmp
ln -s $Agff3 A.gff3
ln -s $Acds A.cds.fasta
ln -s $Apep A.pep.fasta
ln -s $Bgff3 B.gff3
ln -s $Bcds B.cds.fasta
ln -s $Bpep B.pep.fasta
for i in A B; do python -m jcvi.formats.gff bed $i.gff3 -o $i.bed; done
for i in A B; do python -m jcvi.formats.bed uniq $i.bed; done
for i in A B; do cut -f 4 $i.uniq.bed > $i.uniq.list; done
for i in A B; do seqkit grep -r -f $i.uniq.list $i.cds.fasta | seqkit seq -i > $i.cds; done
for i in A B; do seqkit grep -r -f $i.uniq.list $i.pep.fasta | seqkit seq -i > $i.pep; done
for i in A B; do mv -f $i.uniq.bed $i.bed; done
python -m jcvi.compara.catalog ortholog --self_remove 100 --cpus $7 --no_dotplot A B
python -m jcvi.compara.catalog ortholog --self_remove 100 --cpus $7 --no_dotplot B A
python3 ~/script/orthologfilter.py A.B.lifted.anchors B.A.lifted.anchors $9 > ../$8.pair.txt
cat A.cds.fasta B.cds.fasta > A.B.cds.fasta
cat A.pep.fasta B.pep.fasta > A.B.pep.fasta
echo "$7" > proc
ParaAT.pl -h ../$8.pair.txt -n A.B.cds.fasta -a A.B.pep.fasta -m mafft -p proc -f axt -o ./tmp
for i in `find ./tmp -name "*.axt"`;do echo "KaKs_Calculator -i $i -o ${i}.kaks -m YN" >> kaks.fly; done
ParaFly -c kaks.fly -CPU $7
for i in `find ./tmp -name "*.kaks"`;do awk 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' $i >> tmp/all-kaks.txt; done
sort ./tmp/all-kaks.txt | uniq | sed '1i\Seq\tKa\tKs\tKa/Ks' > ../$8.kaks.txt
rm -rf ./tmp