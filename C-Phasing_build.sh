#!/usr/bin/env bash
# mamba activate cphasing
set -euxo pipefail

prefix=$1
chrom=$2
ploidy=$3
cphasing utils assembly2agp $prefix.review.assembly -n $chrom:$ploidy 
cphasing agp2fasta $prefix.review.agp ../02.cphasing/$prefix.dup.tig.fasta > $prefix.review.asm.fasta
for i in $(seq 1 $ploidy); do grep "g$i" -A 1 --no-group-separator $prefix.review.asm.fasta > $prefix.hap$i.genome.fasta; done
for i in $(seq 1 $ploidy); do grep "g$i" --no-group-separator $prefix.review.corrected.agp  > $prefix.hap$i.agp; done
grep "utg" -A 1 --no-group-separator $prefix.review.asm.fasta > $prefix.unassigned.tig.fasta
cphasing plot -a $prefix.review.agp -m $prefix.10k.cool -o $prefix.genome.hic.png --add-hap-border --no-lines --only-chr