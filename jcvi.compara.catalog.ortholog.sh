#!/usr/bin/env sh
# require: jcvi, get_longest_transcript.py, seqkit, anchorfilter.py
# usage: <this script> <prefix1> <prefix2>
# prepare: {prefix}.gff3, {prefix}.cds.fasta, {prefix}.pep.fasta
# don't use . in prefix
for i in $1 $2; do python -m jcvi.formats.gff bed $i.gff3 -o $i.bed; done
for i in $1 $2; do python -m jcvi.formats.bed uniq $i.bed; done
for i in $1 $2; do get_longest_transcript.py $i.cds.fasta > $i.cdslong.fasta; done
for i in $1 $2; do get_longest_transcript.py $i.pep.fasta > $i.peplong.fasta; done
for i in $1 $2; do cut -f 4 $i.uniq.bed > $i.uniq.list; done
for i in $1 $2; do seqkit grep -f $i.uniq.list $i.cdslong.fasta | seqkit seq -i > $i.cds; done
for i in $1 $2; do seqkit grep -f $i.uniq.list $i.peplong.fasta | seqkit seq -i > $i.pep; done
for i in $1 $2; do mv $i.bed $i.unfiltered.bed && mv $i.uniq.bed $i.bed; done
python -m jcvi.compara.catalog ortholog --cpus 0 --no_dotplot $1 $2
anchorfilter.py $1.$2.lifted.anchors > $1.$2.pair