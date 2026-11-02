#!/usr/bin/env bash
set -euxo pipefail

genome=$1
protseq=$2 
rnadir=$3 
rnaseq=$4 
geneprefix=$5
minaalength=$6
threads=$7

# TE annotation: EDTA + LAI
mkdir $genome.EDTA && cd $genome.EDTA && ln -s ../$genome 
EDTA.pl --genome $genome --anno 1 --force 1 --threads $threads
LAI -t $threads -genome $genome -intact $genome.mod.EDTA.raw/LTR/$genome.mod.pass.list -all $genome.mod.EDTA.anno/$genome.mod.out
hard2softmask.py $genome $genome.mod.MAKER.masked > ../$genome.softmasked.fasta
cd ..

# RNA-seq alignment
mkdir $genome.RNAalign
hisat2-build -p $threads $genome $genome.RNAalign/$genome
IFS=',' read -ra id_array <<< "$rnaseq"
for id in "${id_array[@]}"; do
    if [ -f "$rnadir/${id}_1.fastq.gz" ] && [ -f "$rnadir/${id}_2.fastq.gz" ]; then
        hisat2 -x $genome.RNAalign/$genome -1 $rnadir/${id}_1.fastq.gz -2 $rnadir/${id}_2.fastq.gz --dta -p $threads 2>/dev/null | samtools sort -@ $threads -o $genome.RNAalign/$id.bam
    elif [ -f "$rnadir/${id}.fastq.gz" ]; then
        hisat2 -x $genome.RNAalign/$genome -U $rnadir/${id}.fastq.gz --dta -p $threads 2>/dev/null | samtools sort -@ $threads -o $genome.RNAalign/$id.bam
    fi
done

# gene structure annotation: BRAKER3 + PASA
braker.pl --genome $genome.softmasked.fasta --prot_seq $protseq --rnaseq_sets_dir $genome.RNAalign --rnaseq_sets_ids $rnaseq --species $genome --workingdir $genome.braker --gff3 --threads $threads --busco_lineage embryophyta_odb10
mkdir $genome.PASA && cd $genome.PASA && ln -s ../$genome
for bam in ../$genome.RNAalign/*.bam; do
    stringtie -o $genome.$(basename "$bam").stringtie.gtf -l ${$genome}_STRG --conservative -p $threads "$bam"
done
stringtie --merge -o $genome.stringtie.merged.gtf -l ${genome}_MSTRG $genome.*.stringtie.gtf
gffread -w $genome.stringtie.merged.mRNA.fasta -g "$genome" $genome.stringtie.merged.gtf
echo "DATABASE=$(pwd)/$genome.pasa_db.sqlite" >> $genome.pasa_db.config
seqclean $genome.stringtie.merged.mRNA.fasta
Launch_PASA_pipeline.pl -c $genome.pasa_db.config -g "$genome" -t $genome.stringtie.merged.mRNA.fasta.clean \
        -u $genome.stringtie.merged.mRNA.fasta --ALIGNERS blat,gmap,minimap2 --CPU $threads -C -R -T
Load_Current_Gene_Annotations.dbi -c $genome.pasa_db.config -g "$genome" -P "../$genome.braker/braker.gff3"
Launch_PASA_pipeline.pl -c $genome.pasa_db.config -g "$genome" -t $genome.stringtie.merged.mRNA.fasta.clean -A
cd ..
gffilter.py $genome $genome.PASA/*_PASA_updates.*.gff3 $genome.structure.anno --min_pep_length $minaalength --rename-expr "f'${geneprefix}c{geneobject.chrom[3:5]}g{str(genecount).zfill(5)}'" &>> $genome.structure.anno.log

# splice evaluation: Portcullis
gffread $genome.structure.anno.gff3 -T > $genome.structure.anno.gtf
portcullis full -t $threads -o ${genome}.portcullis $genome
junctools gtf markup $genome.structure.anno.gtf -j ${genome}.portcullis/3-filt/portcullis_filtered.pass.junctions.tab -o $genome.structure.anno.markup.gtf &>> $genome.structure.anno.junctools.log

# annotation completeness evaluation: OMArk
mkdir $genome.OMArk
omamer search --db ~/tools/OMArk/Viridiplantae.h5 --query $genome.structure.anno.peplong.fasta --out $genome.OMArk/$genome.structure.anno.peplong.omamer -t $threads
omark -f $genome.OMArk/$genome.structure.anno.peplong.omamer -d ~/tools/OMArk/Viridiplantae.h5 -o $genome.OMArk/
