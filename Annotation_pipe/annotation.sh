#!/usr/bin/env bash
set -euxo pipefail

genome=$1
protseq=$2
rnadir=$3
rnaseq=$4 
geneprefix=$5
minaalength=$6
threads=$7
omamerdb=$8
emapperdbpath=$9

# TE annotation: EDTA + LAI
mkdir $genome.EDTA && cd $genome.EDTA && ln -s ../$genome 
EDTA.pl --genome $genome --anno 1 --force 1 --threads $threads
LAI -t $threads -genome $genome -intact $genome.mod.EDTA.raw/LTR/$genome.mod.pass.list -all $genome.mod.EDTA.anno/$genome.mod.out
python3 ~/script/hard2softmask.py $genome $genome.mod.MAKER.masked > ../$genome.softmasked.fasta # require numpy + Biopython
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

# gene structure annotation: BRAKER3 + ICU
braker.pl --genome $genome.softmasked.fasta --prot_seq $protseq --rnaseq_sets_dir $genome.RNAalign --rnaseq_sets_ids $rnaseq --species $genome --workingdir $genome.braker --gff3 --threads $threads --busco_lineage embryophyta_odb10
mkdir $genome.ICU
ICU.py pipe $genome $genome.braker/braker.gff3 $genome.RNAalign/*.bam --prefix $genome.ICU/$genome.ICU -p $threads # require Portcullis, gffread, StringTie
# change rename-expr if chromosome is not named ChrXX or have non-chromosome annotation. e.g. "f'Achred5c{geneobject.chrom[3:5].zfill(2) if geneobject.chrom[0:3] == \"chr\" else \"00\"}g{str(genecount).zfill(5)}'", "f'Avadeh1c{geneobject.chrom[7:9]}g{str(genecount).zfill(5)}'"
python3 ~/script/GRF.py $genome $genome.ICU/$genome.ICU.integrated.gff3 $genome.structure.anno --min_pep_length $minaalength --rename-expr "f'${geneprefix}c{geneobject.chrom[3:5]}g{str(genecount).zfill(5)}'" &>> $genome.structure.anno.log 

# intron precise evaluation: Portcullis
gffread $genome.structure.anno.gff3 -T > $genome.structure.anno.gtf
junctools gtf markup $genome.structure.anno.gtf -j $genome.ICU/$genome.ICU_portcullis_out/3-filt/portcullis_filtered.pass.junctions.tab -o $genome.structure.anno.markup.gtf &>> $genome.structure.anno.junctools.log

# annotation completeness evaluation: OMArk
mkdir $genome.OMArk
omamer search --db $omamerdb --query $genome.structure.anno.peplong.fasta --out $genome.OMArk/$genome.structure.anno.peplong.omamer -t $threads
omark -f $genome.OMArk/$genome.structure.anno.peplong.omamer -d $omamerdb -o $genome.OMArk/
