#!/bin/bash

set -ex

mail="nicolas.delhomme@umu.se"
proj=b2015052
in=/proj/$proj/nobackup/ConeDev/raw
out=/proj/$proj/nobackup/ConeDev
genome=/proj/b2011227/indices/STAR/Pabies01-genome
gff3=/proj/b2011227/reference/gff3/Eugene.gff3
start=2
end=2
mem=256
kallisto_fasta=/proj/b2011227/reference/fasta/Pabies1.0-all.phase.gff3.CDS.fa
kallisto_index=/proj/b2011227/indices/kallisto/Pabies1.0-all.phase.gff3.CDS.fa.inx


export SORTMERNADIR=/home/delhomme/sortmerna

module load bioinfo-tools FastQC trimmomatic samtools/1.3 star/2.4.0f1 htseq kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem \
-f $kallisto_fasta -K $kallisto_index \
-g $genome -G $gff3 -H $gff3 -t $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
