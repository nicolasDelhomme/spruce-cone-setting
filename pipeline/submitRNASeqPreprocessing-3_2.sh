#!/bin/bash

set -ex

#mail="veronika.nordal@slu.se"
mail="nicolas.delhomme@slu.se"
proj=b2015052
in=/proj/b2015052/nobackup/ConeDev-3_2/raw
out=/proj/b2015052/nobackup/ConeDev-3_2
genome=/proj/b2011227/indices/STAR/v2.5.2b/Pabies01-genome
gtf=/proj/b2011227/reference/gff3/Eugene.gtf
gff3=/proj/b2011227/reference/gff3/Eugene.gff3
start=2
end=2
mem=256
kallisto_fasta=/proj/b2011227/reference/fasta/Pabies1.0-all.phase.gff3.CDS.fa
kallisto_index=/proj/b2011227/indices/kallisto/Pabies1.0-all.phase.gff3.CDS.fa.inx

module load bioinfo-tools FastQC SortMeRNA trimmomatic samtools/1.3 star/2.5.2b htseq kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem \
-f $kallisto_fasta -K $kallisto_index \
-g $genome -G $gtf -H $gff3 -t $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
