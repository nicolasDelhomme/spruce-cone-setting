#!/bin/bash

set -ex

proj=b2015052
mail="veronika.nordal@slu.se"
in=/proj/$proj/nobackup/ConeDev-extension
out=/proj/$proj/nobackup/ConeDev-extension/multiqc

if [ ! -d $out ]; then
	mkdir -p $out
fi

module load bioinfo-tools MultiQC

sbatch --mail-user=$mail -o $in/multiqc.out -e $in/multiqc.err -A $proj $UPSCb/pipeline/runMultiQC.sh $in $out
