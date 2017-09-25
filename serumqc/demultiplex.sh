#!/bin/sh
#SBATCH --mem=12G --time=00-02:00 -p qc -c 10 -J 'NGS_demultiplexing'
dir=$(pwd)
srun bcl2fastq --no-lane-splitting -r 10 -d 10 -p 10 -w 10 -o $dir --sample-sheet $1
srun /tools/bin/renamepl -v 's/((?<!L555))_R([12])_/_L555_R$2_/' $dir/*gz
