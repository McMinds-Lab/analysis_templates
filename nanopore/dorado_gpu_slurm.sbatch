#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --nodes=1
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --gres=gpu:1

module purge
module load apps/cuda/11.3.1
## auto batch size seems to lead to gpu memory issues on rra. unclear if 5mC is somehow worse than 5mCG
~/scripts/dorado-0.5.0-linux-x64/bin/dorado duplex --batchsize 64 --verbose sup,5mC_5hmC,6mA ~/data/Lambda_control_20231130/no_sample/20231130_1722_MN45077_FAX70185_75379f78/pod5/ > /shares/pi_mcmindsr/outputs/nano/duplex.bam

source activate SAMTOOLS
samtools fastq /shares/pi_mcmindsr/outputs/nano/duplex.bam | gzip --best > /shares/pi_mcmindsr/outputs/nano/duplex.fastq.gz