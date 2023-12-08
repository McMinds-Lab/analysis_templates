#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --nodes=1
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --gres=gpu:1

module purge
module load apps/cuda/11.3.1
~/scripts/dorado-0.5.0-linux-x64/bin/dorado duplex --batchsize 64 sup,5mCG_5hmCG,6mA ~/data/Lambda_control_20231130/no_sample/20231130_1722_MN45077_FAX70185_75379f78/pod5/ > /shares/pi_mcmindsr/outputs/nano/duplex.bam

