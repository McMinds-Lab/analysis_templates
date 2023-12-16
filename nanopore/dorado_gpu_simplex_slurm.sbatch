#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --nodes=1
#SBATCH --mem=175G
#SBATCH --time=7-00:00:00
#SBATCH --gres=gpu:1

module purge
module load apps/cuda/11.3.1
## duplex is not ideal for amplicons of any size, because it is hard to know if a chimera is due to the second strand following, or a completely different amplicon
## read splitting with dorado is not customizable and not perfect. downstream tools should be able to take care of chimeras and be careful about demultiplexing settings
## auto batch size seems to lead to gpu memory issues on rra.
~/scripts/dorado-0.5.0-linux-x64/bin/dorado basecaller --batchsize 64 --trim adapters --verbose sup /shares/pi_mcmindsr/raw_data/20231212_1932_MN45077_FAX70185_835ac3e3/pod5/ > /shares/pi_mcmindsr/outputs/20231212_1932_MN45077_FAX70185_835ac3e3_bc20231216/20231216_basecalls.bam
