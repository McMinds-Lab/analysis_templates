module load hub.apps/anaconda3
conda activate sepp


run_sepp.py -t mock/pyrg/sate.tre -r mock/pyrg/sate.tre.RAxML_info -a mock/pyrg/sate.fasta -f mock/pyrg/pyrg.even.fas
