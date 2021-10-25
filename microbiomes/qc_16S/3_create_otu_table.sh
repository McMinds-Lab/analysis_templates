## read first file and make a dictionary where the keys correspond to unique sequence hashes and the values correspond to swarm ids
## read second file to make dictionary that translates original names to hash-names
## read third file that maps original sequences to unique sequences; translate each original sequence's corresponding unique sequence to the hash name and look that up in the swarm dictionary to find the swarm that it belongs to

awk 'BEGIN {nfile=0} \
FNR==1 {nfile+=1} \
nfile==1 && $1~"H" {swarm[$9]=$10} \
nfile==1 && $1~"S" {swarm[$9]=$9} \
nfile==2 {derep[$2]=$1} \
nfile==3 && $1~"H" {split($9,out1,"_"); split(swarm[derep[$10]],out2,";"); print out1[1],out2[1]} \
nfile==3 && $1~"S" {split($9,out1,"_"); split(swarm[derep[$9]],out2,";"); print out1[1],out2[1]}' \
RS='\\n' FS='\\s' /raid1/home/micro/mcmindsr/labhome/ryan/20170228_swarm_plus_MED/swarm/swarms.uc \
RS='>' FS='\\s' /raid1/home/micro/mcmindsr/labhome/ryan/20170228_swarm_plus_MED/seqs_derep.fasta \
RS='\\n' FS='\\s' /raid1/home/micro/mcmindsr/labhome/ryan/20170228_swarm_plus_MED/seqs_derep.uc \
| sort | uniq -c > /raid1/home/micro/mcmindsr/labhome/ryan/20170228_swarm_plus_MED/swarms_summed_long.txt




