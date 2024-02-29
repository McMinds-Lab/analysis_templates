# reverse complement a nucleotide sequence with IAUPAC degeneracy codes
echo ${seq} | tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd | rev

# example of piping between subshells that need different conda environments
# what makes the pipe 'skip' the second 'conda activate'? 
(conda activate samtools; samtools help) | (conda activate minimap2; cat)
