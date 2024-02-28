# reverse complement a nucleotide sequence with IAUPAC degeneracy codes
echo ${seq} | tr ACGTRYSWKMBVDHacgtryswkmbvdh TGCAYRSWMKVBHDtgcayrswmkvbhd | rev
