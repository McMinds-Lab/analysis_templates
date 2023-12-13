conda activate samtools

samtools view -b -h -d dx:1 ~/outputs/nano2/duplex.bam > ~/outputs/nano2/onlyduplex.bam
samtools view -b -h -d dx:0 ~/outputs/nano2/duplex.bam > ~/outputs/nano2/onlysimplex.bam


conda activate nanoplot

NanoPlot --ubam ~/outputs/nano2/onlyduplex.bam -o ~/outputs/nano2/nanoplot_onlyduplex
NanoPlot --ubam ~/outputs/nano2/onlysimplex.bam -o ~/outputs/nano2/nanoplot_onlysimplex


## filter bam for duplex/simplex, delete the last line because the file is incomplete and the last line might be nonsense, convert it to fastq and gzip

conda activate samtools
samtools view -h -d dx:1 ~/outputs/nano2/duplex.bam | sed -e '$ d' | samtools fastq | gzip --best > ~/outputs/nano2/onlyduplex.fastq.gz

conda activate minimap2
minimap2 -ax map-ont /Users/Ryan/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna ~/outputs/nano2/onlyduplex.fastq.gz > ~/outputs/nano2/onlyduplex_lambdamap.sam

conda activate samtools
samtools sort -O BAM ~/outputs/nano2/onlyduplex_lambdamap.sam > ~/outputs/nano2/onlyduplex_lambdamap.bam

conda activate nanoplot
NanoPlot --bam ~/outputs/nano2/onlyduplex_lambdamap.bam -o ~/outputs/nano2/nanoplot_onlyduplex_lambdamap



conda activate samtools
samtools view -h -d dx:0 ~/outputs/nano2/duplex.bam | sed -e '$ d' | samtools fastq | gzip --best > ~/outputs/nano2/onlysimplex.fastq.gz

conda activate minimap2
minimap2 -ax map-ont /Users/Ryan/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna ~/outputs/nano2/onlysimplex.fastq.gz > ~/outputs/nano2/onlysimplex_lambdamap.sam

conda activate samtools
samtools sort -O BAM ~/outputs/nano2/onlysimplex_lambdamap.sam > ~/outputs/nano2/onlysimplex_lambdamap.bam

conda activate nanoplot
NanoPlot --bam ~/outputs/nano2/onlysimplex_lambdamap.bam -o ~/outputs/nano2/nanoplot_onlysimplex_lambdamap


