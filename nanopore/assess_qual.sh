conda activate samtools

samtools view -b -h -d dx:1 ~/outputs/nano/duplex.bam > ~/outputs/nano/onlyduplex.bam
samtools view -b -h -d dx:0 ~/outputs/nano/duplex.bam > ~/outputs/nano/onlysimplex.bam


conda activate nanoplot

NanoPlot --ubam ~/outputs/nano/onlyduplex.bam -o ~/outputs/nano/nanoplot_onlyduplex
NanoPlot --ubam ~/outputs/nano/onlysimplex.bam -o ~/outputs/nano/nanoplot_onlysimplex


## filter bam for duplex/simplex, convert it to fastq and gzip. samtools might have option to do all at once?

conda activate samtools
samtools view -h -d dx:1 ~/outputs/nano/duplex.bam | samtools fastq | gzip --best > ~/outputs/nano/onlyduplex.fastq.gz

conda activate minimap2
minimap2 -ax map-ont ~/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna ~/outputs/nano/onlyduplex.fastq.gz > ~/outputs/nano/onlyduplex_lambdamap.sam

conda activate samtools
samtools sort -O BAM ~/outputs/nano/onlyduplex_lambdamap.sam > ~/outputs/nano/onlyduplex_lambdamap.bam

conda activate nanoplot
NanoPlot --bam ~/outputs/nano/onlyduplex_lambdamap.bam -o ~/outputs/nano/nanoplot_onlyduplex_lambdamap



conda activate samtools
samtools view -h -d dx:0 ~/outputs/nano/duplex.bam | samtools fastq | gzip --best > ~/outputs/nano/onlysimplex.fastq.gz

conda activate minimap2
minimap2 -ax map-ont ~/outputs/lambda/lambda_genome/ncbi_dataset/data/genomic.fna ~/outputs/nano/onlysimplex.fastq.gz > ~/outputs/nano/onlysimplex_lambdamap.sam

conda activate samtools
samtools sort -O BAM ~/outputs/nano/onlysimplex_lambdamap.sam > ~/outputs/nano/onlysimplex_lambdamap.bam

conda activate nanoplot
NanoPlot --bam ~/outputs/nano/onlysimplex_lambdamap.bam -o ~/outputs/nano/nanoplot_onlysimplex_lambdamap


