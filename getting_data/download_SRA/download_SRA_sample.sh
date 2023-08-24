## three or four positional arguments specifying 1) a sample/run identifier starting with SRR, 2) a 'BioSample' identifier starting with SAMN, 3) the output directory, and 4) optional path to a script that sets your PATH variable to include the entrez-direct and sratools packages
sample=$1
biosample=$2
outdir=$3
prepare_path=$4

## if tools are not automatically in path, put them there
if [ $# -eq 4 ]; then
  source ${prepare_path}
fi

mkdir -p ${outdir}/${sample}

## download sample metadata
esearch -db biosample -query ${biosample} | efetch > ${outdir}/${sample}/metadata.txt
esearch -db biosample -query ${biosample} | efetch -format xml > ${outdir}/${sample}/metadata.xml

## download raw sequences from sample
fasterq-dump -e 1 ${sample} -t ${outdir}/${sample} -O ${outdir}/${sample}

## compress sequence file
for file in ${outdir}/${sample}/*.fastq; do
  gzip --best ${file} &
done
wait
