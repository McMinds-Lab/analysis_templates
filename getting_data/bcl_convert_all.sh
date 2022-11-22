## two positional arguments specifying 1) an input basecall directory and 2) the output directory

indir=$1
outdir=$2

mkdir -p ${outdir}

cat <<EOF > ${outdir}/bcl_convert_all.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=6-00:00:00
#SBATCH --mem=60
#SBATCH --job-name=bcl_convert_all
#SBATCH --output=${outdir}/bcl_convert_all.log

/shares/omicshub/apps/bcl-convert/bcl-convert --version

/shares/omicshub/apps/bcl-convert/bcl-convert --no-sample-sheet true --bcl-input-directory ${indir} --output-directory ${outdir}

EOF

if $autorun; then

sbatch ${outdir}/bcl_convert_all.sbatch

fi

