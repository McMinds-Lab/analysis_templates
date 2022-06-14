## three positional arguments specifying 1) the directory containing fastqs, 2) the output analysis directory, and 3) the file path of the config file
indir=$1
outdir=$2
config=$3

mkdir -p ${outdir}/02_DE-kupl-run/logs
cp ${config} ${outdir}/02_DE-kupl-run/config.json

cat <<EOF > ${outdir}/02_DE-kupl-run/02_DE-kupl-run.sbatch
#!/bin/bash
#SBATCH --qos=rra
#SBATCH --partition=rra
#SBATCH --time=7-00:00:00
#SBATCH --mem=187G
#SBATCH --job-name=02_DE-kupl-run
#SBATCH --output=${outdir}/02_DE-kupl-run/logs/02_DE-kupl-run_%a.log
#SBATCH --cpus-per-task=24

module purge
module load hub.apps/anaconda3
source activate dekupl

dekupl-run --configfile ${outdir}/02_DE-kupl-run/config.json -j\${SLURM_CPUS_PER_TASK} --resources ram=\${SLURM_MEM_PER_NODE} -p

EOF

if $autorun; then
    sbatch ${outdir}/02_DE-kupl-run/02_DE-kupl-run.sbatch
fi
