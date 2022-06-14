## two positional arguments specifying 1) the file path of the config file, 2) the output analysis directory
config=$1
outdir=$2

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
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1 

module purge
module load hub.apps/anaconda3
source activate dekupl

dekupl-run --configfile ${outdir}/02_DE-kupl-run/config.json -j\${SLURM_NTASKS} --resources ram=\${SLURM_MEM_PER_NODE} -p

EOF

if $autorun; then
    sbatch ${outdir}/02_DE-kupl-run/02_DE-kupl-run.sbatch
fi
