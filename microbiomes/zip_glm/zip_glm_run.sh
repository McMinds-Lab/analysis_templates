# get local variables
source local.env

asvs=$1
taxid_fp=$2
tags=$3
meta=$4
ids=$5
outdir=$6
K_s=$7
nchains=$8
opencl=$9
algorithm=${10}

mkdir -p ${outdir}/zip_glm/
cp zip_glm_run.r ${outdir}/zip_glm/
cp zip_glm.stan ${outdir}/zip_glm/

cat <<EOF > ${outdir}/zip_glm/zip_glm.sbatch
#!/bin/bash
#SBATCH --job-name=zip_glm
#SBATCH --partition=${partition}
#SBATCH --qos=${qos}
#SBATCH --mail-user=${email}
#SBATCH --mail-type=END,FAIL
#SBATCH --output=${outdir}/zip_glm/zip_glm.log
#SBATCH --ntasks=${nthreads}
#SBATCH --mem=${maxram}
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:1

module load hub.apps/R
Rscript ${outdir}/zip_glm/zip_glm_run.r ${asvs} ${taxid_fp} ${tags} ${meta} ${ids} ${outdir} ${K_s} ${nchains} ${nthreads} ${opencl} ${opencl_device} ${outdir}/zip_glm/ ${algorithm}

EOF

if $autorun; then
    sbatch ${outdir}/zip_glm/zip_glm.sbatch
fi
