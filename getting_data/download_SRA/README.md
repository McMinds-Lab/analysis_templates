download_SRA_project.sh is a script that simplifies the downloading of entire projects from NCBI's SRA using just the project ID. It automatically downloads all the sample (run) IDs, then gathers their metadata and dumps their fastq's, then gzips the data.

download_SRA_sample.sh bundles the downloading of fastqs and their metadata for a single sample

use:

configure slurm_params.sh for the slurm queue, parallelization, etc., and prepare_path.sh for loading packages and then:

download a project on a personal computer that needs the PATH set for entrez-direct and sratools:

`bash download_SRA_project.sh SRP324501 destination_directory false /path/to/prepare_path.sh`

download a project on a server that uses the SLURM scheduler and requires modules to be loaded:

`bash download_SRA_project.sh SRP324501 destination_directory /path/to/slurm_params.sh /path/to/prepare_path.sh`

download the sequences and metadata from a single sample, when the software is already in your PATH:

`bash download_SRA_sample.sh SRR5837082 SAMN07360402 destination_directory`

to keep a log of the script's output, try something like:

`mkdir -p destination_directory/logs'
`bash download_SRA_project.sh (options) | tee destination_directory/logs/download_SRA_project.log`
