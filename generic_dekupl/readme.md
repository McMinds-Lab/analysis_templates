this workflow requires a version of DE-kupl that has been modified. We accomplished this by:

1) installing DE-kupl via conda
2) adding https://github.com/rmcminds/dekupl-run/blob/master/bin/zeroinfl_diff_method.R and https://github.com/rmcminds/dekupl-run/blob/master/bin/binomial_diff_method.R to the dekupl bin directory in the conda environment
3) replacing the standard snakefile in the conda env with https://github.com/rmcminds/dekupl-run/blob/master/Snakefile
4) installing the package 'pscl' in the R installation included in the conda env. (e.g. conda activate dekupl; start an interactive R session; run install.packages('pscl'))
5) installing the packages 'Rmpi' and 'doMPI' (our cluster requires gcc and mpi modules to first be loaded, a variable set, and the installation to occur where mpi can be loaded eg. using srun rather than directly on the head node): module load compilers/gcc/4.8.1; module load mpi/openmpi/1.6.1; export OMPI_MCA_mtl=^psm; srun --pty R --interactive; install.packages(c('Rmpi','doMPI'), configure.args=list(Rmpi=c("--with-mpi=/apps/openmpi/1.6.1/","--with-Rmpi-include=/apps/openmpi/1.6.1/include","--with-Rmpi-libpath=/apps/openmpi/1.6.1/lib","--with-Rmpi-type=OPENMPI")), INSTALL_opts="--no-test-load", dependencies=TRUE, repos = "http://cran.us.r-project.org")
6) The mpi modules must also be loaded when running DE-kupl. Note that each chunk seems to use about 4 cores (maybe pscl::zeroinfl() itself is threaded?) so when deciding on the number of nodes and the number of jobs per node with e.g. slurm, you don't want to overcrowd individual nodes.)

These changes allow the use of a zero inflated negative binomial test, or a binary presence/absence test, by using the key diff_method: zeroinfl or diff_method: binomial

For 'zeroinfl', a gene is considered 'significant' if it's either differentially prevalent OR differentially abundant when present. the 'log2fc' value is not actually a log2-fold-change when the gene is differentially prevalent, so units are mixed up among genes and as a whole are essentially meaningless. Therefore it probably only makes sense to remove an effect size cutoff with the key log2fc_threshold:0 in the config.json.  

These scripts appear to be MUCH slower than the DESeq2 option, and this is why I made the modifications to use MPI - this is the most convenient way to parallelize across more than one machine. 
