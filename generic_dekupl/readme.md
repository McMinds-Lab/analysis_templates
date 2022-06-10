this workflow requires a version of DE-kupl that has been modified. We accomplished this by:

1) installing DE-kupl via conda
2) adding https://github.com/rmcminds/dekupl-run/blob/master/bin/zeroinfl_diff_method.R and https://github.com/rmcminds/dekupl-run/blob/master/bin/binomial_diff_method.R to the dekupl bin directory in the conda environment
3) replacing the standard snakefile in the conda env with https://github.com/rmcminds/dekupl-run/blob/master/Snakefile
4) installing the package 'pscl' in the R installation included in the conda env. (e.g. conda activate dekupl; start an interactive R session; run install.packages('pscl'))

These changes allow the use of a zero inflated negative binomial test, or a binary presence/absence test, by using the key diff_method: zeroinfl or diff_method: binomial

For 'zeroinfl', a gene is considered 'significant' if it's either differentially prevalent OR differentially abundant when present. the 'log2fc' value is not actually a log2-fold-change when the gene is differentially prevalent, so units are mixed up among genes and as a whole are essentially meaningless. Therefore it probably only makes sense to remove an effect size cutoff with the key log2fc_threshold:0 in the config.json.  

These scripts appear to be MUCH slower than the DESeq2 option, but I am working on it...
