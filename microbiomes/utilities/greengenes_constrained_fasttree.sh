## location of reference tree and corresponding sequence. Could just as easily be tree from Silva or other database
gg_reference_tree=./gg_13_8_otus/trees/97_otus_unannotated.tree
gg_reference_seqs=./gg_13_8_otus/rep_set/97_otus.fasta
outdir=./

mkdir ${outdir}/constrained_tree

## using vsearch to filter reference to just seqs that are similar to our seqs (I think I did this to reduce the size of the constraints to make this more practical)
vsearch --derep_fulllength ${outdir}/deblur_out/reference-hit.seqs.fa --output - |
vsearch --usearch_global - --db ${gg_reference_seqs} --id 0.97 --dbmatched ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta

## I believe this was a QIIME script
filter_tree.py -i ${gg_reference_tree} -f ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta -o ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre

## create the contstraints file using perl script that is provided on FastTree website
perl ./TreeToConstraints.pl < ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre > ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt

mafft --maxiterate 1000 --retree 3 <(cat ${outdir}/deblur_out/reference-hit.seqs.fa ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta) > ${outdir}/constrained_tree/aligned_ref_and_deblurred_seqs.fasta

FastTree -constraints ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt -nt < ${outdir}/constrained_tree/aligned_ref_and_deblurred_seqs.fasta > ${outdir}/constrained_tree/gg_constrained_fastttree.tre
