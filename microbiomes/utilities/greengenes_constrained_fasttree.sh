gg_reference_tree=./gg_13_8_otus/trees/97_otus_unannotated.tree
gg_reference_seqs=./gg_13_8_otus/rep_set/97_otus.fasta
outdir=./

mkdir ${outdir}/constrained_tree

vsearch --derep_fulllength ${outdir}/deblur_out/reference-hit.seqs.fa --output - |
vsearch --usearch_global - --db ${gg_reference_seqs} --id 0.97 --dbmatched ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta

filter_tree.py -i ${gg_reference_tree} -f ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta -o ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre

perl ./TreeToConstraints.pl < ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre > ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt

mafft --maxiterate 1000 --retree 3 <(cat ${outdir}/deblur_out/reference-hit.seqs.fa ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta) > ${outdir}/constrained_tree/aligned_ref_and_deblurred_seqs.fasta

FastTree -constraints ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt -nt < ${outdir}/constrained_tree/aligned_ref_and_deblurred_seqs.fasta > ${outdir}/constrained_tree/gg_constrained_fastttree.tre
