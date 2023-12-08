## ran:
/Users/Ryan/scripts/ont-guppy-cpu/bin/guppy_basecaller --trim_adapters --do_read_splitting --compress_fastq -i /Library/MinKNOW/data/Lambda_control_20231130/no_sample/20231130_1722_MN45077_FAX70185_75379f78/pod5 -s /Users/Ryan/outputs/nano -c /Users/Ryan/scripts/ont-guppy-cpu/data/dna_r10.4.1_e8.2_400bps_5khz_fast.cfg

## this should be implied by the choice of config: --flowcell FLO-MIN114 --kit SQK-LSK114 -m /Users/Ryan/scripts/ont-guppy-cpu/data/template_r10.4.1_e8.2_400bps_5khz_fast.jsn

## might want to add --min_qscore 0 because i can do that downstream if needed?

## better to use dorado and duplex basecalling...
