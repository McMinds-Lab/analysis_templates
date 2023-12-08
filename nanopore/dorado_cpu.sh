/Users/rm/scripts/ont-dorado-server/bin/dorado_basecall_server --trim_adapters --do_read_splitting --log_path /Users/rm/outputs/nano/server_logs --config dna_r10.4.1_e8.2_400bps_5khz_duplex_hac.cfg -p 5555

/Users/rm/scripts/ont-dorado-server/bin/ont_basecall_client --trim_adapters --do_read_splitting --compress_fastq --config dna_r10.4.1_e8.2_400bps_5khz_duplex_hac.cfg --input_path /Library/MinKNOW/data/Lambda_control_20231130/no_sample/20231130_1722_MN45077_FAX70185_75379f78/pod5 --save_path /Users/rm/outputs/nano -p 5555

