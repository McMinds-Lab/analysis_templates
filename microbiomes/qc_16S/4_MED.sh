source ~/.bashrc
export PATH=/raid1/home/micro/mcmindsr/labhome/local/bin/anaconda2/bin:$PATH


outdir=/raid1/home/micro/mcmindsr/ryan/20170228_swarm_plus_MED/
medout=${outdir}/MED_swarms30


## for each swarm...
for file in ${outdir}/swarms_fastas/*.fasta; do

	fname=$(basename ${file})
	seqname=${fname/.*/}
	mkdir -p ${medout}/${seqname}

	vsearch --rereplicate ${file} --sizein --xsize --fasta_width 0 --output ${medout}/${seqname}/${seqname}_rereplicated.fasta

	## sequences must all have the same length
	o-pad-with-gaps ${medout}/${seqname}/${seqname}_rereplicated.fasta -o ${medout}/${seqname}/${seqname}_rereplicated_padded.fasta

	## apply MED. Allow lots of variation in each node so no reads will be discarded
	decompose -T -S -V 100 -R --skip-check-input-file ${medout}/${seqname}/${seqname}_rereplicated_padded.fasta -o ${medout}/${seqname}/MED

	seqlist="${seqname}"

	## for each MED node...
	for node in ${medout}/${seqname}/MED/NODES/*.unique; do

		## print one line that starts with the name of the node (renamed to match the name of the node's most abundant sequence) and is followed by each dereplicated sequence name that's in the node (append this within the loop, so that this single file will wind up with every MED node from every Swarm)
		awk 'NR==2 {line=$1} NR>2 {line=line" "$1} END {print line}' RS='>' FS='|' ${node} >> ${outdir}/derep_to_MEDs.txt

		## append the name of the MED node to a variable that contains all the other MED names within the Swarm
		medname=$(awk 'NR==2 {print $1; exit}' RS='>' FS='|' ${node})

		if [[ ! ${seqlist} =~ ${medname} ]]
		then
			seqlist="${seqlist} ${medname}"
		fi

    done

	## print one line that starts with the name of the swarm and is followed by each MED node that is in that swarm
	echo ${seqlist} >> ${outdir}/MEDs_to_swarms.txt


	## clean up
	rm ${medout}/${seqname}/${seqname}_rereplicated.fasta ${medout}/${seqname}/${seqname}_rereplicated_padded.fasta ${medout}/${seqname}/MED/NODES/*.fa ${medout}/${seqname}/MED/MATRIX-COUNT.txt ${medout}/${seqname}/MED/MATRIX-PERCENT.txt ${medout}/${seqname}/MED/READ-DISTRIBUTION.txt
	rm -r ${medout}/${seqname}/MED/FIGURES ${medout}/${seqname}/MED/HTML-OUTPUT ${medout}/${seqname}/MED/OUTLIERS

done



