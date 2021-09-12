
for run in "$@"; do
	echo $run
	# downloands run data
	gsutil -m cp gs://covid_terra/${run}/${run}_terra_data_table.tsv .
	while read name gisaid proj; do
		grep $name ${run}_terra_data_table.tsv >> ${run}_seqs_to_download.tsv
	done < fasta_rename_accession_to_gisaid_id.tsv

	while read p; do
	  echo $p
	  accession=$(echo $p | awk '{print $1}')
	  echo $accession
	  barcode=$(echo $p | awk '{print $2}')
	  path=$(echo $p | awk '{print $3}')
	  echo $barcode
	 	 if gsutil ls ${path}/*.fastq.gz 1> /dev/null 2>&1; then
	 	 	gsutil cat ${path}/*.fastq.gz > CO-CDPHE-${accession}.fastq.gz
	 	 else
	 	 	gsutil cat ${path}/*.fastq > CO-CDPHE-${accession}.fastq
	 	 	gzip CO-CDPHE-${accession}.fastq
	 	 fi
	done < ${run}_seqs_to_download.tsv

done