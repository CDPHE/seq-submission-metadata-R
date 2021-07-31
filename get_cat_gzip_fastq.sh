
for run in "$@"; do
	echo $run
	#mkdir -p ${run}_fastq
	#cd ${run}_fastq
	# downloands run data
	gsutil -m cp gs://covid_terra/${run}/${run}_terra_data_table.tsv .


	while read p; do
	  echo $p
	  accession=$(echo $p | awk '{print $1}')
	  echo $accession
	  barcode=$(echo $p | awk '{print $2}')
	  path=$(echo $p | awk '{print $3}')
	  echo $barcode
	  gsutil cat ${path}/*.fastq > CO-CDPHE-${accession}.fastq
	  gzip CO-CDPHE-${accession}.fastq
	done < <(tail -n +3 ${run}_terra_data_table.tsv)

	#cd ..

done
