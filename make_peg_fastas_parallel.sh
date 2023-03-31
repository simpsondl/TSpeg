#### Script arguments
# $1 indir, not full path
# $2 original sample fastq
# $3 demultiplexed fq directory
# $4 split peg id file

ndx=${4##tmp_}

while read p; 
do 
	if [[ $p == '*' ]]; 
	then 
		continue; 
	else 
		samp=${1%%/}
		grep $p ${samp}_tokeep.txt | cut -f1 > keeptmp_${ndx}; 
		~/tools/seqtk/seqtk subseq \
		${3}/$2 \
		keeptmp_${ndx} > peg_fastas/untrimmed/${p}_${samp}.fq;
		~/tools/seqtk/seqtk trimfq -b 6 -e 14 peg_fastas/untrimmed/${p}_${samp}.fq \
		> peg_fastas/trimmed/${p}_${samp}_trimmed.fq
	fi; 
done < $4
