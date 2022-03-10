while read name; 
	do
		cutadapt -e 0.1 -O 30 --discard-trimmed -b GTGACGGAAAAACGTATCAAAATGTACAGCAGTTCATCGATGAAGGCAACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACT \
			-B TTTGCCTTCATCGATGAACTGCTGTACATTTTGATACGTTTTTCCGTCACGATCGTCGGACTGTAGAACTCTGAACGTGTAGATC -o Trimmed/$name.1.fq.gz -p Trimmed/$name.2.fq.gz $name*R1*.gz $name*R2*.gz 1>Trimmed/$name.log;
		bbmerge.sh rem ecct in1=Trimmed/$name.1.fq.gz in2=Trimmed/$name.2.fq.gz out=merge_fq/$name.fq.gz outu1=merge_fq/$name.1.fq.gz outu2=merge_fq/$name.2.fq.gz 2>merge_fq/$name.log
		zcat merge_fq/$name.fq.gz merge_fq/$name.1.fq.gz |awk 'NR%4==2 {print}'|sed 's/[GT]GTAC[AC]/	/g'|awk '$3!="" {if( $5!="") {print $2"\t"$3}else {print $2}}' FS=\\t OFS=\\t| \
			awk 'length($1)>40 {print}'> RAW_seq/$name.seq
	done < sample.name
