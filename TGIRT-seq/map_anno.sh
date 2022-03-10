i=$1
RAW=$DATA/MDACCSPNGS-HX-1112/
P_ref=$REF/Pseudomonas
TRIM=$NGS/MDACCSPNGS-HX-1112/Trim
if [ ! -d Trim ] ; then mkdir Trim; fi

echo "Processing ${i}"
echo "Start trimming the reads ${1}"
if [ ! -d $i ] ; then mkdir $i;fi
cutadapt -m 15 -O 3 -U 1 --nextseq-trim=20 -j 48 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGA \
	-o $TRIM/$i.1.fq.gz -p $TRIM/$i.2.fq.gz $RAW/$i.1.fq.gz $RAW/$i.2.fq.gz 1>$i/cutadapt.log
cd $i

if [ ! -d AZPAE12409 ] ; then mkdir AZPAE12409;fi
cd AZPAE12409
echo "Start mapping to AZPAE12409 (no ltrB) ${1}"
bowtie2 -p 48 --local -N 1 -D 20 -L 20 -X 1000 --no-mixed --no-discordant -x $P_ref/AZPAE12409 -1 $TRIM/$i.1.fq.gz -2 $TRIM/$i.2.fq.gz 2>bowtie2.log | samtools view -bS - > align.bam
samtools view -H align.bam > filtered.sam
cp filtered.sam unmapped.sam
samtools view align.bam |awk '{CIGAR1=$6;L1=$0;getline;CIGAR2=$6;L2=$0;if ((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)&&($2!=141)) \
	{print L1"\n"L2 >> "filtered.sam"} else {print L1"\n"L2 >> "unmapped.sam"}}' 
samtools view -@ 48 -bS filtered.sam > filtered.bam
samtools view -bS unmapped.sam | bam2fastx -AQPN -o unmapped.fq.gz -
bedtools bamtobed -bedpe -mate1 -i filtered.bam | awk '{if ($2>$5) $2=$5;if ($3<$6) $3=$6;print $1,$2,$3,$7,0,$9}' FS=\\t OFS=\\t | sort -k1,1 -k2,2n -k3,3n | gzip -f  > filtered.bed.gz
bowtie2 -p 48 --local -N 1 -D 20 -L 20 -X 1000 --no-mixed --no-discordant -x $P_ref/pBL1_vector -1 unmapped.1.fq.gz -2 unmapped.2.fq.gz 2>pBL1.log | samtools view -bS - > pBL1.bam
bedtools coverage -s -F 0.1 -counts -a $P_ref/AZPAE12409.bed.gz -b filtered.bed.gz > AZPAE12409.counts
bedtools genomecov -strand + -bg -i filtered.bed.gz -g $P_ref/AZPAE12409.chrom > plus.bg
bedGraphToBigWig plus.bg $P_ref/AZPAE12409.chrom plus.bw
bedtools genomecov -strand - -bg -i filtered.bed.gz -g $P_ref/AZPAE12409.chrom > minus.bg
bedGraphToBigWig minus.bg $P_ref/AZPAE12409.chrom minus.bw	
rm *.sam
cd ..

echo "Start mapping to PAO1 (no ltrB) ${1}"
if [ ! -d PAO1 ] ; then mkdir PAO1;fi
cd PAO1
bowtie2 -p 48 --local -N 1 -D 20 -L 20 -X 1000 --no-mixed --no-discordant -x $P_ref/PAO1 -1 $TRIM/$i.1.fq.gz -2 $TRIM/$i.2.fq.gz 2>bowtie2.log | samtools view -bS - > align.bam
samtools view -H align.bam > filtered.sam
cp filtered.sam unmapped.sam
samtools view align.bam |awk '{CIGAR1=$6;L1=$0;getline;CIGAR2=$6;L2=$0;if ((CIGAR1!~/^[1-9][0-9]S|[1-9][0-9]S$/) && (CIGAR2!~/^[1-9][0-9]S|[1-9][0-9]S$/)&&(chr1==chr2)&&($2!=141)) \
        {print L1"\n"L2 >> "filtered.sam"} else {print L1"\n"L2 >> "unmapped.sam"}}'
samtools view -@ 48 -bS filtered.sam > filtered.bam
samtools view -bS unmapped.sam | bam2fastx -AQPN -o unmapped.fq.gz -
bedtools bamtobed -bedpe -mate1 -i filtered.bam | awk '{if ($2>$5) $2=$5;if ($3<$6) $3=$6;print $1,$2,$3,$7,0,$9}' FS=\\t OFS=\\t | sort -k1,1 -k2,2n -k3,3n | gzip -f  > filtered.bed.gz
bedtools coverage -s -F 0.1 -counts -a $P_ref/PAO1.bed.gz  -b filtered.bed.gz > PAO1.counts
cd ..

echo "Start reAnnotating AZPAE12409 with PAO1 ${1}"
bedtools intersect -s -wao -f 0.1 -a PAO1/filtered.bed.gz -b $P_ref/PAO1.bed.gz | cut -f 4,10,13,14,15 | \
	awk '{
		L1=$0;
		ID1=$1;
		Over1=$5;
		while (getline == 1) {
			L2=$0;
			ID2=$1;
			Over2=$5;
		if (ID1 == ID2) {
			if (Over2 > Over1) {
				L1=L2;
				ID1=ID2;
				Over1=Over2;
				}
	        	} else if (ID1 != ID2) {
			print L1;
			L1=L2;
			ID1=ID2;
			Over1=Over2;
			}
		}
		print L1;
	}' | gzip -f > PAO1.info.gz
bedtools intersect -s -wao -f 0.1 -a AZPAE12409/filtered.bed.gz -b $P_ref/AZPAE12409.bed.gz | cut -f 4,10,14,16,17 | \
	awk '{
		L1=$0;
		ID1=$1;
		Over1=$5;
		while (getline == 1) {
			L2=$0;
			ID2=$1;
			Over2=$5;
		if (ID1 == ID2) {
			if (Over2 > Over1) {
				L1=L2;
				ID1=ID2;
				Over1=Over2;
				} 
	        	} else if (ID1 != ID2) {
			print L1;
			L1=L2;
			ID1=ID2;
			Over1=Over2;
			}
		}
		print L1;
	}' | gzip -f > AZPAE12409.info.gz
Rscript ../combined.r
gzip -f info.comm
zcat info.comm.gz | awk 'NR>1 && !($2=="." && $6==".") {print $1,$3,$6,$7,$8,$2,$4}' FS=\\t OFS=\\t|sed 's/ /_/g'|tr ";" \\t|cut -f 1-6,10|sed 's/_name_//g'| \
	awk '{
		if ($3==".") {
		type = $2
		name = $6
		tag = $6":"$7
		} else {
		type = $4
		name = $3
		tag = $5";"$6":"$7
		}
		print $1"\t"name"\t"type"\t"tag
	}' | gzip -f > simplified.comm.gz 
zcat simplified.comm.gz | cut -f 2,3|cut -f 1 -d ":"|sort|uniq -c|awk '{print $4"\t"$2"\t"$3"\t"$1}' > reAnno.counts
