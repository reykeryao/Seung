cutadapt -m 40 -O 3 --nextseq-trim=20 -j 8  --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCACX  -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAX -o Trim/$1.1.fq.gz -p Trim/$1.2.fq.gz $1.1.fq.gz $1.2.fq.gz 1>log/$1.log
bash bbmerge-auto.sh in1=Trim/$1.1.fq.gz in2=Trim/$1.2.fq.gz out=Merged/$1.merge.fq.gz outu1=Unmerged/$1.1.fq.gz outu2=Unmerged/$1.2.fq.gz 2>log/$1.merge.log
cutadapt -j 8 -O 30 -g XGCAATAATCTATACAATACAACACATACAAACAAATTCTTAAGGTCCCAA --info-file Info/$1.info --action none Merged/$1.merge.fq.gz -o temp 1>log/$1.template.log
awk '($2+1)>0 {print $7}' FS=\\t OFS=\\t Info/$1.info |sort|uniq -c|sort -k2,2|awk '{print $2"\t"$1}' > Summary/$1.summary
rm temp
awk '($2+1)>0 {print length($7)+50}' FS=\\t OFS=\\t Info/$1.info |sort|uniq -c|sort -k2,2n|awk '{print $2"\t"$1}' > Summary/$1.len
