
output_file = open("file_name.tsv","w")
output_file.write('Gene\tA\tA%\tC\tC%\tG\tG%\tT\tT%\tLength\tAT%\tCG%\n')
from Bio import SeqIO

for cur_record in SeqIO.parse("File_name.fasta","fasta"):
    gene_name=cur_record.name
    sequence = cur_record.seq
    A_count = cur_record.seq.count("A")
    C_count = cur_record.seq.count('C')
    G_count = cur_record.seq.count('G')
    T_count = cur_record.seq.count('T')
    length = len(cur_record.seq)
    A_percentage = A_count / length *100
    C_percentage = C_count / length *100
    G_percentage = G_count / length *100
    T_percentage = T_count / length *100
    at_percentage = float(A_count + T_count) / length *100
    cg_percentage = float(C_count + G_count) / length *100
    output_line = '%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\n' % \
                  (gene_name, A_count, A_percentage, C_count, C_percentage, G_count,G_percentage, T_count,T_percentage, length, at_percentage, cg_percentage)
    output_file.write(output_line)



