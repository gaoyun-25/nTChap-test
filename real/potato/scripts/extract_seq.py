import sys
import gzip
# f_in = open('test.txt', 'rb')
# f_out = gzip.open('test.txt.gz', 'wb')
# python extract_seq.py GCA_011022315.1_ASM1102231v1_genomic.fna ChrSc2.fa
def extract_fasta_sequences(input_file, output_file, target_ids):
    if input_file.endswith(".gz"):
        with gzip.open(input_file, 'rb') as infile, open(output_file, 'w') as outfile:
            write_sequence = False
            for line in infile:
                line = line.decode()
                if line.startswith('>'):
                    # 检查ID是否在目标ID列表中
                    sequence_id = line[1:].strip()
                    if sequence_id in target_ids:
                        write_sequence = True
                        outfile.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    outfile.write(line)
    else:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            write_sequence = False
            for line in infile:
                if line.startswith('>'):
                    # 检查ID是否在目标ID列表中
                    sequence_id = line[1:].strip()
                    if sequence_id in target_ids:
                        write_sequence = True
                        outfile.write(line)
                    else:
                        write_sequence = False
                elif write_sequence:
                    outfile.write(line)

if __name__ == "__main__":
    input_fasta_file = sys.argv[1]
    # input_fasta_file = "/home/gaoyun_amd/poly/data/potato/reference/potato/DM_1-3_516_R44_potato_genome_assembly.v6.1.fa.gz"
    output_fasta_file = sys.argv[2]
    # output_fasta_file = "/home/gaoyun_amd/poly/data/potato/reference/potato/ref.fa"
    target_sequence_ids = sys.argv[3:]
    # target_sequence_ids = ["chr02"]
    # target_sequence_ids = ["chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12"]  # 替换为你的目标ID列表

    extract_fasta_sequences(input_fasta_file, output_fasta_file, target_sequence_ids)