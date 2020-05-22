import pandas as pd
import numpy as np
from Bio import SeqIO
import PSSM_mammals
import subprocess

def protein_fasta2peptides(protein_fasta, seq_len, outfile_name):
    #Read protein
    seq_name = []
    seq = []
    fasta_sequences = SeqIO.parse(open(protein_fasta),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seq_name.append(name)
        seq.append(sequence)

    peptide_data = pd.DataFrame({"name":seq_name,"sequence":seq})

    #Split sequence when the protein is "K"

    with open(outfile_name,"w") as w:
        for index, row in peptide_data.iterrows():
            #print(row["sequence"])
            temp_name = row["name"]
            half = int(seq_len/2)
            for pos in range(len(row["sequence"])): 
                if pos<half and row["sequence"][pos]=="K":
                    name_n_position = temp_name + "_" +str(pos+1)
                    w.write(">" + name_n_position)
                    w.write("\n")

                    temp = row["sequence"][:pos]+row["sequence"][pos:pos+half+1]
                    temp = temp.rjust(seq_len,'X')

                    w.write(temp)
                    w.write("\n")
                    continue

                if row["sequence"][pos]=="K":
                    name_n_position = temp_name + "_" +str(pos+1)
                    w.write(">" + name_n_position)
                    w.write("\n")
                    temp = row["sequence"][pos-half:pos+half+1]
                    temp = temp.ljust(seq_len,'X')
                    w.write(temp)
                    w.write("\n")
					
def del_mid_K(in_file, seq_len, out_file):
    f = open(in_file)
    with open(out_file,'w') as w:
        i = 0
        for line in f.readlines():
            if i % 2:
                w.write(line[:int(seq_len/2)] + line[int(seq_len/2)+1:])
            else:
                w.write(line)
            i = i+1
    print("%s is done !" %out_file)
	

def create_features(input_FASTA, out_path):
	# make input file to sequences
	# the best length of feautres:
	# AAINDEX: 33
	protein_fasta2peptides(input_FASTA, 33, out_path+"\\window33.fasta")
	test_file = out_path+"\\window33_no_midK.fasta"
	del_mid_K(out_path+"\\window33.fasta", 33, test_file)
	output_feature = out_path + "\\raw_AAINDEX.txt"
	command = "python iFeature-master\\iFeature.py --file " + test_file + " --type AAINDEX --out " + output_feature
	print(command)
	result = subprocess.check_output(command, shell=True) 
	#retrieve selected
	with open("data\\final_PCP_mammal.txt") as f:
		pcp = f.readline()
	raw_aaindex = pd.read_table(output_feature)
	my_feature = raw_aaindex[raw_aaindex.filter(regex=pcp).columns]
	output_feature = out_path + "\\sel_AAINDEX.txt"
	my_feature.to_csv(output_feature,sep="\t", index=False)	
	
	# BINARY 31
	protein_fasta2peptides(input_FASTA, 31, out_path+"\\window31.fasta")
	del_mid_K(out_path+"\\window31.fasta", 31, out_path+"\\window31_no_midK.fasta")
	test_file = out_path+"\\window31_no_midK.fasta"
	output_feature = out_path + "\\binary.txt"
	command = "python iFeature-master\\iFeature.py --file " + test_file + " --type BINARY --out " + output_feature
	print(command)
	result = subprocess.check_output(command, shell=True)

	# PSSM : 33
	print("Making PSSM......\n")
	protein_fasta2peptides(input_FASTA, 33, out_path+"\\window33.fasta")
	PSSM_mammals.make_PSSM(out_path + "\\window33.fasta", out_path)


	print("ALL_features_DONE!!")

