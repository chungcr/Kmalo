import os
from itertools import islice
from Bio import SeqIO
import subprocess

def make_PSSM(input_FASTA, out_path):
	'''
	 The following codes are used to read FASTA files and generate several PSSM files to the folder
	 
	 fasta_path：the directory of the fasta file (A fasta file can contain several protein sequences)
	 
	 Database_path： the directory and name of the database
	 
	 pssm_path： the directory which used to save the generated PSSM files
	 
	 name_of_fasta： Name of FASTA file
	 
	 feature_path : the directory of the final gererated features
	 
     Program will open the FASTA file to read the sequences accordingly and use psi-blast to generate PSSM
    Finally, the generated feature will be saved in feature_path
	'''
	
	os.makedirs(out_path + "\\local_PSSM")
	
	fasta_path = out_path
	Database_path = os.path.join("PSI_blast","db","PLMD.fasta")
	pssm_path = os.path.join(out_path,"local_PSSM")
	feature_path = out_path
	name_of_fasta = "window33.fasta"
	
	seq_len = 33
	
	with open(feature_path + "\\pssm.txt", "w") as w:
		amino_acids = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
		header = []
		header.append("#")
		for i in range(1,seq_len+1):
			for AA in amino_acids:
				header.append("Pos." + str(i) + "." + AA)
				#print("Pos." + str(i) + "." + AA)
		w.write(" ".join(header) + "\n")

		fasta_sequences = SeqIO.parse(open(fasta_path + "\\" + name_of_fasta),'fasta')
		for fasta in fasta_sequences:
			name, sequence = fasta.id, str(fasta.seq)
			file_name = pssm_path + "\\" + "{}.pssm".format(name)
			print("dealing with %s" %name)
			with open(fasta_path + "\\temp.fa","w") as out_file:
				out_file.write(">" + name + "\n" + sequence)
			psi_command = "PSI_blast\\bin\\psiblast " + " -query " + fasta_path + "\\temp.fa" + " -db " + Database_path  + " -num_iterations " + "3"   + " -out_ascii_pssm " + pssm_path + "\\{}.pssm".format(name)
			result = subprocess.check_output(psi_command, shell=True)
			w.write(name + " ")

			if os.path.isfile(file_name): #if exists PSSM
				with open(file_name, "r") as f:
					new_line = []
					for line in islice(f,3,3+seq_len):
						new_line = new_line + list(map(int,line.strip().split()[2:22]))
					pssm = " ".join(map(str,new_line))
					w.write(pssm)
					w.write("\n")

			else: # otherwise replace them by 0
				w.write(" ".join(map(str,[0]*(seq_len*20))))
				w.write("\n")