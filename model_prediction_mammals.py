import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from keras.models import load_model
import tensorflow as tf

def run_model(cutoff , out_path):
	test_features= {}
	test_features["binary"] = np.array( pd.read_csv( out_path + "\\binary.txt", sep="\t").iloc[:,1:])
	test_features["aaindex"] = np.array( pd.read_csv( out_path + "\\sel_AAINDEX.txt", sep="\t"))
	test_features["pssm"] = np.array( pd.read_csv( out_path + "\\pssm.txt",sep=" ").iloc[:,1:])

	model = load_model("models\\mammals\\model.h5")

	#Data preprocessing

	#binary
	binary_train_raw = test_features["binary"]
	binary_test = binary_train_raw.reshape(binary_train_raw.shape[0], 30, 21)

	# aaindex
	aaindex_train_raw = test_features["aaindex"]
	aaindex_test =  aaindex_train_raw.reshape(aaindex_train_raw.shape[0], 32, 46)

	# pssm
	pssm_train_raw = test_features["pssm"]
	pssm_test = pssm_train_raw.reshape(pssm_train_raw.shape[0], 33, 20)


	score = model.predict([binary_test, aaindex_test, pssm_test])[:,1] #Predicted probabilities of positive case

	# save predict result in csv file
	seq_name = []
	seq = []
	fasta_sequences = SeqIO.parse(open(out_path + "\\window33.fasta"),'fasta')
	for fasta in fasta_sequences:
		name, sequence = fasta.id, str(fasta.seq)
		seq_name.append(name)
		seq.append(sequence)

	data = pd.DataFrame({"name":seq_name,"sequence":seq})
	final_df = data["name"].str.split("_" ,expand = True)
	final_df.columns = ["Protein", "Position"]
	final_df["Sequence"] = data["sequence"]
	final_df["Score"] = score
	final_df = final_df[final_df["Score"]>cutoff]
	final_df["Score"] = np.round(final_df["Score"],3)
	final_df.to_csv(out_path + "\\" + out_path.split("\\")[-1] + ".csv", index=False)


