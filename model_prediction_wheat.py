import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
from keras.models import load_model
import tensorflow as tf
from sklearn.externals import joblib

def run_model(cutoff , out_path):
    test_features = {}
    test_features["aaindex"] = pd.read_table( out_path + "\\sel_AAINDEX.txt", sep="\t")
    test_features["binary"] = np.array( pd.read_csv( out_path + "\\binary.txt", sep="\t").iloc[:,1:])
    test_features["pssm"] = np.array( pd.read_csv( out_path + "\\pssm.txt",sep=" ").iloc[:,1:])
    test_features["aac"] = pd.read_table(out_path + "\\AAC.txt").iloc[:,1:]
    test_features["paac"] = pd.read_table(out_path + "\\PAAC.txt").iloc[:,1:]

    models = {}
    models["aac_rf"] = joblib.load('models\\wheat\\aac_rf_.pkl')
    models["paac_rf"] = joblib.load('models\\wheat\\paac_rf_.pkl')
    models["aaindex_cnn"] = load_model('models\\wheat\\AAINDEX_CNN.h5')
    models["pssm_cnn"] = load_model('models\\wheat\\pssm_CNN.h5')
    models["binary_cnn"] = load_model('models\\wheat\\binary_CNN.h5')
    svm_clf = joblib.load('models\\wheat\\ensemble_svm.pkl')

	#Data preprocessing

    test_results = pd.DataFrame()

    aac_p_scores = models["aac_rf"].predict_proba(test_features["aac"])[:, 1]  #Predicted probabilities of positive case
    test_results["aac"] = aac_p_scores

    paac_p_scores = models["paac_rf"].predict_proba(test_features["paac"])[:, 1]  #Predicted probabilities of positive case
    test_results["paac"] = paac_p_scores

    aaindex_test = np.array(test_features["aaindex"])
    aaindex_test = aaindex_test.reshape(aaindex_test.shape[0], 26, 30)
    aaindex_p_scores = models["aaindex_cnn"].predict_proba(aaindex_test)[:, 1]  #Predicted probabilities of positive case
    test_results["aaindex"] = aaindex_p_scores

    binary_test = np.array(test_features["binary"])
    binary_test = binary_test.reshape(binary_test.shape[0], 32, 21)
    binary_p_scores = models["binary_cnn"].predict_proba(binary_test)[:, 1]  #Predicted probabilities of positive case
    test_results["binary"] = binary_p_scores

    pssm_test = np.array(test_features["pssm"])
    pssm_test = pssm_test.reshape(pssm_test.shape[0], 31, 20)
    pssm_p_scores = models["pssm_cnn"].predict_proba(pssm_test)[:, 1]  #Predicted probabilities of positive case
    test_results["pssm"] = pssm_p_scores

    score = svm_clf.predict_proba(test_results)[:, 1]  #Predicted probabilities of positive case

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



