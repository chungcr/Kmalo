import make_feature_mammals
import make_feature_wheat
import model_prediction_mammals
import model_prediction_wheat

species = "mammals"                     #mammals or wheat
input_FASTA = "test_protein.fasta"      #Input FASTA file
out_path = "output"                     #Output directory

if species=="mammals":
    make_feature_mammals.create_features(input_FASTA, out_path)
    model_prediction_mammals.run_model(0.046, out_path)
else:
    make_feature_wheat.create_features(input_FASTA, out_path)
    model_prediction_wheat.run_model(0.001, out_path)