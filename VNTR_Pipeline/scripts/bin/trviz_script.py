import sys
import os
import glob
from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta


pattern_file = sys.argv[1]
motifs_file=sys.argv[2]


dir_path = os.path.dirname(pattern_file) 
os.chdir(dir_path)

with open(motifs_file, "r") as f:
    motif_list = f.read().splitlines() 

tr_visualizer = TandemRepeatVizWorker()

pattern= f"{pattern_file}"
matching_files = glob.glob(pattern)
for fasta_file in matching_files:
    sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(fasta_file)
#sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(f"{path}*best_hit_combined.fasta")
tr_id = "best_trviz"

tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motif_list)
