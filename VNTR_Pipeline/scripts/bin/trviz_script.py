import sys
import os
from trviz.main import TandemRepeatVizWorker
from trviz.utils import get_sample_and_sequence_from_fasta


path=sys.argv[1]
motifs_file=sys.argv[2]

os.chdir(path)

with open(motifs_file, "r") as f:
    motif_list = f.read().splitlines() 

tr_visualizer = TandemRepeatVizWorker()
sample_ids, tr_sequences = get_sample_and_sequence_from_fasta(f"{path}/best_hit_combined.fasta")
tr_id = "best_trviz"

tr_visualizer.generate_trplot(tr_id, sample_ids, tr_sequences, motif_list)
