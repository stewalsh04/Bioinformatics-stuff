from pymsaviz import MsaViz
import sys

forward_aligned = sys.argv[1]
sample_ID = sys.argv[2]
path = sys.argv[3]

def get_alignment_figure (path, alignment_file, sample_ID):
    msa_file = (alignment_file)
    figure = MsaViz(msa_file, wrap_length=60, show_count=True,show_consensus=True)
    figure.savefig("{}/{}_multiple-alignment.png".format(path, sample_ID))

get_alignment_figure (path, forward_aligned, sample_ID)
