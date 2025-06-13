import variant_processing
import matches
import fragment_length
import sys

"""
This script will call all the python scripts in one go for the pipeline, the input are as follows;

filepath: Location of input csv file
filepath_matches: Location of ouput from bioalcidaejdk jvarkit
filepath_fragmentlength: Output taken from fastp
sam_data: Location of samfile
synDNA: Synthetic sequence we synthesized with 4base barcodes
save_path: Location where all output should go
run_id: ID for run this is used for all output file names
ref_seq: The reference sequence itself (rather than ID)
thresh: The optional threshold for the variant discovery script

"""

filepath_vp=sys.argv[1]
filepath_matches= sys.argv[2]
filepath_fragmentlength = sys.argv[3]
sam_data=sys.argv[4]
synDNA=sys.argv[5]
save_path=sys.argv[6]
run_id= sys.argv[7]
ref_seq= sys.argv[8]
thresh= sys.argv[9]

"""

Use below to debug

filepath_vp="C:\\experiment_data\\EVX-45mer-160922\\variants\\EVX-IC45mer\\EVX-IC45mer.var.45mer.csv"
filepath_matches= "C:\\experiment_data\\EVX-45mer-160922\\matches\\EVX-IC45mer.matches.csv"
filepath_fragmentlength = "C:\\experiment_data\\EVX-45mer-160922\\fastp_output\\EVX-45mer-160922.insert_size.csv"
sam_data="C:\\experiment_data\\EVX-45mer-160922\\EVX-IC45mer\\EVX-IC45mer.sam"
synDNA="AATACATGACACGCGCCTGAATTCACAAAGCAAAACTACACCT"
save_path="C:\\experiment_data\\"
run_id= "EVX_29"

"""

synDNA = synDNA.replace('\r', '')
synDNA = synDNA.replace('\n', '')
variant_processing.read_sequence_data(
    filepath=filepath_vp,
    sam_data=sam_data,
    synDNA=synDNA,
    save_path=save_path,
    run_id= run_id,
    ref=ref_seq,
    thresh=thresh
)

matches.read_bin_and_plot_matches(
    filepath=filepath_matches,
    run_id=run_id,
    save_path=save_path,
    remove_first=2
)

fragment_length.plot_and_save_fragment_lengths(
    filepath=filepath_fragmentlength,
    save_path=save_path,
    run_id=run_id,
    read_length=150
)
