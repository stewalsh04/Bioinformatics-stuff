# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio import SeqIO
import os
import pandas as pd
import numpy as np
import re
import argparse
import warnings

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')

# !pip install openpyxl

def rename_references(config_df):
    '''
    function to rename references to be unique, does so inplace to dataframe
    args:
    onfig_df (pandas dataframe): dataframe of config.xlsx
    '''
    unique_refs = {i:0 for i in set(config_df["name"])}

    edited_references = []

    for ref in config_df["name"]:
        if unique_refs[ref] == 0:
            print(ref)
            edited_references.append(ref)
            # +2 for logical naming
            unique_refs[ref] += 2
        else:
            edit_ref = "{}_{}".format(ref, unique_refs[ref])
            edited_references.append(edit_ref)
            unique_refs[ref] += 1

    return edited_references


def find_file(file_list, sample_name):
    """Wrapper for a regex to find if a given sample name is present in a list of fastq files
    Args:
        file_list (list): list of fastq files
        sample_name (str): name of sample
    Returns:
        present (bool): true if sample present, else false"""
    
    present = False
    file_pattern = f'{sample_name}_L001_R(1|2)_001.fastq.gz'
    
    # checks there are two unique files following the accepted naming pattern
    matches = []
    for file in file_list:
        match = re.match(file_pattern, file)  #, re.IGNORECASE):
        if match:
            matches.append(match.group())

    if len(set(matches)) == 2:
        present = True
    
    return present
        

def check_sequence_validity(seq):
    """quick function to check sequence is valid
    Args:
        seq (str): DNA sequence
    """
    for i in seq.upper():
        if i not in ["A", "C", "G", "T"]:
            raise ValueError(f"Invalid character found in sequence: {i}")
        

def check_for_barcodes(
    syndna='AAATACAACTGGCCGTCGTTTTACGAAGACCT', 
    fixed_barcode="ACCT", 
    variable_barcodes=[
        "AAAT", 
        "AAAA", 
        "AAAC",
        "AAAG",
        "ATAA",
        "ATAC",
        "ATAG",
        "CAAT",
        "CAAA",
        "CAAC",
        "CAAG",
        "TTCC",
        "TTCG",
        "TTGA",
        "TTGC",
        "TTGG",
        "AGAT",
        "AGAA",
        "AGAC",
        "AGAG",
        "TAGT",
        "TAGC",
        "TAGG",
        "TCGT",
        "TCGA",
        "TGGT",
        "TGGA",
        "TATT",
        "TATG",
        "TATC"
    ]
):
    """This is a quick function to check if a syndna sequence has barcodes at start and end of the sequence
    Args:
        syndna (str): syndna sequence
        fixed_barcode (str): fixed barcode the reference sequence may end in
        variable_barcodes (list): list of barcodes the reference sequence may start with
    Returns:
        start_present (bool): True if variable barcode present at start of synDNA
        end_present (bool): True if fixed barcode present at end of synDNA
    """
    # make sure reference is uppercase
    reference = syndna.upper()

    # find if starting barcode is present
    start_present = False     
    for i in variable_barcodes:
        if reference.startswith(i):
            start_present = True
            break
            
    end_present = reference.endswith(fixed_barcode)
            
    return start_present, end_present
 

def add_vector(
    syndna, 
    start_vector="CAGGCTCTGCTCTTC", 
    end_vector="TTGTGACTCAGGATGCTGT",
    fixed_barcode="ACCT", 
    variable_barcodes=[
        "AAAT", 
        "AAAA", 
        "AAAC",
        "AAAG",
        "ATAA",
        "ATAC",
        "ATAG",
        "CAAT",
        "CAAA",
        "CAAC",
        "CAAG",
        "TTCC",
        "TTCG",
        "TTGA",
        "TTGC",
        "TTGG",
        "AGAT",
        "AGAA",
        "AGAC",
        "AGAG",
        "TAGT",
        "TAGC",
        "TAGG",
        "TCGT",
        "TCGA",
        "TGGT",
        "TGGA",
        "TATT",
        "TATG",
        "TATC"
    ]
):

    """Function to 
insert synthetic DNA into a vector. Will only add if start and end barcodes are 
    both present in syna sequence
    Args: 
        syndna (str): syndna sequence
        start_vector (str): 5' end of vector
        end_vector (str): 3' end of vector
        fixed_barcode (str): fixed barcode the reference sequence may end in
        variable_barcodes (list): list of barcodes the reference sequence may start with
    Returns:
        (str): syndna inserted into vector sequence, in uppercase
    """
    
    start_present, end_present = check_for_barcodes(
        syndna=syndna, fixed_barcode=fixed_barcode, variable_barcodes=variable_barcodes
    )
    
    # check barcodes are present
    if start_present + end_present < 2:
        raise ValueError(
            f'Barcodes missing from synDNA, cannot add vector.\nStart present: {start_present}\nEnd present: {end_present}'
        )
        
    # insert into vector
    return(start_vector.upper()+syndna.upper()+end_vector.upper())


def generate_reference_fasta(
    sequence, name="ref", outdir="./"
):
    """Function to save a sequence into a reference fasta file
    Args:
        sequence(str): sequence to save into fasta file
        name (str): name of sequence
        outdir(str): name of out directory"""

    # ref_obj = SeqRecord(
    #     Seq(sequence),
    #     name=name,
    #     id=name,
    #     description=""
    # )
    # SeqIO.write(ref_obj, os.path.join(outdir, f'{name}.fasta'), "fasta")

    with open(os.path.join(outdir, f'{name}.fasta'), 'w') as f:
        f.write(f'>{name}\n')
        f.write(sequence)


def save_references(config_df, outdir):
    """Function to save references from configuration dataframe and output any failed references
    Args: 
        config_df (pandas dataframe): dataframe containing references and names
        outdir (str): save directory
    return:
     config_clean)df (pandas dataframe): dataframe with invalid references removed"""

    for i in config_df[["name", "synDNA"]][config_df["reference"].isna()].values:
        print(f'Invalid reference, missing one or more barcodes and/or invalid characters: {i[0]}: {i[1]}')

    config_clean_df = config_df.dropna(subset=["reference"]).copy()

    for pair in config_clean_df[["name", "reference"]].values:
        generate_reference_fasta(sequence=pair[1], name=pair[0], outdir=outdir)
    
    return config_clean_df

def generate_references_from_config(
    config_path="./config.xlsx", 
    start_vector="CAGGCTCTGCTCTTC", 
    end_vector="TTGTGACTCAGGATGCTGT",
    fixed_barcode="ACCT", 
    variable_barcodes=[
        "AAAT", 
        "AAAA", 
        "AAAC",
        "AAAG",
        "ATAA",
        "ATAC",
        "ATAG",
        "CAAT",
        "CAAA",
        "CAAC",
        "CAAG",
        "TTCC",
        "TTCG",
        "TTGA",
        "TTGC",
        "TTGG",
        "AGAT",
        "AGAA",
        "AGAC",
        "AGAG",
        "TAGT",
        "TAGC",
        "TAGG",
        "TCGT",
        "TCGA",
        "TGGT",
        "TGGA",
        "TATT",
        "TATG",
        "TATC"
    ],
    outdir="./",
    fastq_dir="./"
):
    """
    Function which reads a config file containing synDNA and names of the synDNA sequences
    Where possibile (i.e. barcodes present) these sequences will then be inserted into the vector
    to form a reference sequence, else np.nan will be used. The config file also contains the sample
    the reference should be paired with and checks for the presence of this sample in the fastq folder
    and eliminates missing samples from further consideration
    Args: 
        config_path (str): path to configuration excel containing synda name, sequence, sample name and 
            wheter or not the syndna should be inserted into a vector
        start_vector (str): 5' end of vector
        end_vector (str): 3' end of vector
        fixed_barcode (str): fixed barcode the reference sequence may end in
        variable_barcodes (list): list of barcodes the reference sequence may start with
        outdir (str): directory to save reference fasta files to
        fastq_dir (str): directory containing fastq files
    Returns:
        config_df (pandas dataframe): dataframe containing name, synDNA sequence and reference sequence
        failure_df (pandas dataframe): dataframe of samples which cannot be run due to absence or failed reference"""
    
    #Extract pipeline mode and seed from config file
    config_df = pd.read_excel(config_path,header=None)
    mode = str(config_df.iloc[0,1])
    seed = int(config_df.iloc[1,1])
    config_df.columns = config_df.iloc[2]
    config_df = config_df.iloc[3:]
    file = open(os.path.join(os.path.dirname(args.config), "mode.txt"), 'a')
    file.write(f'{mode}\n')
    file.write(f'{seed}')
    file.close

    config_df["name"] = rename_references(config_df)
    
    references = []
    for seq, seq_insert in zip(config_df["synDNA"], config_df["PCR"]):
        seq = seq.rstrip()
        if seq_insert != True:
            try:
                check_sequence_validity(seq)     
                references.append(
                    add_vector(
                        syndna=seq,
                        start_vector=start_vector,
                        end_vector=end_vector,
                        variable_barcodes=variable_barcodes,
                        fixed_barcode=fixed_barcode
                    )
                )
            except: 
                references.append(np.nan)
        else:
            try:
                check_sequence_validity(seq)
                references.append(seq)
            except:
                references.append(np.nan)

    config_df["reference"] = references
    
    # find is sample siles are present
    file_found = []
    file_list = [i for i in os.listdir(fastq_dir) if "fastq.gz" in i.lower()]
    
    # put sample names to uppercase because that is how fastqs are
    #config_df['sample'] = config_df["sample"].str.upper()
    
    for sample in config_df["sample"]:
        file_found.append(find_file(file_list, sample))
        
    config_df["sample_found"] = file_found
    
    # print out which files have not been found
    failed_samples = config_df["sample"][config_df["sample_found"] == False].values
    if len(failed_samples) != 0:
        print(
            'The following samples are not present, please check fastq folder:\n{}'.format(
            failed_samples
            )
        )
    
    config_clean_df = config_df[config_df["sample_found"] == True].copy()
    config_clean_df = save_references(config_df=config_clean_df, outdir=outdir)

    sample_fail = config_df[~config_df["sample"].isin(config_clean_df["sample"])].copy()
    ref_fail = config_df[~config_df["name"].isin(config_clean_df["name"])].copy()
    failure_df = pd.concat([sample_fail, ref_fail]).drop_duplicates()
        
    return config_clean_df, failure_df

if __name__ == '__main__':

    parser=argparse.ArgumentParser()

    parser.add_argument("--config", help="path to config file")
    parser.add_argument("--repo", help="path to reference reposititory")
    parser.add_argument("--fastqs", help="path to folder containing fastq files")

    args=parser.parse_args()

    config_df, failure_df = generate_references_from_config(config_path=args.config, outdir=args.repo, fastq_dir=args.fastqs)

    # save edited config_df as a ref.txt which we already have the code to parse
    config_df[["sample", "name", "synDNA"]].to_csv(os.path.join(os.path.dirname(args.config), "ref.txt"), sep='\t', header=False, index=False)
    # save failed_samples dataframe
    failure_df.to_csv(os.path.join(os.path.dirname(args.config), "failures.csv"), header=True, index=False)
