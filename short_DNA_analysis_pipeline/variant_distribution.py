"""
Created on 19/04/2022

This module will plot aligned read distribution; the number % of filtered reads distribution 
that were successfully aligned to the ref using BWA-MEM

.. sectionauthor:: Supuli Jayaweera
"""
# import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
import seaborn as sns
from difflib import SequenceMatcher
# from pathlib import Path
import re
import os
import base64
from PIL import Image
from io import BytesIO
import io
import heatmap

def get_number_perfect(ref,counts_table,percent):
    """Function to get the number of perfect sequences
    The function has to take account of the maximum read length of the sequencer, if the read length is exceeded but paired end sequncing is used to double the read length, 
    the perfect sequence will be found accross two variants, once in the forward and once in the reverse (see note below). If double the read length is exceeded the number 
    of perfect matches cannot be calculated.
    Args:
        ref = Reference seuqnece
        counts_table = Pandas dataframe with the counts for each sequence variant
        ('max_read' will be added in the future when we get a NextSeq2000)
    Returns:
        The number of perfect matches found

    Note: When DNA length is over the maximum read length but under double the read length (i,e paired end max) there is a level of uncertainty, so the followed is done:     
    In this case a perfect match is split across the forward and reverse reads, to get the number of perfect matches, the number in the forward and reverse must be determed,
    then the lower of the two is taken as you must have a perfect read in the forward and reverse to get a full read (this is were the uncertainty is as the reads most 
    likely will overlap across the two reads), the lower number is then doubled to take account for the fact we are not couting the reads in the other strand. This is currently
    ONLY set up to calculate perfect macthes for NGS generated from PCR generated libraries, it's assumed longer lengths of DNA will not be ligated, so this feature has not been added
    yet. To update the code you have to take account for the vector and barcode sequences when calculating the maximum read length because part of the vector and barcode is sequenced,
    therefore you have to determine exactly how much vector is sequenced and if this is consistent across runs to update the code.
    """

    print("Getting the number of perfect sequences")
    total_perfect_counts = 0
    counts_table = counts_table
    percent = percent
    max_read = 151 # This can be made a pipeline option when we get a NextSeq2000, i.e it will need to increase to ~300

    #Determine if we want to add up the percentages or counts
    if percent == True:
        increaseby = 2
    else:
        increaseby = 1

    #If the reference exceeds double the maximum read length, i.e forward and reverse, the number of perfect cannot be found
    if len(ref)> max_read*2:
        print("Maximum sequncer length exceeded, cannot determine if perfect match is present")
        return 0

    #Below is code to determine the number of perfect matches when the DNA length is over the maximum read length but under double the read length (i,e paired end max)
    elif len(ref) > max_read <= (max_read*2):
        print("Read legnth in one direction exceeded, number of perfect matches may not be 100% accuracte (But very close)")
        print("NOTE Ligated libraries are not currently supported at this read length when calculating the number of perfect matches")
        forwardmax = ref[0:max_read]
        reversemax = ref[(len(ref)-max_read):len(ref)]
        justforward_perfect=(counts_table[counts_table['synths'].str.contains(forwardmax,True)])
        justreverse_perfect=(counts_table[counts_table['synths'].str.contains(reversemax,True)])
        if len(justreverse_perfect) > 1 or len(justforward_perfect) > 1:
            print("ERROR! Multiple perfect sequences in different rows, check the maximum sequencing length")
            return 0
        elif len(justreverse_perfect) == 0 and len(justforward_perfect) == 0:
            print ("No perfect found")
            return 0
        elif len(justreverse_perfect) == 1 and len(justforward_perfect) == 0:
            print("Perfect only found in one direction")
            return 0
        elif len(justreverse_perfect) == 0 and len(justforward_perfect) == 1:
            print("Perfect only found in one direction")
            return 0
        elif len(justreverse_perfect) == 1 and len(justforward_perfect) == 1:
            if int(justforward_perfect.iloc[0,increaseby]) < int(justreverse_perfect.iloc[0,1]):
                total_perfect_counts += (int(justforward_perfect.iloc[0,increaseby]))*2
            elif int(justreverse_perfect.iloc[0,increaseby]) < int(justforward_perfect.iloc[0,1]):
                total_perfect_counts += (int(justreverse_perfect.iloc[0,increaseby]))*2
            elif int(justforward_perfect.iloc[0,increaseby]) == int(justreverse_perfect.iloc[0,1]):
                total_perfect_counts += (int(justreverse_perfect.iloc[0,increaseby]))*2
    
        return(total_perfect_counts)

    #The ref is less than the maximum read length
    else:
        just_perfect=(counts_table[counts_table['synths'].str.contains(ref,True)])

        if len(just_perfect) > 1:
            print("ERROR! Multiple perfect sequences in different rows, check the maximum sequencing length")
            return 0
        elif len(just_perfect) == 0:
            print ("No perfect found")
            return 0
        elif len(just_perfect) == 1:
            total_perfect_counts += int(just_perfect.iloc[0,increaseby])

        return(total_perfect_counts)

def deletion_pie_chart(seqframe, max_deletions=2, count_col="count"):
    """Function to plot a pair of pie
    pie charts showing number of reads with deletions and a breakdown of those deletions
    by number of deletions present in each read
    Args:
        seqframe (pandas dataframe): dataframe containing sequence variant information
        max_deletions (int): maximum number of deletions to give their own section of the pie chart
            sequences with more deletions in this will be rolled into a single category of > max_deletions
        count_col (str): column of seqframe containing readcount data
    Returns: 
        fig (matplotlib figure): figure showing pair of pie charts
        deletion_frame (pandas dataframe): dataframe associated with piechart
    """
    #copy dataframe
    seqframe = seqframe.copy()
    # calculate number of deletions in each variant
    seqframe["dels"] = [
        i.count("_") if i.count("_") <= max_deletions else ">{}".format(max_deletions) for i in seqframe["synths"]
    ]
    # create a summary deletions frame with readcount for each number of deletions
    read_count = []
    del_nos = list(seqframe["dels"].unique())
    for i in del_nos:
        read_count.append(seqframe[count_col][seqframe["dels"] == i].sum())
    deletion_frame = pd.DataFrame({"read_count": read_count, "del_no": del_nos})
    # generate plots
    fig, axes = plt.subplots(ncols=2, figsize=(10,5))
    fig.suptitle("Pie charts showing breakdowns of deletions. Total read count: {}".format(sum(read_count)))
    axes_flat=axes.flatten()
    # first pie: deletion vs no deletion:
    no_dels = [
        deletion_frame["read_count"][deletion_frame["del_no"]!=0].sum(),
    deletion_frame["read_count"][deletion_frame["del_no"]==0].sum()
    ]
    axes_flat[0].pie(
        x=no_dels, 
        explode=[i*0.05 for i in range(len(no_dels))], 
        labels=["Deletion", "No Deletion"],
        autopct='%1.1f%%'
    )
    axes_flat[0].axis("equal")
    axes_flat[0].set_title("Fraction with/without deletions")
    # second pie: deletion subset
    trimmed_frame = deletion_frame[deletion_frame["del_no"] != 0].copy()
    axes_flat[1].pie(
        x=trimmed_frame["read_count"], 
        explode=[i*0.05 for i in range(len(trimmed_frame))], 
        labels=trimmed_frame["del_no"],
        autopct='%1.1f%%'
    )
    axes_flat[1].axis("equal")
    axes_flat[1].set_title("Breakdown of sequences with a deletion")
    fig.tight_layout(w_pad=5.0)
    readstotal = deletion_frame["read_count"].sum()
    deletion_frame['% Reads'] = round((deletion_frame['read_count']/readstotal)*100,2)

    pic_IObytes = io.BytesIO()
    fig.savefig(pic_IObytes,  format='png')
    pic_IObytes.seek(0)
    pic_hash = base64.b64encode(pic_IObytes.read())
    piechartb64 = pic_hash.decode('utf-8')

    return fig, deletion_frame, piechartb64

def trim_barcodes(
    reference='AAATACAACTGGCCGTCGTTTTACGAAGACCT', 
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
    """This is a quick function to return the trim indices for removing the barcodes from the start and ends
    of a reference sequence. These indices are None
    Args:
        reference (str): reference sequence
        fixed_barcode (str): fixed barcode the reference sequence may end in
        variable_barcodes (list): list of barcodes the reference sequence may start with
    Returns:
        start_trim (None or int): index to trim start from
        end_trim (None or int): index to trim end to
    """
    # make sure reference is uppercase
    reference = reference.upper()
    
    start_seq = reference[:4]
    end_seq = reference[-4:]
    
    # if you use None as an index in a string slice it acts like no character
    start_trim = None
    end_trim = None
    
    # using lengths for trimming index in case we change the length of barcodes in the future
    if reference.endswith(fixed_barcode):
        end_trim = -len(fixed_barcode)
        
    for i in variable_barcodes:
        if reference.startswith(i):
            start_trim = len(i)
            break
            
    return start_trim, end_trim
    
    

def calculate_per_base_error(
    variant_results, 
    reference="ACAACTGGCCGTCGTTTTACGAAG", 
    fixed_barcode="ACCT",
    trim=True,
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
    """Function to calculate fraction of expected synthesised based which are errors
    Uses the syndna sequence supplied in the NGS pipeline rather than the whole synthesised
    sequence (which would include vector and barcodes). Offers option to remove barcodes 
    prior to error calculation
    Args:
        variant_results (pandas dataframe): dataframe output of variant_distribution.py
        reference (str): reference sequence in uppercase
        trim (bool): If True trim off barcodes prior to error calculation
        reference (str): reference sequence
        fixed_barcode (str): fixed barcode the reference sequence may end in
        variable_barcodes (list): list of barcodes the reference sequence may start with
    Returns:
        variant_df (pandas dataframe): variant_results frame with added columns for errors and base counts
        error_fraction (float): fraction of expected number of synthesised bases that are a synthesis error
    
    """
    variant_df = variant_results.copy()
    
    if not reference.isupper():
        raise ValueError("Reference sequence {} entered must be uppercase".format(reference))

    
    # handle barcode trimming
    if trim:
        start_trim, end_trim = trim_barcodes(
            reference=reference,
            fixed_barcode=fixed_barcode,
            variable_barcodes=variable_barcodes
        )
    else:
        start_trim = None
        end_trim = None
        

    # create a trimmed reference sequence, nb is the same as reference if no barcodes found
    trimmed_reference = reference[start_trim:end_trim]
    
    # handle error counting
    insertions = []
    deletions = []
    mismatches = []
    
    for seq in variant_df["synths"]:
        # extract synthesised sequence without insertions and calculate number of insertions
        insert_free = "".join([base for base in seq if not base.islower()])
        
        # handle terminal deletions by padding right ride of string with deltion characters
        if len(insert_free) < len(reference):
            insert_free.ljust(len(reference), "_")
                  
        insertions.append(len(seq) - len(reference))
        
        # trim sequence
        insert_free = insert_free[start_trim:end_trim]
        
        # dealing with other errors
        deletion_count = 0
        mismatch_count = 0
        for synbase, base in zip(insert_free, trimmed_reference):
            if synbase != base:
                if synbase == "_":
                    deletion_count += 1
                else:
                    mismatch_count += 1
                    
        deletions.append(deletion_count)
        mismatches.append(mismatch_count)
        
    # adding the errors of the variant frame
    variant_df["insertions"] = insertions
    variant_df["deletions"] = deletions
    variant_df["mismatches"] = mismatches
    variant_df["error_sum"] = variant_df[["insertions", "deletions", "mismatches"]].sum(axis=1)
    
    # calculating base errors and expected bases
    variant_df["total_base_errors"] = variant_df["error_sum"] * variant_df["count"]
    variant_df["expected_bases"] = variant_df["count"] * len(trimmed_reference)
    
    # calculating error fraction
    error_fraction = variant_df["total_base_errors"].sum()/variant_df["expected_bases"].sum()
    
    return variant_df, error_fraction


def plot_cumsum(
    base_counts,
    len_col="length",
    red_level=1000,
    green_level=10000,
    expected_length=150,
    alpha=0.08,
    title="Culmulative sum of reads below a given length",
    ylabel="Cumulative sum of reads",
    xlabel="Number of bases in reads",
    savepath="./",
    filename="cumulative_graph.png",
    left_len=None,
    right_len=None
    ):
    """
    Function which plots the cumulative sum of number of reads against read length, flagging red, amber
    and green regions (i.e. warning on number of reads), as well as the expected read length. This plot is saved
    Args:
        base_counts (pandas dataframe): dataframe associating number of bases in a given read with the 
            number of reads with this count. 
        len_col (str): column containing read length data
        green_level (int): minimum number of reads to consider an experiment without issue in
            terms of read number
        red_level (int): number of reads below which to consider the results concerning
        expected_length (int): number of bases expected in the sequence
        alpha (float): transparency for background regions
        title (str): title for plot
        ylabel (str): label for y axis
        xlabel (str): label for x axis
        savepath (str): directory to save file in
        file_name (str): name of save file
        left_len (int): if not None mark on graph the length of the left side of read 
            (for assembly of longer sequences)
        right_len (int): if not None mark on graph the length of the left side of read 
            (for assembly of longer sequences)
    """
    # copying and sorting frame by read length
    cumsum_frame = base_counts.copy()
    cumsum_frame.sort_values(by=len_col, inplace=True)
    cumsum_frame["cumsum"] = cumsum_frame["count"].cumsum()
    
    # generating a plot
    fig, ax = plt.subplots(figsize=(12,8))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # adding the culmulative frequency
    sns.lineplot(data=cumsum_frame, x=len_col, y="cumsum", ax=ax, label= "Cumulative number of reads", ci=None)
    
    # generating flagging regions of shaded background
    ax.axhspan(
        min(ax.get_ylim()), 
        red_level, color="red", 
        alpha=alpha,
        label="Total reads: Red Warning"
    )
    ax.axhspan(
        red_level, 
        green_level, 
        color="orange", alpha=alpha,
        label="Total reads: Amber Warning"
        )
    ax.axhspan(
        green_level, 
        max(ax.get_ylim()), 
        color="green", 
        alpha=alpha,
        label="Total reads: Green (no warning)")
    
    # add expected read length    
    ax.axvline(expected_length, linestyle=(0, (5, 5, 1 ,5)), label="Expected length of read", color="red", linewidth=1)
    # add expected read lengths for left and right if used
    if left_len is not None:
        ax.axvline(
            left_len, linestyle=(0, (1, 3)), label="Expected length of left read", color="darkgrey", linewidth=1
        )
    if right_len is not None:  
        ax.axvline(
            right_len, linestyle=(0, (5, 3)), label="Expected length of right read", color="darkgray", linewidth=1
        )

    # add legend
    ax.legend(bbox_to_anchor=(1.05, 1))
    
    plt.tight_layout()
    
    plt.savefig(os.path.join(savepath, filename))


def get_clipped_reads(read, cigar):
    """get the clipped reads- clipped reads do not have the soft clipped bases in the return sequnces
    Args:
        read (str): a sequence
        cigar (str): a string representing the alignment; matches, mismatches, deletions, insertions
            and whether or not to use soft clipping. Soft clipping means the bases are retained but not 
            scored as part of the alignment
    Returns
        clipread (str): a clipped sequence string"""
    vals = re.findall(r"(\d+)S", cigar)
    if len(vals)==2:
        # read = preread[int(vals[0]):]
        # read = read[:-int(vals[1])]
        clipread = read[int(vals[0]):-int(vals[1])]
    elif len(vals)==1:
        if cigar[-1]=='S':
            clipread = read[:-int(vals[0])]
        else:
            clipread = read[int(vals[0]):]
    elif len(vals)==0:
        clipread = read
    else:
        raise ValueError(f"CIGAR string has more than 2 soft clippings.")
    return clipread    


def get_final_reads(clipread, cigar, st_pos, ref_length):
    """This will get the soft clipped read and add the deletions/insertions based on the cigar.
    Finally it will add the start /end deletion bases to match to the length of the reference.
    Args:
        clipread (str): a clipped sequence string
        cigar (str): a string representing the alignment; matches, mismatches, deletions, insertions
                and whether or not to use soft clipping. Soft clipping means the bases are retained but not 
                scored as part of the alignment
        st_pos (int): start position for sequence
        ref_length (int): length of reference sequence
    Returns:
        clipread (str): clipped sequence string with insertions and deletions added
    """
    vals = re.findall(r"(\d+)(\D)", cigar)
    # eg vals = [('100', 'S'), ('9', 'M'), ('2', 'D'), ('1', 'I'), ('10', 'M'), ('12', 'S')]
    count = 0
    for item in vals:
        if item[1] != 'S':
            if item[1] =='M':
                count += int(item[0])
            elif item[1] =='D':
                clipread = clipread[:count]+'_'*int(item[0])+clipread[count:]
                count += int(item[0])
            elif item[1] =='I':
                insert_temp = clipread[count:count+int(item[0])].lower()
                clipread = clipread[:count]+insert_temp+clipread[count+int(item[0]):]
                count += int(item[0])
            else:
                print(f'CIGAR: not implemented chacters-{item[1]}')
                count += int(item[0])

    #add start/end shifts to the string
    if st_pos != 1:
        clipread = '_'*(st_pos-1)+clipread
    if len(re.findall(r"[A-Z_]", clipread))<ref_length:
        add_ = ref_length-len(re.findall(r"[A-Z_]", clipread))           
        clipread = clipread + '_'*add_
    return clipread

def get_synthesised_reads(clipread='________________TTCTGCTAGTCTAGACGTACACCTTTGTGACTCAGGATGC', 
                          ref_seq='TTCTGCTAGTCTAGACGTACACCT',
                          syn_seq='TTCTGCTAGTCTAGACGTACACCT'):
    """
    extract the synthesised region from the clipped read.
    Args:
        clipread (str): clipped sequence string with insertions and deletions added
        ref_seq (str): reference sequence
        syn_seq (str): sythetic sequence to extract from the reference sequence
    Returns:
        synth_reg (str): sequence extracted from the clipped sequence (read)
    """
    match = SequenceMatcher(None, ref_seq, syn_seq).find_longest_match(0, len(ref_seq), 0, len(syn_seq))
    if match[0] + match[2] > len(ref_seq):
        raise ValueError('Synthesis Region is not fully covered in the reference')
    syn_st=match[0]
    
    # if inserts are there then adjustments have to be made to the starting and end regions of the synthesised read
    before_insert = re.findall(r"[a-z]", clipread[0:syn_st])
    during_insert = re.findall(r"[a-z]", clipread[syn_st+len(before_insert):syn_st+len(syn_seq)+len(before_insert)])
    synth_reg = clipread[syn_st+len(before_insert):syn_st+len(syn_seq)+len(before_insert)+len(during_insert)]    
    return synth_reg      

def get_synthesised_size(read):
    """Return the bases (not counting the deletions)
    Args: 
        read (str): a sequence, can have deletions
    Returns:
        (int): length of sequence not including deletions
    """
    return len(re.findall(r"[ A-Za-z]", read)) 


def get_aligner_score(clipread='ACAGGCTCTGCTCTTCTTCTGCTAAATG_', ref_seq='ACAGGCTCTGCTCTTCTTCTGCTAGTCTAGACGTACACCTTTGTGACTCAGGATGC'):
    """This will generate a score sequence that can be mapped to the synthesis ref and the reads bases 
    that falls in this region. The array will be as the same length as the synthesis region.
    Score 1 = perfect match base
    Score -1 = Deletion
    Score 0.5 = Mismatch
    Penalty of -0.3 is added for all adjecnet bases that has an insert
    Args:
        clipread (str): clipped sequence with insertions and deletions added
        ref_seq (str): reference sequence
    Returns:
        score_array (numpy array): array of scores for the clipped sequence
    """
    score_array = np.zeros(len(ref_seq))
    count = 0
    
    for i, char in enumerate(list(map(str,ref_seq))):
        if clipread[count].islower():
            score_array[i-1] -= 0.3
            score_array[i] -= 0.3
            count +=1
            
        if clipread[count]==char:
            score_array[i] += 1
            count +=1
        else:
            if clipread[count]=='_':
                score_array[i] -= 1
                count +=1
            elif clipread[count] != char:
                score_array[i] += 0.5
                count +=1
    if np.where(score_array==0)[0].size !=0:  #this is just to assert the read scoring array is correct
        print('Error in aligner score array')    
    return score_array


def get_aligner_df(result_df, ref_seq, col_key='clipped'):    
    """
    create a dataframe with a score array for each read type. The dataframe should be the output of get_result_df

    Args:
        result_df (pandas dataframe): output of get_result_df, which is processed reads from a SAM output file along with
        the metrics calculated from the functions in this script
        ref_seq (str): reference sequence
        col_key (str): column in result_df containing clipped sequences
    Returns:
        score_df (pandas dataframe): A dataframe with a score array for each read type
    
    """
    # get the align score array
    val = result_df.apply(lambda x: get_aligner_score(x[col_key], ref_seq=ref_seq),axis=1)
    score_df = pd.DataFrame(val)
    
    # create the score df
    score_df = pd.DataFrame(score_df[0].to_list(), columns = list(map(str,ref_seq)), index=score_df.index)    
    return score_df


def unampped_aligner(read='TTCTGCACCTTTGTGACCCGAGCCCACGT', ref_seq='ACAGGCTCTGCTCTTCTTCTGCTAGTCTAGACGT'):
    """get the unammped reads aligned to reference
    Args:
        read (str): sequence
        ref_seq (str): reference sequence
    Return:
        read (str): read sequence aligned to reference sequence
    """
    match = SequenceMatcher(None, ref_seq, read).find_longest_match(0, len(ref_seq), 0, len(read))
    read = read[match[1]:match[1]+match[2]]
    read = '_'*int(match[0])+read
    read = read+ '_'*int(len(ref_seq)-len(read))
    return read


def get_result_df(filepath='output_short.txt', 
                  syn_seq='TTCTGCTAGTCTAGACGTACA', 
                  ref_seq='TTCTGCTAGTCTAGACGTACACCT', 
                  mapped = True,
                  pp_only=True):
    """
    Process all reads in a SAM output file and give the reads aditribution results for plotting
    for mapped version = all mapped sequence of quality score more than 30 is considederd. 
                         both reads(R1 and R2) shoulded pass the quality to be selected.
    Args:
        filepath (str): path to SAM file
        syn_seq (str): synthetic sequence
        ref_seq (str): reference sequence
        mapped (bool): if True get clipped reads, else sequence is aligned first
        pp_only (bool): if true the reads are properly paired, i.e. the alignment found its forwards and reverse reads
    Returns:
        read_dist_df (pandas dataframe): dataframe containing reads distribution results for plotting

    
    """
    #read the aligner output file and rename columns
    read_dist_df = pd.read_csv(filepath, usecols=range(10), sep='\t', header=None,skiprows=2)
    read_dist_df = read_dist_df[[0, 1, 2, 3, 4, 5, 9]] 
    read_dist_df = read_dist_df.rename(columns={0: "name", 1: "flag", 2:"ref", 3:"st_pos", 4: "score", 5: "cigar"})
    
    if mapped:    
        read_dist_df = read_dist_df.drop(read_dist_df[read_dist_df['cigar'] == '*'].index)     
        # drop the reads that have less then 30 in quality score
        read_dist_df = read_dist_df.drop(read_dist_df[read_dist_df['score'] < 30].index)
        # get the clipped bases and the length
        read_dist_df['clipped'] = read_dist_df.apply(lambda x: get_clipped_reads(x[9], x['cigar']), axis=1)
        read_dist_df['length']=read_dist_df.apply(lambda x: len(x['clipped']), axis=1)
    else:
        read_dist_df = read_dist_df.loc[read_dist_df['cigar'] == '*']
        read_dist_df['clipped'] = read_dist_df.apply(lambda x: unampped_aligner(read=x[9], ref_seq=ref_seq), axis=1)
        read_dist_df['length']=read_dist_df.apply(lambda x: get_synthesised_size(x['clipped']), axis=1)
        
    read_dist_df = read_dist_df.drop(columns=9)  
    read_dist_df = read_dist_df.sort_values(by=['length', 'score'], ascending=False)
        
    # consider only the _pp reads
    if pp_only:
        read_dist_df = read_dist_df[read_dist_df.groupby('name').name.transform('count')==2]    
   
    if mapped:
        #update reads with the missing/del/insert characters
        read_dist_df['clipped'] = read_dist_df.apply(lambda x: get_final_reads(clipread=x['clipped'], 
                                                           cigar=x['cigar'], 
                                                           st_pos=x['st_pos'], 
                                                           ref_length=len(ref_seq)), axis=1)
    
    #get reads that falls into the synthesis region
    read_dist_df['synths'] = read_dist_df.apply(lambda x: get_synthesised_reads(clipread= x['clipped'],
                                                            ref_seq=ref_seq,
                                                            syn_seq=syn_seq), axis=1)
    read_dist_df['synth_len'] = read_dist_df.apply(lambda x: get_synthesised_size(x['synths']), axis=1)
    return read_dist_df

def plot_read_fraction(
    percent_frame,
    percent_col="percent", 
    count_col="count",
    red_level=5000,
    green_level=10000,
    top_x=20,
    alpha=0.1,
    title="Read count against percentage of reads for the top {} variants",
    ylabel="number of reads",
    xlabel="percentage of reads represented by each variant"
    ):
    """
    Function which plots the number of reads represented by each variant against the percentage of reads with
    that variant. This plot is saved
    Args:
        percent_frame (pandas dataframe): dataframe containing percentage of read counts represented by a variant
            and the number of reads that variant has
        percent_col (str): name of column in dataframe containing percentage data 
            (percentage of reads that have that variant)
        count_col (str): name of column in dataframe containing count data
        green_level (int): minimum number of reads to consider an experiment without issue in
            terms of read number
        red_level (int): number of reads below which to consider the results concerning
        top_x (int): the number of variants to plot (i.e. takes the top x variants representing the highest 
            percentage of reads)
        alpha (float): transparency for background regions
        title (str): title for plot
        ylabel (str): label for y axis
        xlabel (str): label for x axis
    Returns:
        fig: matplot lib figure object"""
    if "{}" in title:
        try:
            title = title.format(top_x)
        except:
            title = title
    # copying and sorting frame by read length
    plot_frame = percent_frame.copy()
    plot_frame.sort_values(by=count_col, inplace=True, ascending=False)
    # generating a plot
    fig, ax = plt.subplots(figsize=(12,8))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # adding the culmulative frequency
    sns.barplot(
        data=plot_frame.iloc[:top_x],
        x=percent_col, y=count_col,
        ax=ax,
        label="Read count of variant",
        color="darkgrey"
    )
    # generating flagging regions of shaded background
    ax.axhspan(
        min(ax.get_ylim()), 
        red_level, color="red", 
        alpha=alpha,
        label="Total reads: Red Warning"
    )
    ax.axhspan(
        red_level, 
        green_level, 
        color="orange", 
        alpha=alpha,
        label="Total reads: Amber Warning"
        )
    ax.axhspan(
        green_level, 
        max(ax.get_ylim()), 
        color="green", 
        alpha=alpha,
        label="Total reads: Green (no warning)")
    # add legend
    ax.legend(bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    
    return fig


def plot_read_dist(
    read_dist_df, 
    ref, 
    col_key='clipped', 
    thresh_reads=[400, 800], 
    run_name='Run10-001', 
    save_path='/mnt/c/Swap/', 
    red_level=1000,
    green_level=10000,
    # expected_length=150
    ):
    """Function which plots the distribution of the reads and saves them and associated CSVs
    ***********Lots of outputs results from this function!!************
    Args:
        df (pandas dataframe): dataframe containing reads distribution results for plotting
            (output of get_result_df)
        ref (str): reference sequence
        col_key (str): column of df containing clipped read sequence
        thresh_reads (list of ints): list of thresholds for read counts to include in variant distribution plots
        run_name (str): name of run
        save_path (str): path to where csvs and plots should be saved
        green_level (int): minimum number of reads to consider an experiment without issue in
        terms of read number
        red_level (int): number of reads below which to consider the results concerning
    Returns:
        variant_fig_output (list): list of tuples. bp and the plot showing a visualisation of reads with deletions/insertions/mismatches relative to
            the reference and the percentage of reads in which they occur
        count_fig_output (list): list of tuples, bp and plot showing the count agains the top 20 variants
"""
    if col_key =='synths':
        # this will restric the reads to the syntheses area with max 5 mutation
        read_dist_df = read_dist_df.loc[read_dist_df['synth_len']>len(ref)-150]
    read_dist_count_df = pd.DataFrame(read_dist_df[col_key].value_counts())
    read_dist_count_df.index.name = col_key
    read_dist_count_df.columns = ['count']
    read_dist_count_df = read_dist_count_df.reset_index(drop=False)
    read_dist_count_df = read_dist_count_df.sort_values(by=['count'], ascending=False) 
    # handle all counts   
    all_reads_count_df = read_dist_count_df.copy()
    totalcount = all_reads_count_df['count'].sum()
    all_reads_count_df['% Reads'] = round((all_reads_count_df['count']/totalcount)*100,2)
    # add in sequence length
    all_reads_count_df["sequence_length"] = [len(re.findall("[a-zA-Z]", seq)) for seq in all_reads_count_df["synths"]]
    all_reads_count_df.to_csv(os.path.join(save_path, "{}_allreads_variant_dist.csv".format(run_name)), sep=',',index=False)

    #Get the number of perfect
    perfect = get_number_perfect(ref,all_reads_count_df,False)
    
    #Use all_reads_count_df to generate the pie chart for the number of deletions from the pie chart function
    piechart = deletion_pie_chart(all_reads_count_df)
    piechart[0].savefig(os.path.join(save_path, "{}_deletions_piechart.png".format(run_name)))
    piechart[1].to_csv(os.path.join(save_path, "{}_deletions_variant_dist.csv".format(run_name)), sep=',',index=False)

    # use all_read_count_df to get culmulative number of reads
    plot_cumsum(
        base_counts=all_reads_count_df,
        len_col="sequence_length",
        red_level=red_level,
        green_level=green_level,
        expected_length=len(ref),
        title="Culmulative sum of reads below a given length: {}".format(run_name),
        savepath=save_path,
        filename="cumulative_graph_{}.png".format(run_name)
    )

    # create flags for variant distibution plots:
    if totalcount < red_level:
        flag = ": Red warning - count number low"
    elif red_level <= totalcount < green_level:
        flag = ": Amber warning - count number low"
    else:
        flag = ""
    # plot variant distribution
    variant_fig_output = []
    count_fig_output = []
    error_fractions = []
    thresholded_reads = []
    thresholded_match_percentages = []
    
    #This for loop uses the threshold given to create the thresholded output including the error rates.
    #Output must be stored as an array if more than one threshold is given (normally 2 are given), otherwise only 
    # the output from the last threshold is passed back
    print('\n')
    print('**************************************')
    print('Calculating the error rates')
    print('**************************************')
    for threshold in thresh_reads:
        # handle counts over threshold
        threshold_count_df = read_dist_count_df.loc[read_dist_count_df['count']>threshold].copy()
        thresholdtotalcount = threshold_count_df['count'].sum()
        threshold_count_df['% Reads'] = round((threshold_count_df['count']/thresholdtotalcount)*100,2)
        threshold_count_df.to_csv(os.path.join(save_path, "{}_{}_threshold_variant_dist.csv".format(run_name,str(threshold))),sep=',',index=False)
        thresholded_reads.append((threshold,thresholdtotalcount))

        percent_matches = get_number_perfect(ref,threshold_count_df,True)
        thresholded_match_percentages.append((threshold,percent_matches))

        #Get Heatmap for each threshold
        heatmapobjecct = heatmap.getheatmap(threshold_count_df,ref,save_path,threshold,run_name)

        #Calculate the error rates for each threshold, output an error csv and return the overall error rate
        all_errors_df, all_errors_fraction = calculate_per_base_error(read_dist_count_df,ref)
        all_errors_df.to_csv(os.path.join(save_path, "{}_all_error_rates.csv".format(run_name)),sep=',',index=False)
        error_df, error_fraction = calculate_per_base_error(threshold_count_df,ref)
        error_fraction = round(error_fraction,4)
        error_df.to_csv(os.path.join(save_path, "{}_{}_threshold_error_rates.csv".format(run_name,str(threshold))),sep=',',index=False)
        print(f"The error rate with a threshold of {threshold} reads is: {error_fraction}")
        error_fractions.append((threshold, error_fraction))

        # plot threshold count dataframe
        y_labs = np.round(threshold_count_df['count']*100/threshold_count_df['count'].sum(),1).astype('str')
        align_df=get_aligner_df(result_df=threshold_count_df, ref_seq=ref, col_key=col_key)
    
        # set color bar definitions
        cmap = colors.ListedColormap(['#ff1a1a', '#000000', '#ff9999', '#b3ffb3', '#e6ff99', '#e6b3ff'])
        bounds = [-1.5, -1.1, -0.5, 0.4, 0.6, 0.8, 1.1]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        sm = plt.cm.ScalarMappable(cmap=cmap)
        variant_fig, ax = plt.subplots(figsize=(15,45))
        cbar = plt.colorbar(sm, ticks=[0.12, 0.27, 0.47, 0.6, 0.78, 0.93],fraction=0.1)
        cbar.ax.set_yticklabels(['D-Insertion','Deletion','M-M-Insertion','Mis-match','P-M-Insertion','Match'], rotation=90,fontsize=18)
        cbar.ax.axes.tick_params(length=0)   
        res = sns.heatmap(align_df.transpose(), ax=ax, cmap=cmap, norm=norm, cbar=False,xticklabels=y_labs+'%',
                            yticklabels=align_df.columns, cbar_kws={"shrink": .82},linewidths=0.1, linecolor='gray')
        res.set_yticklabels(labels=res.get_yticklabels(), va='center',fontsize=18)
        res.set_xticklabels(labels=y_labs+'%',fontsize=18,rotation=90)
        ax.set_ylabel("Synthesised Sequence (3\' -> 5\')", fontsize=20,labelpad=20)        
        plt.title(f'Variant Distribution: {run_name} \n {align_df.shape[0]} Variants - {read_dist_count_df["count"].sum()} reads{flag}', fontsize=20, pad=20)
        ax.set_xlabel('% Fraction of Reads', fontsize=20,labelpad=20)

        variant_fig_output.append((threshold, variant_fig))

        count_fig_output.append(
            (
                threshold,
                plot_read_fraction(
                    threshold_count_df,
                    percent_col="% Reads", 
                    count_col="count",
                    red_level=red_level,
                    green_level=green_level,
                    top_x=20,
                    alpha=0.1,
                    title="Read count against percentage of reads for the top {} variants" + " {} thresh".format(threshold),
                    ylabel="number of reads",
                    xlabel="percentage of reads represented by each variant"
                )
            )
        )
    
    print('\n')
    print(thresholded_match_percentages)
    return variant_fig_output, count_fig_output,piechart,error_fractions,heatmapobjecct,thresholded_reads,perfect,thresholded_match_percentages



def get_all_results(datafile, save_path, ref_seq, syn_seq, run_id, thresh, red_level=1000, green_level=10000):
    """Main Function which will call the df generating function and the plotting
    functions. results will be saved in the given location
    Args:
        datafile (str): path to SAM file
        save_path (str): path to where files are to be saved
        ref_seq (str): reference sequence
        syn_sequence (synthetic sequence)
        thresh: Threshold for the optional threshold for variant inclusion
        run_name (str): name of run number/experimen (run_id)
        green_level (int): minimum number of reads to consider an experiment without issue in
        terms of read number
        red_level (int): number of reads below which to consider the results concerning"""
  
    thresh=int(thresh)
    df=get_result_df(filepath=datafile, syn_seq=syn_seq, ref_seq=ref_seq)
    variant_fig_output, count_fig_output,piechart,error_fractions,heatmapobjecct,thresholded_reads,perfect,thresholded_match_percentages = plot_read_dist(
        read_dist_df=df, 
        col_key='synths', 
        ref=syn_seq, 
        thresh_reads=[thresh, 800], 
        run_name=run_id,
        save_path=save_path,
        red_level=red_level,
        green_level=green_level
    )
    
    # fig2 = plot_read_dist(read_dist_df=df, col_key='synths', ref=syn_seq, thresh_read=800, run_name=run_number,save_path=save_path)

    # save the plots of the aligned reads
    variantfigs = []

    # save the plots of the aligned reads
    for (thresh, fig) in variant_fig_output:
        fig.savefig(os.path.join(save_path, "{}_variant_dist_{}_thresh.png".format(run_id,thresh)))
        pic_IObytes = io.BytesIO()
        fig.savefig(pic_IObytes,  format='png')
        pic_IObytes.seek(0)
        pic_hash = base64.b64encode(pic_IObytes.read())
        variantfigs.append(pic_hash.decode('utf-8'))

    for (thresh, fig) in count_fig_output:
        fig.savefig(os.path.join(save_path, "{}_variant_count_{}_thresh.png".format(run_id,thresh)))
        pic_IObytes = io.BytesIO()
        fig.savefig(pic_IObytes,  format='png')
        pic_IObytes.seek(0)
        pic_hash = base64.b64encode(pic_IObytes.read())
        thresh2b64 = pic_hash.decode('utf-8')

    return variantfigs, piechart[2],error_fractions,heatmapobjecct,thresholded_reads,perfect,thresholded_match_percentages
