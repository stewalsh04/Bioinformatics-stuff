#!/usr/bin/env python
#from curses.ascii import NAK
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import seaborn as sns
import variant_distribution
from copy import copy
import base64
from PIL import Image
from io import BytesIO
import io
import report


def read_sequence_data(filepath, sam_data, save_path, thresh, synDNA="REF", run_id="EVX", ref="refseq"):
    """Function to read in reference sequence data
    Args:
        filepath (str): path to a csv contain reference and sequence information
        synDNA (str): synthetic sequence, REF means use whole reference sequence
        sam_data (str): path to .sam data file for vaiant processing
        save_path (str): path to directory to save variant plots into
        run_id (str): id for run    
    ref (str) file path for reference sequence, used for variant distribution"""
    data = pd.read_csv(filepath, header=0, delimiter="\t") 
    # extract reference sequence from data   
    reference = ''.join(map(str, data["ref"].copy().values)) #reference taken from csv in case pysamstats has trimmed it due to no reads
    with open(ref) as f:
        lines = f.readlines()
        fasta_reference=lines[1].upper()
    if reference != fasta_reference:
        print("*****Reference in CSV file does not match fasta reference! Consider checking the CSV!*****")
    # processing synthetic DNA
    synDNA = synDNA.upper()
    if synDNA=="REF":
        synDNA = copy(reference) # replace with whole reference if whole reference sequence to be examined
        synDNAdata = data.iloc[:,2:21].reset_index(drop=True)
        synDNAlength = len(synDNA)
        print("******** The whole reference has been analyzed ********")
        
    else:
        if synDNA==reference:
            synDNAdata = data.iloc[:,2:21].reset_index(drop=True)
            synDNAlength = len(synDNA)
            print("******** The whole reference has been analyzed********")
        else:
            synDNAstart=(reference.find(synDNA))
            if synDNAstart < 0:
                print ("Sequence entered not found in reference!")
                # sys.exit()
                raise ValueError("Sequence entered not found in reference: {}".format(synDNA))
            synDNAlength=len(synDNA)
            synDNAend=synDNAstart+synDNAlength
            synDNAdata = data.iloc[synDNAstart:synDNAend,2:21].reset_index(drop=True) # added .values
   
    synDNAdata["match_percent"] =(synDNAdata.iloc[:,4]/synDNAdata.iloc[:,2])*100
    synDNAdata["seq_pos"] = synDNAdata.index
    bases = ["A", "C", "G", "T"]
    for base in bases:
        synDNAdata[base + "_mismatch"] = [
            j if i!= base else 0 for i, j in zip(synDNAdata["ref"], synDNAdata[base]) 
        ]
    synDNAdata["seq_pos"] = synDNAdata.index
    synDNAdata["mismatch_percent"] = (synDNAdata.iloc[:,6]/synDNAdata.iloc[:,2])*100
    synDNAdata["deletion_percent"] = (synDNAdata.iloc[:,8]/synDNAdata.iloc[:,2])*100
    synDNAdata["insertion_percent"] = (synDNAdata.iloc[:,10]/synDNAdata.iloc[:,2])*100

    synDNAdata = synDNAdata.fillna(0)

    # adding plots
    plotwidth = (len(synDNAdata))
    if plotwidth > 12 and plotwidth < 50:
        plotwidth = plotwidth / 3
    elif plotwidth > 50:
        plotwidth = plotwidth / 10


    stacked = ["match","mismatch","deletion","insertion"]       
    x_pos = np.arange(synDNAlength)

    # Stacked bar, match percent with deletions etc.
    match_fig, match_ax = plt.subplots(1,figsize=(plotwidth,4))
 
    bottom = np.zeros(len(synDNAdata))
    for stacked, color in zip(stacked, ["darkred", "darkorange", "gold","tan"]):
        match_ax.bar(x_pos, synDNAdata[stacked + "_percent"], bottom=bottom, color=color, label=stacked)
        bottom += synDNAdata[stacked + "_percent"].values

    match_ax.legend(bbox_to_anchor =(1.13,0.8), prop={'size': 6})
    match_ax.set_title('Matches')
    match_ax.set(xlabel="Synthesised Sequence", ylabel="Percentage")
    match_ax.set_xticks(x_pos, synDNAdata["ref"])
    match_fig.savefig(os.path.join(save_path, "{}matches.png".format(run_id)))
    
    #Write figure to base64 for the html report
    pic_IObytes3 = io.BytesIO()
    match_fig.savefig(pic_IObytes3,  format='png')
    pic_IObytes3.seek(0)
    pic_hash3 = base64.b64encode(pic_IObytes3.read())
    match64 = pic_hash3.decode('utf-8')

    #Stacked bar for what the mismatches are
    mismatch_fig, mismatch_ax = plt.subplots(1,figsize=(plotwidth,4))
    bottom = np.zeros(len(synDNAdata))
    for base, color in zip(bases, ["skyblue", "purple", "blue", "slateblue"]):
        mismatch_ax.bar(x_pos, synDNAdata[base + "_mismatch"], bottom=bottom, color=color, label=base)
        bottom += synDNAdata[base + "_mismatch"].values
    mismatch_ax.legend(bbox_to_anchor =(1.08, 0.8), prop={'size': 6})
    mismatch_ax.set_title('Mismatches')
    mismatch_ax.set(xlabel="Synthesised Sequence", ylabel="Reads")
    mismatch_ax.set_xticks(x_pos, synDNAdata["ref"])
    mismatch_fig.savefig(os.path.join(save_path, "{}_mismatches.png".format(run_id)))

    pic_IObytes = io.BytesIO()
    mismatch_fig.savefig(pic_IObytes,  format='png')
    pic_IObytes.seek(0)
    pic_hash = base64.b64encode(pic_IObytes.read())
    mismatch64 = pic_hash.decode('utf-8')

    # barplot of total reads
    deletion_fig, deletion_ax = plt.subplots(1,figsize=(plotwidth,4))
    sns.barplot(data=synDNAdata, x="seq_pos", y="reads_pp", color="green", ax=deletion_ax)
    deletion_ax.set_title("Total Number of Reads Aligned")
    deletion_ax.set(xlabel="Synthesised Sequence", ylabel="Reads")
    deletion_ax.set_xticklabels(synDNAdata["ref"])
    deletion_fig.savefig(os.path.join(save_path, "{}_totalreads.png".format(run_id)))

    pic_IObytes2 = io.BytesIO()
    deletion_fig.savefig(pic_IObytes2,  format='png')
    pic_IObytes2.seek(0)
    pic_hash2 = base64.b64encode(pic_IObytes2.read())
    deletion64 = pic_hash2.decode('utf-8')

    synDNAdata.to_csv(os.path.join(save_path,"{}_analyzed.csv".format(run_id)))

    #Call variant distribution to get variant plots (purple ones), and thresholding stats, base64 figures are returned for html report
    variantplots = variant_distribution.get_all_results(
        datafile=sam_data,
        save_path=save_path,
        ref_seq=fasta_reference,
        syn_seq=synDNA,
        run_id=run_id,
        thresh=thresh)

    #Variant plots tuple contains an array of the variant distrubtion plots and a string with the base64 pie chart and other stats from variant distribution
    variantplotsarray = variantplots[0]
    piechart = variantplots[1]
    error_fractions = variantplots[2]
    heatmap = variantplots[3]
    reads_after_thresholding = variantplots[4]
    perfect = variantplots[5]
    perfect_percentages = variantplots[6]

    #Get path of stats text file and csv file to get stats so far and update
    statspath = (os.path.join(save_path,"{}_stats.txt".format(run_id)))
    statscsvpath = (os.path.join(save_path,"{}_stats.csv".format(run_id)))

    #Open csv stats file to get stats so far
    read_dist_df = pd.read_csv(statscsvpath, sep=',', header=None,skiprows=2)
    aligned = read_dist_df.iloc[3:4]

    #Add stats from variant distribution to text stats file
    file = open(statspath, 'a')
    file.write(f'Number of perfect matches identified:\n{perfect}\n')
    for percentages in perfect_percentages:
        file.write(f'Percentage of perfect matches after a threshold of {percentages[0]}:\n{percentages[1]}\n')
    for reads in reads_after_thresholding:
        file.write(f'Reads remaining after a threshold of {reads[0]}:\n{reads[1]}\n')
        reads_lost_after_thresholding = int(aligned[1] - reads[1])
        percentage_lost_after_thresholding = int((reads_lost_after_thresholding / aligned[1])*100)
        file.write(f'Reads lost after a threshold of {reads[0]}:\n{reads_lost_after_thresholding}\n')
        file.write(f'Percentage of reads lost after a threshold of {reads[0]}:\n{percentage_lost_after_thresholding}%\n')
    for errors in error_fractions:
        file.write(f'The error rate with a threshold of {errors[0]} reads is {errors[1]}\n')
    file.close()

    #Add stats from variant distribution to csv stats file
    file = open(statscsvpath, 'a')
    file.write(f'Perfect-matches, {perfect}\n')
    for percentages in perfect_percentages:
        file.write(f'Perfect match percentage {percentages[0]} threshold,{percentages[1]}\n')
    for reads in reads_after_thresholding:
        file.write(f'Reads remaining {reads[0]},{reads[1]}\n')
        reads_lost_after_thresholding = int(aligned[1] - reads[1])
        percentage_lost_after_thresholding = int((reads_lost_after_thresholding / aligned[1])*100)
        file.write(f'Reads lost {reads[0]} threshold,{reads_lost_after_thresholding}\n')
        file.write(f'Percentage lost {reads[0]} threshold,{percentage_lost_after_thresholding}\n')
    for errors in error_fractions:
        file.write(f'errors-{errors[0]} threshold,{errors[1]}\n')
    file.close()

    #Grab data from stats file for the html report
    statsarray= []
    with open(statspath) as my_file:
        for line in my_file:
            statsarray.append(line)

    my_file.close()

    #Generate the html report
    report.getreport(
        statsarray,
        save_path,
        run_id,
        match64,
        mismatch64,
        deletion64,
        str(variantplotsarray[0]),
        str(variantplotsarray[1]),
        piechart,
        heatmap)
