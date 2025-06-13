#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
import csv
import os

def plot_and_save_fragment_lengths(filepath, save_path, run_id, read_length=150):
    """Function which plots and saves figures showing the length of insert fragments
    Args:
        filepath (str): path to a CSV of fragment insert lengths, output will be saved 
            in same directory as this
        read_length (int): maximum length of fragment"""
    
    # read in the data
    with open(filepath,"r") as f: 
        # the first entry of this file format is text and so is removed, leaving only data
        a = list(csv.reader(f, delimiter=","))[1:] # list of lists without the first elements which is text

    # trim the data, removing text from each row and cut off at a given read length
    b=[] # This becomes a list of lists, one for each sample
    for x in a:
        b.append([int(i) for i in x[2:read_length+2]])

    # plot the data
    counter = 1
    for fragment in b: # loops through each list in list of lists
        # define x axis
        fragrange= list(range(min(fragment),(len(fragment)))) 
        lookuptable = dict(zip(fragrange,fragment)) # Dictionary for binned data
        topvalue = list(lookuptable.keys())[list(lookuptable.values()).index(max(fragment))] # Get bin with max value
        trimmed_fragment = fragment[:topvalue+20] # Keep bins with max value and the next 20 
        # the above code seems to actually take all bins before the max value and the next 19
        xaxis= list(range(min(trimmed_fragment),(len(trimmed_fragment)))) # Amend X-axis for just bins left
        # plotting
        plt.ylabel("Number of Reads")
        plt.xlabel("Insert Length (bp)")
        plt.title("Insert Lengths Sample "+ str(counter))
        plt.bar(xaxis, trimmed_fragment)
        #plt.show()
        outfile = (os.path.join(save_path, "{}_fragment_length.png".format(run_id)))
        plt.savefig(outfile) 
        plt.clf()
        counter = counter + 1 # Counter for title
