#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import os

def read_in_matches(filepath, remove_first=2):
    """Function to read in matches and convert into usable data
    Args:
        filepath (str): path to a matches.csv file
        remove_first (int): number of items to remove from start of data
    Return:
        (list): data, a list of ints
    """
    with open(filepath,"r") as f: # $dir.matches.csv
        # take in file, keep only the first item and remove the initial items from the list
        raw_data = f.read().splitlines()[remove_first:] # usually has a 0 index post splitlines

    return sorted([int(i) for i in raw_data])


#The following section gets the binned matches as a table
def splitIntoBins(arr, nbins, minval=None, maxval=None):
    """Function to convert binned matches into a matrix
    Args:
        arr (array): array of numbers
        nbins (int): number of bins to use
        minval (int/float): minimum value to use (if none uses minimum of array)
        maxval (int/float): minimum value to use (if none uses maximum of array)
    Returns:
        allbins (array): matrix of binned matches
    """
    minval = min(arr) if minval is None else minval # Select minval if specified, otherwise min of data
    maxval = max(arr) if maxval is None else maxval # Same for maxval
    
    binwidth = (maxval - minval) / nbins # Bin width
    allbins = [[] for _ in range(nbins)] # Pre-make a list-of-lists to hold values

    for elem in arr:
        binnum = int((elem - minval) // binwidth) # Find which bin this element belongs in
        binindex = min(nbins-1, binnum) # To handle the case of elem == maxval
        allbins[binindex].append(elem) # Add this element to the bin
    return np.array(allbins, dtype=object)


def read_bin_and_plot_matches(filepath, run_id, save_path, remove_first=2):
    """Function to read in matches, convert it to usable data, bin matches and turn it into a matrix
    Also saces a histogram of the data and associated .csv files
    Args:
        filepath (str): path to a matches.csv file
        remove_first (int): number of items to remove from start of data"""

    results = read_in_matches(filepath=filepath,remove_first=remove_first)
    mindata = int(min(results))
    maxdata = int(max(results))

    binamounts = splitIntoBins(
        arr=results, 
        nbins=maxdata-mindata, 
        minval=mindata, 
        maxval=maxdata
    ) # Use above function to calculate amount int bins
    binnumbers = np.array(range(mindata,maxdata)) # Get an np array of bin numbers
    # binamountslist = [len(i) for i in binamounts] #convert to list rather than object
    binamountsnp= np.array([len(i) for i in binamounts]) #convert list to np array
    merged = np.vstack((binnumbers, binamountsnp)) #Merge bin numbers and amount in each bin

    #samplename = samplename.replace(".csv", "")
    # saving the merged binned data as a csv
    np.savetxt(os.path.join(save_path, "{}.alignment_matches.csv".format(run_id)), (merged.astype(str)), delimiter="\t",fmt='%s')

    #Get the histogram
    plotwidth = (len(binamounts) / 4)
    if plotwidth < 5:
        plotwidth = 5
    plt.figure(figsize=(plotwidth,4))
    plt.bar(merged[0],merged[1])
    plt.xticks(np.arange(min(merged[0])-1, max(merged[0])+1, 1.0))
    plt.ylabel("Number of Reads")
    plt.xlabel("Bases aligned per Read (Matches)")
    plt.title("Bases aligned per Read for " + run_id)
    plt.tight_layout()

    plt.savefig(os.path.join(save_path, "{}_alignment_matches.png".format(run_id)))
    