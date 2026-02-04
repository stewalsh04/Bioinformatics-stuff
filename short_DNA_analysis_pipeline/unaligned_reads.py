from Bio import SeqIO, bgzf
import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re
import glob
import sys

outputpath = sys.argv[1]
syn = sys.argv[2]
sample_ID = sys.argv[3]
try: 
    top_vals = sys.argv[4] # currently not specified in pipeline so have used a hacky fix
except:
    top_vals = 15

def gen_reverse_complement(sequence):
    """quick funcation to take a sequence and calculate its reverse complement
    Args:
        sequence (str): sequence
    Returns:
        (str): reverse complement of sequence
    """
    mapping = {"A":"T", "C":"G", "G": "C", "T":"A"}
    
    return "".join([mapping[i] for i in sequence.upper()[::-1]])



def read_zipped_fastq_to_combined_count( 
    forward = "NULL",
    reverse = "NULL",
    top_vals = 15
):
    """Funtion to read in a .gz fastq file and convert it into a dataframe of sequence counts.
    Assumes 5' to 3' direction.
    Args: 
        forward: R1 fastq file
        reverse: R2 fastq file
    Returns:
        count_df (pandas dataframe): dataframe containing sequences and how many times they appear"""

    # make faster! Less memory! See: https://www.biostars.org/p/317524/
    def _process_entry(lines):
        """Function from https://www.biostars.org/p/317524/ to process a fastq file entry
        Args:
            lines: the 4 lines of a fast q file
        Returns: A sequence"""

        if len(lines)!=4:
            raise ValueError(f"Number of lines supplies is not 4, instead is {len(lines)}")

        return "".join(lines[1].splitlines())
    
    sequence_dict = {}
    with gzip.open(forward, "rt") as handle:
        lines = []
        sequence_dict = {}
        for line in handle:
            lines.append(line)
            if len(lines) == 4:
                seq = _process_entry(lines)
                if seq in sequence_dict.keys():
                    sequence_dict[seq] += 1
                else:
                    sequence_dict[seq] = 1
            lines=[]
    with gzip.open(reverse, "rt") as handle:
        lines = []
        sequence_dict = {}
        for line in handle:
            lines.append(line)
            if len(lines) == 4:
                seq = gen_reverse_complement(_process_entry(lines))
                if seq in sequence_dict.keys():
                    sequence_dict[seq] += 1
                else:
                    sequence_dict[seq] = 1
                lines=[]
        
    # plut in a dataframe
    combined_count_df = pd.Series(sequence_dict).to_frame()
    # rename the count column for ease of use
    combined_count_df.columns = ["count"]
    combined_count_df.sort_values(by="count", ascending=False, inplace=True)
    # add order to dataframe
    combined_count_df["order"] = [i+1 for i in range(len(combined_count_df))]
    
    return combined_count_df


def read_zipped_fastq_to_count_old(
    filepath, reverse_complement=False, 
    barcode_five="ACCT", 
    barcode_three=None,
    top_vals=15,
    savepath=outputpath,
):
    """Funtion to read in a .gz fastq file and convert it into a dataframe of sequence counts.
    *****This function is used when treating the forward and reverse fastq files seperately and is not currently used*******
    
    **Note: The function outputs a filtered fastq file which may be useful in the future**

    Assumes 5' to 3' direction. It will then filter to the top X variants and save these as .gz fastq files
    Args: 
        file_path (str): path t fastq file
        reverse_complement (bool): if True, transform sequences to reverse complement
        savepath (str): file path to save filtered fastq.gz into
        top_vals (int): save top X variants as a fastq file
        barcode_five (str): 5' barcode,
        barcode_three (str): 3' barcode, if not none seqence will be restricted between the barcodes five and three
    Returns:
        count_df (pandas dataframe): dataframe containing sequences and how many times they appear"""

    # make faster! Less memory! See: https://www.biostars.org/p/317524/


    with gzip.open(filepath, "rt") as handle:
        raw_fastq = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
        
    sequences = [str(i.seq).upper() for i in raw_fastq.values()]
        
    if reverse_complement:
        sequences = [gen_reverse_complement(i) for i in sequences]
        
    # grab those sequences which are between the 5' and 3' barcode
    if barcode_three:
        # make sure they are uppercase
        barcode_three = barcode_three.upper()
        barcode_five = barcode_five.upper()
        
        filtered_seqs = []
        for seq in sequences:
            try:
                filtered = re.search(r'{}(.*?){}'.format(barcode_five, barcode_three), seq).group(1)
                filtered_seqs.append(seq)
            except:
                continue
        if len(filtered_seqs) > 0: 
            sequences = filtered_seqs
        else:
            raise ValueError('Barcodes not present in any sequences\n5-prime barcode:{}\n3-prime barcode{}'.format(
                barcode_five, barcode_three
                )
            )
        
    # convert to dataframe of count values
    count_df = pd.Series(sequences).value_counts().to_frame()
    # rename the count column for ease of use
    count_df.columns = ["count"]
    # add order to dataframe
    count_df["order"] = [i+1 for i in range(len(count_df))]
    
    # create a filtered fastq file
    if reverse_complement:
        filtered_fastq = [gen_reverse_complement(i) for i in count_df.index[:top_vals]]
    else:
        filtered_fastq = count_df.index[:top_vals]
        
    filtered_fastq_dict = {key:val for key, val in raw_fastq.items() if val.seq.upper() in filtered_fastq}
    
    if "R1" in filepath:
        savename = "R1_filtered.fastq.gz"
    elif "R2" in filepath:
        savename = "R2_filtered.fastq.gz"
    else:
        savename = "unknown_type.fastq.gz"
    
    with bgzf.BgzfWriter(os.path.join(savepath, savename), "wb") as outgz:
        SeqIO.write(sequences=filtered_fastq_dict.values(), handle=outgz, format="fastq")
    
    return count_df

def read_zipped_fastq_to_combined_count_old(
    forward = "NULL",
    reverse = "NULL",
    top_vals = 15
):
    """Funtion to read in a .gz fastq file and convert it into a dataframe of sequence counts.
    Assumes 5' to 3' direction. It will then filter to the top X variants and save these as .gz fastq files
    Args: 
        top_vals (int): save top X variants as a fastq file
        forward: R1 fastq file
        reverse: R2 fastq file
    Returns:
        count_df (pandas dataframe): dataframe containing sequences and how many times they appear"""

    with gzip.open(forward, "rt") as handle:
        forward_raw_fastq = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
    
    with gzip.open(reverse, "rt") as handle:
        reverse_raw_fastq = SeqIO.to_dict(SeqIO.parse(handle, "fastq"))
        
    forward_sequences = [str(i.seq).upper() for i in forward_raw_fastq.values()]
    reverse_sequences = [str(i.seq).upper() for i in reverse_raw_fastq.values()]
        
    reverse_sequences = [gen_reverse_complement(i) for i in reverse_sequences]
        
    # convert to dataframe of count values
    combined_sequences = forward_sequences + reverse_sequences
    combined_count_df = pd.Series(combined_sequences).value_counts().to_frame()

    # rename the count column for ease of use
    combined_count_df.columns = ["count"]

    # add order to dataframe
    combined_count_df["order"] = [i+1 for i in range(len(combined_count_df))]
    
    return combined_count_df

def plot_fastq_variants(
    syn,
    count_df,
    top_vals,
    title = "Plot of variants present in raw data",
    xlabel = "variant",
    ylabel = "count",
    x_var="order",
    y_var="count",
    alpha=0.08,
    red_level=1000,
    green_level=10000,
    savepath=outputpath,
    filename="fastq_variants.png"
):
    """Function to plot the counts of the top variants present in the unaligned data
    Args:
        count_df (pandas dataframe): output of read_zipped_fastq_to_count. This data may have been subject to filtering
            by barcode, or be transformed to reverse complemnt. Data is in descending order of counts
        top_vals (int): number of variants to plot (e.g. 10 highest counts)
        title (str): title for plot
        xlabel (str): name for x axis
        ylabel (str): name for y axis
        x_var (str): column in count_df to use as x axis
        y_var (str): column in count_df to as y axis
        green_level (int): minimum number of reads to consider an experiment without issue in
            terms of read number
        red_level (int): number of reads below which to consider the results concerning
        alpha (float): transparency for background regions
        savepath (str): folder to save image in
        filename (str): name of file
    """

    fig, ax = plt.subplots(figsize=(12,6))
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # adding a bar plot
    sns.barplot(data=count_df[:top_vals], x=x_var, y=y_var, ax=ax, label= "Variant count", color="#1A7E9f")

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

    ax.legend(bbox_to_anchor=(1.1, 1))
    
    out_data = count_df[:top_vals].copy()
    out_data["% total reads"] = out_data["count"]/sum(out_data["count"])*100
    out_data.to_csv(os.path.join(savepath, filename.split(".")[0] + ".csv"))
    
    fastaout = []

    fastaout.append("> reference")
    fastaout.append(syn)

    counter = 1
    for row in out_data.index.values:
        fastaout.append("> " + str(counter))
        fastaout.append(row)
        counter +=1


    fastaoutpd = pd.DataFrame(fastaout)

    fastaoutpd.to_csv(os.path.join(savepath, filename.split(".")[0] + ".fasta"),index=False, header=False)

    # table removed from plot as visually breaks for long sequences
#     table = plt.table(cellText=count_df.values[:top_vals],
#               rowLabels=count_df.index[:top_vals],
#               colLabels=count_df.columns,
#               cellLoc = 'center',
#               rowLoc = 'center',
#               transform=plt.gcf().transFigure,
#               bbox = ([0.4, -0.15, 0.5, 0.3]))

#     table.auto_set_font_size(False)
#     table.set_fontsize(7)

    plt.tight_layout()

    plt.savefig(os.path.join(savepath, filename), bbox_inches='tight')

def read_and_plot_fastqs(
    syn= None,
    input_dir = None,
    savepath=None,
    top_vals = 15,
    barcode_three=None,
    alpha=0.08,
    red_level=1000,
    green_level=10000,
):
    """Function to read in fastq files and plot the counts of the top variants present in the unaligned data. Automatically
    calculates the reverse complement of R2 files
    Args:
        input_dir (str): directory containign fast q files
                savepath (str): folder to save image in
        top_vals (int): number of top variants to plot (e.g. 10 highest counts)
        alpha (float): transparency for background regions
        barcode_three (str)
        green_level (int): minimum number of reads to consider an experiment without issue in
            terms of read number
        red_level (int): number of reads below which to consider the results concerning
        barcode_three (str): 3' barcode, if not none seqence will be restricted between the barcodes five and three
            five is the universal barcode "ACCT" """

    """
    Code below is to find all fsp.fq.gz files, currently these are hard coded for the pipeline
    file_list = glob.glob(os.path.join(input_dir, "*.fsp.fq.gz"))
    title_stem = "top {} fastq variant counts".format(top_vals)
    allreads = []
    """

    #Adding code to merge forward and reverse reads
    forwardfile = os.path.join(input_dir, "R1.fsp.fq.gz")
    reversefile = os.path.join(input_dir, "R2.fsp.fq.gz")

    combined_count = read_zipped_fastq_to_combined_count(
        forward = forwardfile,
        reverse = reversefile,
        top_vals = top_vals
    )

    filename = sample_ID + "_top_" + str(top_vals) + "_reads"

    """
    Code below treats the forward and reverse separately, this is not currently implemented but could be re-activated if needed

    for file in file_list:
        title = title_stem
        # automatically reverse complment r2
        if "R2" in file: 
            reverse_complement = True
            title += sample_ID + ": Reverse Reads"
            filename = sample_ID + "_R2_top15_reads"
        else:
            reverse_complement = False
            title += sample_ID + ": Forward Reads"
            filename = sample_ID + "_R1_top15_reads"

        if barcode_three:
            title += " Between barcodes 5'{} 3'{}".format("ACCT", barcode_three)
            filename += "_{}_{}".format("ACCT", barcode_three)


        count_df = read_zipped_fastq_to_count(
            filepath=file,
            reverse_complement=reverse_complement, 
            barcode_five="ACCT", 
            barcode_three=barcode_three,
            savepath=savepath,
            top_vals=top_vals
        )
        allreads.append(count_df)

        plot_fastq_variants(
            syn=syn,
            count_df=count_df,
            top_vals=top_vals,
            title = title,
            xlabel = "variant",
            ylabel = "count",
            x_var="order",
            y_var="count",
            alpha=alpha,
            red_level=red_level,
            green_level=green_level,
            savepath=savepath,
            filename=filename + ".png"
        )

    """

    plot_fastq_variants(
    syn=syn,
    count_df=combined_count,
    top_vals=top_vals,
    title = filename,
    xlabel = "variant",
    ylabel = "count",
    x_var="order",
    y_var="count",
    alpha=alpha,
    red_level=red_level,
    green_level=green_level,
    savepath=savepath,
    filename=filename + ".png"
    )
read_and_plot_fastqs(syn, outputpath,outputpath,int(top_vals))
