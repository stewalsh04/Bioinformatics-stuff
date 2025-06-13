import pandas as pd # 1.5.3
import numpy as np # 1.24.3
import plotly.graph_objects as go # 5.14.1
from collections import OrderedDict
import os

def align_sequences(ref, seq):
    """Quick function to align sequences to reference - simply adds a * to the reference
    at the index of base insertion. Also right pads the seq to deal with terminal deletions 
    (we really need to fix this bug)
    Args:
        ref (str): reference sequence
        seq (str): sequence
    Returns:
        aligned_ref (str): aligned reference sequence
        aligned_seq (str): aligned sequence"""
    
    ref_list = list(ref)
    insert_indices = [i for i, j in enumerate(seq) if j.islower()]
    
    # in a loop because adding inserts changes the list length...
    for i in insert_indices:
        # using "*" to denote insert
        ref_list.insert(i, "*")
    
    aligned_ref = "".join(ref_list)
    aligned_seq = seq.ljust(len(aligned_ref), "_")
    
    return aligned_ref, aligned_seq
    

def calculate_positional_error(ref_seq, variant_df):
    """Function to calculate fraction of bases synthesised at each position in the reference sequence that
    are errors (deletions, mismatches, insertions) and the fraction of each base (ACTG) present of all bases synthesised
    at each position. 
    Args:
        ref_seq(str): reference sequence
        variant_df (pandas dataframe): dataframe containing synthetic sequences and read counts
    Returns
        errors_df (pandas dataframe): dataframe containign positional errors for reference sequence"""
        
    # check for formatting of reference sequence
    if not ref_seq.isalpha or not ref_seq.isupper():
        raise ValueError("Please check reference sequence, must be all caps: {}".format(ref_seq))

    # store errors in an ordered dict due to needing precise order for plotly
    errors_dict = OrderedDict()
    for i, base in enumerate(ref_seq):
        errors_dict[i] = OrderedDict([
            ("ref_base", base),
            ("deletions", 0),
            ("mismatches", 0),
            ("insertions", 0),
            ("A", 0),
            ("C", 0),
            ("G", 0),
            ("T", 0)
        ])

    # iterating through synthetic sequences, sligning them with the references and calculting errors at
    # each position in the reference sequence
    for synth, count in zip(variant_df["synths"], variant_df["count"]):
        # alignment step
        aligned_ref, aligned_seq = align_sequences(ref_seq, synth)
        # set base index to 0, this advanced by 1 for each non insert character in aligned ref/synth    
        base_index = 0
        # iterate through synth and ref and compare bases
        for ref_base, seq_base in zip(aligned_ref, aligned_seq):
            # add to count of each type of base
            if seq_base.isalpha():
                errors_dict[base_index][seq_base.upper()] += count

            # classify types of mismatches
            if ref_base != seq_base:
                # check for inserts and do not advance base idex if found
                if ref_base == "*": 
                    errors_dict[base_index]["insertions"] += count
                elif seq_base == "_":
                    # check for deletions
                    errors_dict[base_index]["deletions"] += count
                    base_index +=1
                else:
                    # mismatch
                    errors_dict[base_index]["mismatches"] += count
                    base_index +=1
            else: 
                base_index +=1

    # convert to dataframe 
    errors_df = pd.DataFrame.from_dict(errors_dict).T
    # calculate total error count
    errors_df["error_total"] = errors_df[["mismatches", "deletions", "insertions"]].sum(axis=1)
    # extract columns with count data
    count_cols = [col for col in errors_df.columns if col != "ref_base"]
    # convert errors count to fraction
    errors_df[count_cols] = errors_df[count_cols]/variant_df["count"].sum()
    # calculate position from right end of sequence (starts at 1)
    errors_df["Position"] = list(errors_df.index +1)[::-1]
    
    return errors_df

def positional_error_heatmap(errors_df, read_thresh=400):
    """Function to create and interactive heatmap shoing the fraction of base reads at
    each position that are errors. Includes hover data with additional information about bases
    present at each location and the types of errors
    Args:
        errors_df (pandas dataframe): dataframe output of calculate_positional_error - has the errors
            present at each position in the reference sequence
        read_thresh (int): optional, read count threshhold in errors_df
    Returns: 
        fig (plotly graph object): heatmap showing positional errors relative to the reference sequence
            additional error breakdowns are present as hover info"""

    # NB: hover info must be accessed positionally from errors_df
    # initialise figure
    fig = go.Figure()
    # add heatmap. Hover template accessed positionally froms errors df
    fig.add_trace(
        go.Heatmap(
            z=[errors_df["error_total"]],
            x=errors_df.index,
            y=[""],
            colorscale="viridis",
            colorbar=dict(title='Error fraction'),
            customdata = np.dstack([errors_df[col] for col in errors_df.columns]),
            hovertemplate=
                "<b>Base position: %{customdata[9]}</b><br>"
                "<b>Intended base: %{customdata[0]}</b><br><br>" +
                "<b><i>Error fraction: %{customdata[8]:.5f}</i></b><br>"
                "<i>Breakdown (fraction of observed bases):</i><br>"
                "    Deletions: %{customdata[1]:.5f}<br>" +
                "    Mismatches: %{customdata[2]:.5f}<br>" +
                "    Insertions: %{customdata[3]:.5f}<br><br>" +
                "<b><i>Observed base (fraction)</b></i><br>" +
                "    A: %{customdata[4]:.5f}<br>" +
                "    C: %{customdata[5]:.5f}<br>" +
                "    G: %{customdata[6]:.5f}<br>" +
                "    T: %{customdata[7]:.5f}<br>" +
                "<extra></extra> "
        )
    )

    # make a title:
    title = "Plot of errors by position."
    # optionally add read threshold info
    if read_thresh:
        title += " Read threshold {}".format(read_thresh)
    
    fig.update_xaxes(ticktext=errors_df["ref_base"], tickvals=errors_df.index)
    fig.update_layout(title_text=title, title_x=0.5, xaxis_title="synDNA")
    
    return fig

# set up variables

def getheatmap(variant_df, ref_seq, save_path,threshold,run_name):

    # create figures and data frame
    errors_df = calculate_positional_error(ref_seq, variant_df)
    heatfig = positional_error_heatmap(errors_df, threshold)

    # save output
    heatfig.write_html(os.path.join(save_path, "{}_{}_threshold_heatmap.html".format(run_name,str(threshold))))
    errors_df.to_csv(os.path.join(save_path, "{}_{}_threshold_heatmap.csv".format(run_name,str(threshold))),sep=',',index=False)

    return heatfig