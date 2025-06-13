import pandas as pd
import os
import matplotlib.pyplot as plot
import argparse

"""
This script will take the stats.csv files from a number of NGS samples that have been run throught the pipeline
and summarise the statistics into 3 summary plots:
Read numbers
Match percentages
Percentages of reads lost by thresholding
Author Steve Walsh 
"""

def process_stats(path,mode):
    files = []
    path = path
    mode=mode

    df=[]

    for file in os.listdir(path):
        if file.endswith("stats.csv"):
            files.append(path + '/' + file)

    for file in files:
        read_dist_df = pd.read_csv(file, sep=',', header=None)
        df.append(read_dist_df.iloc[:,1])

    df = pd.DataFrame(df)

    df.columns=['Run ID',
                'Seq ID',
                'fastqR1',
                'fastqR2',
                'Total Reads',
                'Aligned Reads',
                'Perfect matches',
                'Perfect match percentage 1% threshold',
                'Perfect match percentage 800 threshold',
                'Reads remaining 1% threshold',
                'Reads lost 1% threshold',
                'Percentage lost 1% threshold',
                'Reads remaining 800 threshold',
                'Reads lost 800 threshold',
                'Percentage lost 800 threshold',
                'errors-5156 threshold',
                'errors-800 threshold']

    if mode == "multi":
        df.index = df['Seq ID']
    else:
        df.index = df['Run ID']
    
    df = df.drop(df.columns[0],axis=1)
    df = df.drop(df.columns[0],axis=1)
    df = df.apply(pd.to_numeric)

    df_reads = df[['Total Reads','Aligned Reads','Perfect matches',
                'Reads remaining 1% threshold',
                'Reads remaining 800 threshold']].copy()

    df_percentage_match = df[['Perfect match percentage 1% threshold','Perfect match percentage 800 threshold']]

    df_percentage_lost = df[['Percentage lost 1% threshold','Percentage lost 800 threshold']]

    df_percentage_match.plot.line(linewidth = 2,ylabel="Percentage",figsize=(15, 8),title='Percentage Matches')
    plot.legend(bbox_to_anchor =(0.7, -0.11), prop={'size': 9})
    plot.subplots_adjust(bottom=0.15, wspace=0.1)
    plot.savefig(os.path.join(path, "Percentage_matches.png"))
    df_reads.plot.line(linewidth = 2,ylabel="Read Count",figsize=(15, 8),title="Sequencing Reads")
    plot.legend(bbox_to_anchor =(0.6, -0.11), prop={'size': 8})
    plot.subplots_adjust(bottom=0.2, wspace=0.15)
    plot.savefig(os.path.join(path, "Sequencing Reads.png"))
    df_percentage_lost.plot.line(linewidth = 2,ylabel="Percentage",figsize=(15, 8),title="Percentage Lost")
    plot.legend(bbox_to_anchor =(0.7, -0.11), prop={'size': 9})
    plot.subplots_adjust(bottom=0.15, wspace=0.1)
    plot.savefig(os.path.join(path, "Percentage_lost.png"))

if __name__ == '__main__':

    parser=argparse.ArgumentParser()

    parser.add_argument("--path", help="--path= path to stas csv files")
    parser.add_argument("--mode", help="--multi= for pipeline use only")
    args=parser.parse_args()
    mode = args.mode

    if mode == "true":
        mode = "multi"
    else:
        mode = "single"

    process_stats(args.path, mode)
