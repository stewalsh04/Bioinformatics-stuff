import pandas as pd
import os
import matplotlib.pyplot as plot
import argparse

def process_stats(path):
    files = []
    path = path

    for file in os.listdir(path):
        if file.endswith("numbers.csv"):
            files.append(path + "/" + file)

    read_numbers_df = pd.read_csv(files[0], sep=',', header=None)
    read_numbers_df = read_numbers_df.rename({0: 'Sequence', 1: 'Reads'}, axis='columns')
    read_numbers_df.index = read_numbers_df['Sequence']
    read_numbers_df = read_numbers_df.drop(read_numbers_df.columns[0],axis=1)

    total_reads = sum(read_numbers_df.iloc[:,0])
    read_numbers_df['percentage'] = round((read_numbers_df.iloc[:,0]/total_reads)*100,2)

    read_numbers_df.plot.pie(y='percentage',labels=None,ylabel='',autopct='%1.1f%%', title = 'Proportion of Reads', figsize=(20, 15))
    plot.subplots_adjust(bottom=0.15, wspace=0.1)
    plot.legend(prop={'size': 10},labels=read_numbers_df.index,ncol = 4,bbox_to_anchor=[0.8, 0.01])
    plot.subplots_adjust(bottom=0.15, wspace=0.1)
    plot.savefig(os.path.join(path, "Read_proportions_piechart"))
    read_numbers_df.to_csv(os.path.join(path,'Reads_proportions.csv'))

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument("--path", help="path to stas csv files")
    args=parser.parse_args()

    process_stats(args.path)
