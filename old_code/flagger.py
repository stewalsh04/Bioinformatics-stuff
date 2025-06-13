"""
This script handles functions which allow automatic flagging of NGS pipeline runs
This defines whether a given NGS sample is pass/fail or soft pass/soft fail

"""
import pandas as pd
import os
import argparse

def extact_perfect_percentages(name, input_dir, syndna, approx_id=0.9):
    """Function to automatically flag pass, soft pass, soft fail and fail using the
    CSV output of the variant distribution script
    Args:
        name (str): name of experiment (derived from fastq file)
        input_dir (str): name of directory containing variant distribution csvs
        syndna (str): synthetic dna
        approx_id: threshold for minimum approximate sequence id - i.e. soft fail
            a sequence if there are no perfect sequences but there is a near perfect
            sequence present
    """


    
    # thresh_400 = pd.read_csv(os.path.join(input_dir, f"{name}_400_threshold_error_rates.csv"))
    thresh_800 = pd.read_csv(os.path.join(input_dir, f"{name}_800_threshold_error_rates.csv"))
    all_reads = pd.read_csv(os.path.join(input_dir, f"{name}_all_error_rates.csv"))
    thresh_800["approx_id"] = [(len(syndna) - i)/len(syndna) for i in thresh_800["error_sum"]]

    try:
        percent_800 = 100*thresh_800["count"][thresh_800["total_base_errors"] == 0].sum()/thresh_800["count"].sum()
    except:
        percent_800 = 0
    try:
        percent_all = 100*all_reads["count"][all_reads["total_base_errors"] == 0].sum()/all_reads["count"].sum()
    except: 
        percent_all = 0

    if percent_800 >= 0.25:
        flag = "pass"
    elif 0.05 <= percent_800 < 0.25:
        flag = "soft pass"
    elif percent_all > 0:
        flag = "soft fail"
    elif thresh_800["approx_id"].max() >= approx_id:
        flag = "soft fail"
    else:
        flag = "fail"


    with open(os.path.join(input_dir, "flags.txt"), "w") as f:
        f.write(f'{name}\n')
        f.write(f'percent perfect at 800 threshold :\t{percent_800}\n')
        f.write(f'percent perfect at all reads :\t{percent_all}\n')
        f.write(f'Max approximate identity to reference at 800 threshold (threshold for soft fail: {approx_id}):\t{thresh_800["approx_id"].max()}\n')
        f.write(f"flag: {flag}. Nb manually check fails and soft fails")


if __name__ == '__main__':
    
    parser=argparse.ArgumentParser()

    parser.add_argument("--input_dir", help="to variant files")
    parser.add_argument("--syndna", help="synthetic dna sequence")
    parser.add_argument("--name", help="Experiment name")
    args=parser.parse_args()

    extact_perfect_percentages(name=args.name, input_dir=args.input_dir, syndna=args.syndna)
