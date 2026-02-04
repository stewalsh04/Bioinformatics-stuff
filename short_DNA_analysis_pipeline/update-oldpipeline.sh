#!/bin/sh -x

conda env create -f environment.yml
mkdir /mnt/c/ngs_repository
mkdir /mnt/c/ngs_experiment_data
mkdir /mnt/c/ngs_fastq
cp ref.txt /mnt/c/ngs_fastq
