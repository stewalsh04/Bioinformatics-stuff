# README #
### What is this repository for? ###

NGS pipeline to anlayze sequencing data, current version November 2022.


### Setup ###

- Clone repository into Windows

- Navigate into repository

- add jvarkit submodeule using the following git command:
'''git submodule update --init'''

- Run the NGS Pipeline installer:
```
. ./NGSPipeline-Install-WSL.sh (. ./ is correct not a typo)
```
During installation, the installer will install conda, this will prompt you to accept the licence agreement, 
simply type yes or press enter if prompted.


### How to run ###

**Running the Code**

The pipeline can be run in 2 different modes:

Mode 1: Single pair of fastq files against one reference sequence.

Mode 2: Multiple pairs of fastq files against one or multiple reference sequences.

Depending on which mode is run the pipeline has different input requirements, the most efficient way to run the pipeline is to use mode 2.

**Raw data location for modes 1**

Data should be placed into:

C:\ngs_experiment_data (Within windows not WSL)

The data must be in it own folder with a name starting EVX, e.g EVX-Run10.

Each pair of fastq files should then be placed in a folder called inputs.

Example directory structure below:


									ngs_experiment_data  
				
			    EVX-Run10									EVX-Run11

				inputs										inputs
				
	EVX-Run10-FastqR1   EVX-Run10-FastqR2		EVX-Run11-FastqR1	EVX-Run11-FastqR2
	

The pipeline expects a forward and reverse fastq.gz file in each inputs directory with the following suffix:

Forward: *R1_001.fastq.gz

Reverse: *R2_001.fastq.gz

**Raw data location for mode 2**

Simply drop all fastq files into:

C:\ngs_fastq

The fastq files must all share the same group of references sequences.

**Reference file location**

In all modes the reference .fasta file should be placed in the ngs_repositories directory:

C:\ngs_repository (Within windows not WSL)

The reference file should have the following format:

\>ref1

CAACAAGTCGATGCA

For mode 2, a reference text file containing the references to align against is also required, this should have the following format:

<reference-sequence>    /t     <synthetic-sequence>

<reference-sequence>    /t     <synthetic-sequence>

<reference-sequence>    /t     <synthetic-sequence>


This file should already be present in the ngs_fastq folder called ref.txt, simply edit this file.

**Running the Pipeline in mode 1**

The pipeline is run from the terminal in WSL from from the ngs_analysis directory:

- Open up a terminal

- cd ngs_analysis

- Use the following command:

 ```
 source ngs_pipeline.bash --id <Run-ID> --ref <reference> --seq <Synthesised-sequence> --debug
 ```
 
 ``` 
 e.g source ngs_pipeline.bash --id EVX_ID --ref reference --seq ATTTCG --debug
```

if a Synthesised sequence is not supplied the whole reference sequence will be used instead, with no use of 4 base flanking regions

* --id: the name of the sample directory containing the samples (e.g. EVX_run10)
* --ref: the name of the reference fasta file (do no include .fasta file extension)
* --syndna: the synthetic DNA sequence you want to analyse. Optional - leave blank for whole reference
* --debug: use debug mode which leaves temporary files in the working directory. Optional, automatically triggered by pipeline errors


**Running the Pipeline in mode 2**

The pipeline is run from the terminal in WSL from from the ngs_analysis directory:

- Open up a terminal

- cd ngs_analysis

- Simply use the following command:

```
source multifastq_ngs_pipeline.bash 
```

**N.B. Ensure the ref.txt file is up to date in the folder ngs_fastq.**

In all modes the pipeline will create an output directory (name: <ref_id>_outputs) with a folder called "tmp" within it. "tmp" is the working 
directory and is deleted at the end of the pipeline unless in debug mode. The pipeline automatically saves relevant output from this 
folder prior to deletion

**N.B The Run-ID is the directory name in experiment_data where the data is located**


### Running the Example ###

**Mode 1:**

- Move the EVX-Example folder into the ngs_experiment_data folder

- Move the reference sequence example_reference.fasta into the ngs_repository folder

- Run the pipeline with the following command:

```
source ngs_pipeline.bash --id EVX-Example --ref example_reference.fasta --seq ATAATTTCTATCTTGCCCACCCTACTCGACACAGAGCAAAAATCCAACACTCCCAATATTGCCGTGGCTTCGACCTCTTGCTCAGATTTTCTTGTTACCT
```

**Mode 2:**

- Move the two fastq files from the EVX-Example folder into the ngs_fastq folder

- Update the ref.txt file in ngs_fastq so that it is the same as reflist-example.txt 

- copy and rename the example fasta file in the repository to M13-bc6

- Run the pipeline with the following command:

```
source multifastq_ngs_pipeline.bash
```

### Threshold for Variant Distribution plot ###

In both modes the threshold for the variant distribution plot can be set with the additional flag:

```
--threshold <200>
```

There is also the option to use an auto threshold, this is activated with the following flags:

```
--threshold AUTO --error <0.003>
```

The sequencing error rate must also be set in decimal form, this is the error rate for the sequencing run the sample originated from, the error rate can be obtained from the Illumina sequence analysis viewer (SAV).

**The threshold is calculated as follows:**

No. aligned reads x error rate

So if a sample had 200,000 aligned reads and the error rate was 1.2%, the threshold would be calculated as follows:

200,000 x 0.012 = 2,400

### Log Files ###

Every pipeline run in both modes will generate a log file, in mode 1 the log files are located in the log folder in the ngs_analysis folder, in mode 2 the log files are located in the log folder in the ngs_fastq folder.

### Re-running the same data ###

In all modes if the same fastq files are used as input with the same reference, additional folders called SYN_1, SYN_2, SYN_3……. will be created in the output folder.

### Running in Phi-X mode ###

When the sequencer is run, the PhiX genome is spiked in, this gives an opportunity to check the run quality against a
known good sample. If aligment to the PhiX genome is good, the sequencer has likely run well. The PhiX reads are in
the undertermined reads, to use PhiX mode create a directory in the run folder called 'Phix' and place both the R1
and R2 undetermined fastq.gz files inside. Then to activate PhiX mode set the reference to 'phix' as shown below and run in the cloned repository:

source ngs_pipeline.bash --id <Run-ID> --ref phix

The undetermined read file should have the following file name structure:

undetermined<runID>R1_001.fastq.gz
undetermined<runID>R2_001.fastq.gz