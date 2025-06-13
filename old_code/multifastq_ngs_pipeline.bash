#!/bin/bash

######################################################################################
# NGS pipeline intiator script for multi fastq and reference input
# This script will take input from /mnt/c/ngs_fastq
# All fastq files in the correct format in this folder will be sorted
# References and synthetic sequences are taken from ref.txt
# Log files are output to /mnt/c/ngs_fastq/logs
######################################################################################

###################################################################################### 
# Clear all variables, the script is run with source so if the same terminal is used
# variables may be set from a previous run, this creates a problem with the
# sequencing file array, but also clear MULTIRUNID just in case
###################################################################################### 

unset FILES
unset MODE
unset MULTIRUNID
unset THRESH

DIRNAME=$PWD

options=$(getopt -o '' -l threads: -l path: -l id: -l threshold: -l holdinput -l error: -l debug -- "$@") || exit_with_bad_args

#Default arguments
THRESH=400
QUAL=0.03
DEBUG="false"
THREADS=1
HOMEPATH="/mnt/c"
HOLDINPUT="false"
MULTIRUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"

#handling the commandline options

eval set -- "$options"
while true; do
    case "$1" in
        --threads)
            shift
            THREADS="$1"
            ;;
        --path)
            shift
            HOMEPATH="$1"
            ;;
        --threshold)
            shift
            THRESH="$1"
            ;;
        --id)
            shift
            MULTIRUNID="$1"
            ;;
        --error)
            shift
            QUAL="$1"
            ;;
        --holdinput)
           HOLDINPUT="true"
           echo "Holding input fastq files in inputs folder"
           ;;
        --debug)
            DEBUG="true"
        echo "Debugging enabled"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

MODE="multi"
REFERENCEFILE="$HOMEPATH/ngs_pipeline/inputs/ref.txt" # Get reference file
LOGFOLDER="$HOMEPATH/ngs_pipeline/inputs/logs"
echo "Sample ID,Reference,Status,Info" > $LOGFOLDER/$MULTIRUNID.csv
echo $THRESH

cd $HOMEPATH/ngs_pipeline/inputs/

######################################################################################
# Sort out the fastq file into their pairs and put them into an array
# All fastq files should be in the standard iseq format: 
#  Forward fastq: sampleID_L001_R1.fastq.gz
#  Reverse fastq: sampleID_L001_R2.fastq.gz
# L001 stands for Lane 1, the iseq only has 1 lane so this is expected to be present
######################################################################################


declare -A FILES                   # associative array to tie basename with files

for f in *fastq.gz; do                  # search the files with the suffix
    base=${f%_L001_*}                        # remove after "_L001_" To make sample ID the hash key
    if [[ $f == $base* ]] && [[ $f == *"R1"* ]]; then    # if the variable is the current sample ID and is forward
        FILES[$base]=$f                  # then store the filename
    elif [[ $f == $base* ]] && [[ $f == *"R2"* ]]; then # if the variable is the current sample and is reverse
        FILES[$base]+=" $f"
    fi
done

######################################################################################
# Set up the required folders based on the array above
# Move the files to where they need to be for the main NGS pipeline
# If the sample has been previously analyzed the user is told and the pipeline continues
######################################################################################

for base in "${!FILES[@]}"; do           # loop over the hash keys (the sample ID's)
    if [ -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base} ]; then
        echo "Adding to folder of previous analysis for $base, maybe check if it's the same pair of fastq files"
        if [ ! -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq
        fi
    elif [ -d /$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID ]; then
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq
    else
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq
    fi
        # move the fastqs over unles holdinput flag is true
    if [ $HOLDINPUT == "false" ]; then
        mv ${base}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq
    else
        cp ${base}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${base}/fastq
    fi
done


cd $DIRNAME

######################################################################################
# Extract and check the reference file
# Any empty lines will be removed and also end of line charcters /r and /n
######################################################################################

sed -i '/^[[:space:]]*$/d' $REFERENCEFILE # Remove any empty lines from ref file

REFERENCES=($(awk '{print $1}' $REFERENCEFILE))
SYNSEQS=($(awk '{print $2}' $REFERENCEFILE))

######################################################################################
# Initiate the pipeline
# Loop over the sample ID's first and for each sample ID loop through the reference array
# Each individual pair of fastqs are run through ngs_pipeline.bash as a new process with 
#. ./ (This allows transfer of variables to child process)
######################################################################################

echo $DEBUG

function kill_children {
    echo "Forwarding SIGINT to children..."
    jobs -p | xargs kill
    exit 1
}

# ensure that children processes are killed if we receive a SIGINT
trap kill_children SIGINT

#Counter to monitor how many threads are being used
COUNTER=0

for base in "${!FILES[@]}"; do
    for ((i=0; i<${#REFERENCES[@]}; i++)); do
        if [[ $COUNTER <$THREADS ]] &&  [ ! -z "$base" ] ; then
            if [[ $DEBUG == "true" ]]; then
                . ./ngs_pipeline.bash --fastq EVX-$base --ref "${REFERENCES[$i]}" --seq "${SYNSEQS[$i]}" --thresh $THRESH --error $QUAL --debug &
            else 
                . ./ngs_pipeline.bash --fastq EVX-$base --ref "${REFERENCES[$i]}" --seq "${SYNSEQS[$i]}" --thresh $THRESH --error $QUAL &
            fi
            COUNTER=$(( COUNTER + 1 ))
            
            #Maximum number of threads reached, so wait for them to finish then re-set the counter to 0, then it will re-start
            if [ "$COUNTER" -ge "$THREADS" ]; then
                echo "Maximum number of pipelines are running ($THREADS), waiting for them to finish"
                wait
                unset "${COUNTER}"
                echo "Running the next group of pipelines now"
                COUNTER=0
            fi
        fi
    done
done

# wait for all subprocesses to finish
wait

###################################################################################### 
# Clear all variables, the script is run with source so if the same terminal is used
# variables may be set from a previous run, this creates a problem with the
# sequencing file array, but also clear MULTIRUNID just in case
###################################################################################### 

unset FILES
unset MULTIREADID
unset MODE
unset THRESH
unset QUAL

echo "**********Multi NGS $MULTIRUNID Complete!**********"
