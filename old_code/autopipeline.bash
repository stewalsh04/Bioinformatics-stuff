##!/bin/bash

######################################################################################
# NGS pipeline intiator script for automated NGS pipeline running
# This script will take input from /mnt/c/ngs_fastq/input
# All fastq files in the correct format in this folder will be sorted
# synthetic DNA is processed from config.xlsx int he same folder as the fastq files
# Log files are output to /mnt/c/ngs_fastq/logs
######################################################################################

###################################################################################### 
# Clear all variables, the script is run with source so if the same terminal is used
# variables may be set from a previous run

# usage:
# source autopipeline.bash --folder input
# where input is the name of the subdirectory within ngs_fastq containing the relevant
# config and fastqs
###################################################################################### 

unset LINE
unset THRESH

DIRNAME=$PWD

options=$(getopt -o '' -l threads: -l path: -l id: -l threshold: -l holdinput -l error: -l debug -l phix -l noqc -l multi-site -l summary-stats -l unaligned-thresh: -- "$@") || exit_with_bad_args

#Default arguments
THRESH="AUTO"
QUAL=0.03
DEBUG="false"
THREADS=1
HOMEPATH="/mnt/c"
MULTIRUNID="PipelineRun-$(date '+%Y-%m-%d-%R')"
DATEID=$MULTIRUNID
HOLDINPUT="false"
PHIX="false"
QC="true"
MULTI_SITE="false"
SUMMARY_STATS="false"
SEED=20 # This is the seed sequence for bwa-mem, the default is 20
UNALIGNEDTHRESH=15

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
        --id)
            shift
            MULTIRUNID="$1"
            ;;
        --threshold)
            shift
            THRESH="$1"
            ;;
        --holdinput)
           HOLDINPUT="true"
           echo "Holding input fastq files in inputs folder"
           ;;
        --error)
            shift
            QUAL="$1"
            ;;
        --debug)
            DEBUG="true"
        echo "Debugging enabled"
            ;;
        --phix)
            PHIX='true'
        echo "PhiX scanned enabled"
            ;;
        --noqc)
            QC='false'
        echo "QC checks switched off" # This is only applicable if the pipeline in in multi-site mode
            ;;
        --multi-site)
            MULTI_SITE='true'
            POOL_NAME='Pooled-summary'
        echo "Pipeline in multi site mode"
        echo "****Have you set a seed sequence? Default is 20****"
        sleep 2
            ;;
        --summary-stats)
            SUMMARY_STATS='true'
            POOL_NAME='Summary-stats'
        echo "Summarising stats for all samples"
            ;;
        --unaligned-thresh)
            shift
            UNALIGNEDTHRESH="$1"
         echo "Setting unaligned threshold to $UNALIGNEDTHRESH"
            ;;
         --)
            shift
            break
            ;;
    esac
    shift
done

sleep 3
DIRNAME=$PWD
RUNFOLDER="$HOMEPATH/ngs_pipeline/inputs/"
CONFIGFILE="$RUNFOLDER/config.xlsx" # Get reference file
LOGFOLDER="$HOMEPATH/ngs_pipeline/inputs/logs"
REPOSITORY="$HOMEPATH/ngs_pipeline/references"
POOLSTATS="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/$POOL_NAME/"
echo "Sample ID,Reference,Status,Info" > $LOGFOLDER/$MULTIRUNID.csv
echo $THRESH

if [ ! -f $CONFIGFILE ]; then
    echo "Could not find the configuration file confit.xlxs"
    echo "Exiting now"
    sleep 5
    exit 127
fi

######################################################################################
# process configuration file into ref.txt (for running the pipeline)
######################################################################################

conda activate NGS-Pipeline    

echo "*******************************************"
echo "Checking config file for references, invalid synDNA will not be run"
echo "*******************************************"
sleep 2
python automation_utils.py --config $CONFIGFILE --repo $REPOSITORY --fastqs $RUNFOLDER
REFERENCEFILE="$RUNFOLDER/ref.txt"

if [ $HOLDINPUT == "false"  ] && [ $PHIX == "false" ]; then
    mv $CONFIGFILE $LOGFOLDER/configfiles/$MULTIRUNID.config.xlsx
fi

CONFIG_MODE=$(head -n 1 $RUNFOLDER/mode.txt)
CONFIG_SEED=$(tail -n 1 $RUNFOLDER/mode.txt)

if [ $CONFIG_MODE == 'multi-site' ]; then
    MULTI_SITE='true'
    echo "Pipeline in multi site mode"
fi

SEED=$CONFIG_SEED

rm $RUNFOLDER/mode.txt

if [ $MULTI_SITE == 'true' ]; then
    echo "Seed taken from configuration as $SEED"
fi

#When multi-site is switched on summary stats are switched on by default and handled 
# differently by the multi-site script, so must be turned off here if switched on.
if [ $MULTI_SITE == 'true' ] && [ $SUMMARY_STATS == 'true ' ]; then
    echo "***Summary stats are switched on by default when running in multi site mode***"
    sleep 2
    SUMMARY_STATS='false'
    POOL_NAME='Pooled-summary' # I differentiate between summary stats across samples and within a pool
fi

######################################################################################
# Extract and check the reference file
# Any empty lines will be removed and also end of line charcters /r and /n
######################################################################################

sed -i '/^[[:space:]]*$/d' $REFERENCEFILE # Remove any empty lines from ref file

######################################################################################
#**************PhiX Check***************
#This carries out the PhiX check if the PhiX flag is used, this aligns the undetermined reads against the PhiX genome, a small amount of PhiX is spiked into each
#sequencing run as a QC check for the sequencing process. The correct undetermined reads files must be present in the inputs folder along with the sequencing files
#from THE SAME sequencing run. Undetermined reads CANNOT be mixed between sequening runs. 
#The PhiX check is carried out prior to processing the remaining sequeincing samples and is handled independent of the config input file, it is therefore outside
#the while loop.
######################################################################################

#Handle the directories for PhiX
if [[ $PHIX == "true" ]]; then 
    echo "Running PhiX error check"
    sleep 2
    if [ -d /$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined ]; then
        echo "PhiX check may have already been carried out for this sequencing run"
        echo "The previous PhiX check will be overwritten"
	sleep 2
	#Remove previous undetermined files to re-set analysis
	if  [ -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq ]; then
	    rm $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq/*
	fi
        if [ ! -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
        fi
    elif [ -d /$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID ]; then
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
    else
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
    fi

    #Handle the PhiX files
    if [[ $HOLDINPT == "false" ]]; then
        mv ${RUNFOLDER}/Undetermined_S0_L001_R1_001.fastq.gz $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
        mv ${RUNFOLDER}/Undetermined_S0_L001_R2_001.fastq.gz $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
    elif [[ $HOLDINPUT == "true" ]]; then
        cp ${RUNFOLDER}/Undetermined_S0_L001_R1_001.fastq.gz $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
        cp ${RUNFOLDER}/Undetermined_S0_L001_R2_001.fastq.gz $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Undetermined/fastq
    fi

    . ./ngs_pipeline.bash --fastq Undetermined --ref phix
fi

COUNTER=0
QCCOUNTER=0
while IFS= read -r LINE
do
    # split line into array using tab delimitator - 0: sample 1: reference fasta 2: synDNA
    ARRAYLINE=(${LINE// / })
    SAMPLE=${ARRAYLINE[0]}
    REFFASTA=${ARRAYLINE[1]}
    SYNDNA=${ARRAYLINE[2]}

    # ######################################################################################
    # # Set up the required folders based on the sample names above
    # # Move the files to where they need to be for the main NGS pipeline
    # # If the sample has been previously analyzed the user is told and the pipeline continues in 
    # # the same folder
    # ######################################################################################

    if [ -d /$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE} ]; then
        if [ $SUMMARY_STATS == "true" ] && [ ! -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/summary-stats ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/summary-stats
        fi
        echo "Adding to folder of previous analysis for $SAMPLE"
        if [ ! -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
        fi
    elif [ -d /$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID ]; then
        if [ $SUMMARY_STATS == "true" ] && [ ! -d $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/summary-stats ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/summary-stats
        fi
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
        #sleep 5
    else
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID
        if [ $SUMMARY_STATS == "true" ]; then
            mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/summary-stats
        fi
        #If running in multi site mode make directory to gather stats on different sequences, this should only be made when the parent directory is created
        if [ $MULTI_SITE == 'true' ]; then
            mkdir $POOLSTATS
        fi
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}
        mkdir $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
    fi

    # move the fastqs over to working output directory unless holdinput flag is true or the files are already there
    if [ $HOLDINPUT == "false" ] && [ ! -f $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/${SAMPLE}_L001_R1_001.fastq.gz ]; then
	mv ${RUNFOLDER}/${SAMPLE}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
    elif [ $HOLDINPUT == "false" ] && [ ! -f $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/${SAMPLE}_L001_R2_001.fastq.gz ]; then
        mv ${RUNFOLDER}/${SAMPLE}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
    elif [ $HOLDINPUT == "true"  ] && [ ! -f $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/${SAMPLE}_L001_R1_001.fastq.gz ]; then
        cp ${RUNFOLDER}/${SAMPLE}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
    elif [ $HOLDINPUT == "true"  ] && [ ! -f $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/${SAMPLE}_L001_R1_001.fastq.gz ]; then
        cp ${RUNFOLDER}/${SAMPLE}_* $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq
    fi

    ######################################################################################
    # Carry out the QC for one pair of fastq files that contian the pool of different DNA sequences
    ######################################################################################

    function create_folder {
        mkdir -p "$1"
        chmod 700 "$1"
    }

    QC_DIR="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/QC"

    if [[ $QCCOUNTER < 1 ]] && [ $QC == 'true' ] && [ $MULTI_SITE == 'true' ]; then

        echo "**************************************************************"
        echo "Carrying out QC for the single pair of fastq files"
        echo "**************************************************************"

        #Check if QC directory is there and create one if needed
        if [[ -d $QC_DIR ]]; then
            echo "QC directory already present, overwriting previous QC checks"
            sleep 2
        else
            create_folder "$QC_DIR"
        fi

        fastqc $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/*.fastq.gz --outdir $QC_DIR

        fastp --thread 16 --trim_poly_x --in1 $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/*R1*.fastq.gz \
        --in2 $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/EVX-${SAMPLE}/fastq/*R2*.fastq.gz --out1 $QC_DIR/R1.fsp.fq.gz \
        --out2 $QC_DIR/R2.fsp.fq.gz --json $QC_DIR/${SAMPLE}_fastp.json \
        --html $QC_DIR/${SAMPLE}_fsp_data.html \
        --report_title ${SAMPLE}_data_fsp 1 > $QC_DIR/${SAMPLE}_fsp.log

        # concatenate json files into one .csv
        cat $QC_DIR/*json | head -41 | tail -1 >> $QC_DIR/${SAMPLE}_insert_size.csv

        ######################################################################################
        # Carry out the unaligned analysis on just the one pair of fastq files
        ######################################################################################
        
        echo "**************************************************************"
        echo "Carrying out Unaligned-analysis for the single pair of fastq files"
        echo "**************************************************************"

        UNALIGNED_DIR="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Unaligned_Analysis"

        #Check if Un-aligned analysis directory is there and create one if needed
        if [[ -d $UNALIGNED_DIR ]]; then
            echo "Un-aligned analysis directory already present, overwriting previous Un-aligned analysis checks"
            sleep 2
        else        
            create_folder "$UNALIGNED_DIR"
        fi

        python unaligned_reads.py $QC_DIR $SYNDNA EVX-${SAMPLE} $UNALIGNEDTHRESH

        mv $QC_DIR/*_top_${UNALIGNEDTHRESH}_reads* $UNALIGNED_DIR

        echo "*******************************************"
        echo "Carrying out multiple alignment with Muscle for unaligned reads"
        echo "*******************************************"
        #Multiple aligner Muscle https://www.drive5.com/muscle/
        muscle -in $UNALIGNED_DIR/*_top_${UNALIGNEDTHRESH}_reads.fasta -out $UNALIGNED_DIR/EVX-${SAMPLE}_top_${UNALIGNEDTHRESH}_reads.aligned.fasta.tmp

        #Re-order the aligment so the most abundant variants are at the top (Using amended python script from Muscle author)
        python stable.py $UNALIGNED_DIR/*_top_${UNALIGNEDTHRESH}_reads.fasta $UNALIGNED_DIR/EVX-${SAMPLE}_top_${UNALIGNEDTHRESH}_reads.aligned.fasta.tmp > $UNALIGNED_DIR/EVX-${SAMPLE}_top_$UNALIGNEDTHRESH.aligned.fasta

        #Remove tempory alignment files
        rm $UNALIGNED_DIR/EVX-${SAMPLE}_top_{$UNALIGNEDTHRESH}_reads.aligned.fasta.tmp

        #Make the multiple alignment plots using pymsaviz (https://pypi.org/project/pymsaviz/)
        python unaligned_reads_plot.py $UNALIGNED_DIR/EVX-${SAMPLE}_top_{$UNALIGNEDTHRESH}_reads.aligned.fasta EVX-${SAMPLE} $UNALIGNED_DIR

        QCCOUNTER=$(( QCCOUNTER + 1 ))    
    fi

    # ######################################################################################
    # # Initiate the pipeline
    # ######################################################################################

    function kill_children {
        echo "Forwarding SIGINT to children..."
        jobs -p | xargs kill
        exit 1
    }

    # ensure that children processes are killed if we receive a SIGINT
    trap kill_children SIGINT

        if [[ $DEBUG == "true" ]]; then
            . ./ngs_pipeline.bash --fastq EVX-$SAMPLE --ref $REFFASTA --seq $SYNDNA --thresh $THRESH --eror $QUAL --debug &
        elif  [ $MULTI_SITE == 'true' ]; then
            . ./ngs_pipeline_multi-site.bash --fastq EVX-$SAMPLE --ref $REFFASTA --seq $SYNDNA --thresh $THRESH --error $QUAL &
        else
            . ./ngs_pipeline.bash --fastq EVX-$SAMPLE --ref $REFFASTA --seq $SYNDNA --thresh $THRESH --error $QUAL &
        fi

    #Maximum number of threads reached, so wait for them to finish then re-set the counter to 0, then it will re-start
    COUNTER=$(( COUNTER + 1 ))
    echo $COUNTER
    if [ "$COUNTER" -ge "$THREADS" ]; then
        echo "Maximum number of pipelines are running ($THREADS), waiting for them to finish"
        wait
        unset COUNTER
        echo "Running the next group of pipelines now"
        COUNTER=0
    fi

done < "$REFERENCEFILE"

wait
if [ $SUMMARY_STATS == "true" ] || [ $MULTI_SITE == 'true' ]; then
    STATS_DIR="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/$POOL_NAME/Summary_Stats_${DATEID}"
    mkdir $STATS_DIR
    mv $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/$POOL_NAME/*.csv $STATS_DIR/
    python run_summary.py --path $STATS_DIR  --mode $MULTI_SITE
    cp $CONFIGFILE $STATS_DIR/
fi

if [ $MULTI_SITE == "true" ]; then
    python multi_site_summary_stats.py --path $STATS_DIR
    rm $STATS_DIR/*numbers.csv
fi

#Remove any remaining fastq file from input unless holdinput is true
if [ $HOLDINPUT == "false"  ] && [ -f ${RUNFOLDER}/${SAMPLE}_L001_R1_001.fastq.gz ]; then
    rm ${RUNFOLDER}/${SAMPLE}_*
fi

# ###################################################################################### 
# # Clear all variables, the script is run with source so if the same terminal is used
# ###################################################################################### 
unset THRESH
unset QUAL
unset LINE
unset SAMPLE
unset REFASTA
unset SYNDNA
unset QCCOUNTER
unset MULTI_SITE
echo "**********Pipeline Complete! All samples run**********"
conda deactivate
