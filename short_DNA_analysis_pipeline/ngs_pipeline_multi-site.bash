##]#!/bin/bash

#################################################################################################
#This is the main aligning pipeline, my current understanding of what it does is as follows:	#
# - Pre filter the .fastq with fastp															#
# - Align with burrows wheeler aligner (bwa) an aligner suited to the type of data we have		#
# - Filer using samtools for low quality, and then aligned or not aligned						#
# - It outputs a number of .csv files I'm still figuring out									#
# - It uses a bioinformatic java package called jvarkit to pull out aligned reads				#
# - It uses python package pysamstats to extract stats from the aligned reads	

#################################################################################################

# example command:
# source ngs_pipeline.bash --id EVX_ID --ref reference.fasta --syndna ATTTCG --debug

conda activate NGS-Pipeline

#################################################################################################### 
# Check if variable MULTIRUNID is available from parent multi input bash script
# If it is there the pipeline is being run in multi mode
# If it's not there the pipeline is in single mode so set it
# MODE defaults to multi unless the if statment changes it to single
#################################################################################################### 

####################################################################################################
# commandline argument handling
####################################################################################################


function exit_with_bad_args {
    echo "Usage: source ngs_pipeline.bash --id <EV_ID> --ref <REF_ID> [--seq <SYNDNA>] [--thresh <THRESHOLD>] [--error <Error Rate>] [--debug]"
    echo "Invalid arguments provided" >&2
    return # this stops the terminal closing when run as source
}

options=$(getopt -o '' -l id: -l fastq: -l ref: -l seq: -l thresh: -l error: -l debug -- "$@") || exit_with_bad_args

# default arguments for debug and threshold
DEBUG="false"
THRESH="400"
QUAL=0.003

echo $MULTIRUNID
if [[ "$MODE" != "multi" ]]; then
    LOGFOLDER="$HOMEPATH/ngs_pipeline/inputs/logs"
    echo "Sample ID,Reference,Status,Info" > $LOGFOLDER/$MULTIRUNID.csv
    MODE="single"
fi

echo "NGS Pipeline is being run in $MODE mode"

# handling the commandline options
eval set -- "$options"
while true; do
    case "$1" in
        --id)
            shift
            MULTIRUNID="$1"
            ;;
        --fastq)
            shift
            EVX_RUNID="$1"
            ;;
        --ref)
            shift
            REF_ID="$1"
            REFASTA="$1"
            ;;
        --seq)
            shift
            SYNDNA="$1"
            ;;
        --thresh)
           shift
           THRESH="$1"
           ;;
        --error)
           shift
           QUAL="$1"
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

# check for required args
if [ -z "$EVX_RUNID" ] || [ -z "$REFASTA" ]; then
    echo "--fastq and --ref  require arguments"
    exit_with_bad_args
fi

#Check for SYN DNA if empty use the whole refernce
if [ -z "$SYNDNA" ]; then
    SYNDNA=$"REF"
    echo "WARNING: No Synthetic sequence supplied. Using whole reference for $EVX_RUNID with $REF_ID"
fi

refID=$REFASTA # Set ID for reference used in log file and output folder (Without .fasta)

#Check if the reference has the correct suffix .fasta for the pipeline
if [[ $REFASTA != *.fasta ]]; then
    REFASTA+=".fasta"

#For the log file remove .fasta if present (Think it looks tidier)
elif [[ $REFASTA == *.fasta ]]; then
    refID="${REFASTA%%.*}"
fi

if [[ $THRESH == "AUTO" ]]; then
    THRESHMODE="AUTO"
fi

LOGFILE=$LOGFOLDER/$MULTIRUNID.csv

echo "Run ID: $MULTIRUNID"
echo "Input data: $EVX_RUNID"
echo "Reference fasta: $REFASTA"
echo "Synthetic sequence: $SYNDNA"
echo "Debug mode: $DEBUG"


####################################################################################################
# directory handling
# includes options for deleting old folders
####################################################################################################

EXP_DIR="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/$EVX_RUNID/fastq"
OUT_DIR="$HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/$EVX_RUNID/${refID}_outputs"
TMP_DIR="$OUT_DIR/tmp"
REF_FILE="$HOMEPATH/ngs_pipeline/references/$REFASTA"

#################################################################################################### 
#check reference sequence, input directory and input fastq files exist
#Check if the fastq files have previously been analyzed against reference
#Make requires folders
#################################################################################################### 

#Check for input directory
if [ ! -d $EXP_DIR ]; then
    echo "Error!! Input directory not found: $EXP_DIR"
    conda deactivate
    return # file not found error
fi


#Check for fastq files:
if [ ! -f $EXP_DIR/*R1*.fastq.gz ]; then
    echo "ERROR!! Forward fastq file not found for $EVX_RUNID with $REF_ID"
    sleep 3
    conda deactivate
    echo "If running a multi analysis the pipeline will resume shortly"
    echo "$EVX_RUNID,$refID,FAIL,Forward fastq not found" >> $LOGFILE
    return
elif [ ! -f $EXP_DIR/*R2*.fastq.gz ]; then
    echo "Error!! Reverse fastq file not found"
    sleep 3
    conda deactivate
    echo "If running a multi analysis the pipeline will resume shortly"
    echo "$EVX_RUNID,$refID,FAIL,Reverse fastq not found" >> $LOGFILE
    return
fi

#Check for reference sequence:
if [ ! -f $REF_FILE ] && [ "$REF_ID" != "phix" ]; then
    echo "ERROR!! Reference file not found at $REF_FILE: $REFASTA"
    echo "If running a multi analysis the pipeline will resume shortly"
    echo "$EVX_RUNID,$refID,FAIL,Reference file not found" >> $LOGFILE
    conda deactivate
    sleep 5
    return  
fi

#Check for previous analysis
if [ -d $OUT_DIR ] && [ "$REF_ID" != "phix" ]; then
    echo "Previous analysis found for $EVX_RUNID with $refID. Adding to folder for $EVX_RUNID"   
fi

# make directories
function create_folder {
    mkdir -p "$1"
    chmod 700 "$1"
}
create_folder "$OUT_DIR"
create_folder "$TMP_DIR"


#####################################################################################################
# Phix
#####################################################################################################
if [ "$REF_ID" = "phix" ]; then
    echo This is the PhiX QC option, this will align the undetermined reads to the PhiX genome, \
    if they align well the sequencer is likely working correctly.

    if [ ! -d $EXP_DIR/PhiX ]; then
    echo "$EXP_DIR/PhiX not found checking for $EXP_DIR/Undetermined"
    elif [ ! -d $EXP_DIR/Undetermined ]; then
        echo "PhiX folder not found in experimental data. Files found:"
        return
    fi
    
    cp $HOMEPATH/ngs_pipeline/references/PhiX* $TMP_DIR
    
    fastp --thread 16 --trim_poly_x --in1 $EXP_DIR/Undetermined*R1_001.fastq.gz --in2 $EXP_DIR/Undetermined*R2_001.fastq.gz \
    --out1 $TMP_DIR/R1.fsp.fq.gz --out2 $TMP_DIR/R2.fsp.fq.gz || pipeline_error
    bwa mem $TMP_DIR/PhiX.fasta $TMP_DIR/*fsp.fq.gz > $TMP_DIR/${EVX_RUNID}_PhiX.sam || pipeline_error
    samtools view -bSq 2 $TMP_DIR/${EVX_RUNID}_PhiX.sam > $TMP_DIR/${EVX_RUNID}_x1.bam || pipeline_error
    samtools view -q 30 -b $TMP_DIR/${EVX_RUNID}_x1.bam > $TMP_DIR/${EVX_RUNID}_x2.bam || pipeline_error
    samtools view -F 0x04 -b $TMP_DIR/${EVX_RUNID}_x2.bam > $TMP_DIR/${EVX_RUNID}_PhiX.bam || pipeline_error
    samtools sort $TMP_DIR/${EVX_RUNID}_PhiX.bam > $TMP_DIR/${EVX_RUNID}_sort_PhiX.bam || pipeline_error
    samtools index $TMP_DIR/${EVX_RUNID}_sort_PhiX.bam || pipeline_error

    mv $TMP_DIR/* $OUT_DIR
    rm -rf $TMP_DIR

    echo PhiX aligment is complete, check the results in $OUT_DIR, open the sort.PhiX.bam file with IGV to see if it aligned correctly.
    return

fi
#####################################################################################################
# NGS pipeline
#####################################################################################################

# Function to capture pipeline errors, INFO captured when function called, input fastq returned to inputs upon error
PIPELINE_ERROR="false"
function pipeline_error {
    PIPELINE_ERROR="true"
    DEBUG="true"
    INFO=$1
    cp $EXP_DIR/*R* $HOMEPATH/ngs_pipeline/inputs
}

# copying the refasta_file into the tmp file

echo $REF_FILE
echo $TMP_DIR/$REFASTA
echo $REFASTA

sleep 3



cp $REF_FILE $TMP_DIR/$REFASTA

echo "**************************************************************"
echo "Starting NGS Pipeine run for $EVX_RUNID with reference $REF_ID"
echo "**************************************************************"

############################################
#QC transferred to autopipeline.bash for multi-site
###########################################


# run fastqc
#echo "**************************************************************"
#echo "Starting fastq qc checks with fastqc for $EVX_RUNID with $REF_ID"
#echo "**************************************************************"
#fastqc $EXP_DIR/*.fastq.gz --outdir $TMP_DIR || pipeline_error "Fastqc Error"

# FASTQ preprocessing using fastp, saves preprocessed files
#echo "**************************************************************"
#echo "Starting pre-processing of fastq files with fastp for $EVX_RUNID with $REF_ID"
#echo "**************************************************************"
#fastp --thread 16 --trim_poly_x --in1 $EXP_DIR/*R1*.fastq.gz \
#--in2 $EXP_DIR/*R2*.fastq.gz --out1 $TMP_DIR/R1.fsp.fq.gz \
#--out2 $TMP_DIR/R2.fsp.fq.gz --json $TMP_DIR/${EVX_RUNID}_fastp.json \
#--html $TMP_DIR/${EVX_RUNID}fsp_data.html \
#--report_title ${EVX_RUNID}_data_fsp 1 > $TMP_DIR/${EVX_RUNID}_fsp.log \
#|| pipeline_error "Fastq Error"

# align to reference: burrows wheeler aligner, written to handle low-divergent sequences
# index first then process.

####################################################################
#Copy trimmed fastq files from QC directory
####################################################################

cp $QC_DIR/*fsp.fq.gz $TMP_DIR

####################################################################
#Begin the main pipeline, multiple pipeline runs against 1 pair of fastq files
####################################################################

echo "**************************************************************"
echo "Aligning the processed fastq files for $EVX_RUNID to the reference $REF_ID"
echo "**************************************************************"
bwa index $TMP_DIR/$REFASTA || pipeline_error  "Reference indexing Error"
# bwa mem is the recommended aligner out of 3 possibles
bwa mem -k $SEED $TMP_DIR/$REFASTA $TMP_DIR/*fsp.fq.gz > $TMP_DIR/$EVX_RUNID.sam || pipeline_error  "Alignment Error"

# samtools: 
echo "**************************************************************"
echo "Filtering alignment files with SAMTOOLS for $EVX_RUNID with $REF_ID"
echo "**************************************************************"
#faidx indexes reference fasta to .fasta.fai
samtools faidx $TMP_DIR/$REFASTA || pipeline_error "Samtools Error"
# convert .sam to .bam
#Filter SAM file to only contain aligned reads
samtools view -h $TMP_DIR/$EVX_RUNID.sam -F 0x04 > $TMP_DIR/${EVX_RUNID}_filtered.sam
samtools view -bSq 2 $TMP_DIR/${EVX_RUNID}_filtered.sam > $TMP_DIR/${EVX_RUNID}_x1.bam || pipeline_error "Samtools Bam convert Error"
# samtools: read Q30 reads only (above q30 quality) and pipe to a new .bam - essentially filtering on quality
samtools view -q 2 -b $TMP_DIR/${EVX_RUNID}_x1.bam >  $TMP_DIR/${EVX_RUNID}.bam || pipeline_error "Samtools Filter Error"
# samtools: remove unmapped reads -F is used to only exclude certain reads, 0x04 means mapped
#samtools view -F 0x04 -b $TMP_DIR/${EVX_RUNID}_x2.bam > $TMP_DIR/$EVX_RUNID.bam || pipeline_error "Samtools Filter Error"
# samtools: sorts the bam file and creates a new sorted .bam
samtools sort $TMP_DIR/$EVX_RUNID.bam > $TMP_DIR/${EVX_RUNID}_sort.bam  || pipeline_error "Samtools Soring Error"
# samtools: index sorted bam file
samtools index $TMP_DIR/${EVX_RUNID}_sort.bam  || pipeline_error "Samtools Indexing Error"


# Summary stats for sequencing file and alignment, if there are insufficient reads the pipeline will exit here#

# Carry out counts - the counting of sam alignment takes a while so only do it once
echo "**************************************************************"
echo "Checking aligment and recording the stats for $EVX_RUNID with $REF_ID"
echo "**************************************************************"
R1count=$(( $(gunzip -c $EXP_DIR/*R1_*.fastq.gz|wc -l)/4|bc ))
R2count=$(( $(gunzip -c $EXP_DIR/*R2_*.fastq.gz|wc -l)/4|bc ))

prefilter=$(samtools flagstat $TMP_DIR/$EVX_RUNID.sam)
prefilterlist=$(awk '{print $1;}' <<< "$prefilter")
prefilternumber=$(head -n 1 <<< $prefilterlist)
postfilter=$(samtools flagstat $TMP_DIR/${EVX_RUNID}_sort.bam)
postfilterlist=$(awk '{print $1;}' <<< "$postfilter")
postfilternumber=$(head -n 1 <<< $postfilterlist)

#initiate summary stats file and add read counts

STATS=$TMP_DIR/${EVX_RUNID}_stats.txt
echo Sample ID: ${EVX_RUNID} >> ${STATS}
echo $EVX_RUNID >> $STATS
echo Fastq files last modified: >> $STATS
stat $EXP_DIR/*R1*fastq.gz | grep Modify >> $STATS
echo Reference: >>$STATS
echo $REFASTA >> $STATS
sed -n '2p' < $TMP_DIR/$REFASTA >> $STATS
echo >> $STATS
echo Synthetic DNA: >> $STATS
echo $SYNDNA >> $STATS
echo Reads from forward fastq file: >> $STATS
echo $R1count >> $STATS
echo Reads from reverse fastq file: >> $STATS
echo $R2count >> $STATS
echo Total reads: >> $STATS
echo $prefilter >> $STATS
echo Reads that aligned: >> $STATS
echo $postfilter >> $STATS
echo $refID,$postfilternumber >> $POOLSTATS/pooled_stats_read_numbers.csv

#initiate csv summary stats
STATSCSV=$TMP_DIR/${EVX_RUNID}_stats.csv
echo \#Run ID,${EVX_RUNID} >> $STATSCSV
echo Seq ID,${refID} >> $STATSCSV
echo Fastq-R1 , $R1count >> $STATSCSV
echo Fastq-R2 , $R2count >> $STATSCSV
echo Total-reads , $prefilternumber >> $STATSCSV
echo Total-aligned , $postfilternumber >> $STATSCSV

#If there are insufficient aligned reads, currently set to 100, just carry out the unaligned analysis and then exit copying files to ERROR_Analysis
if (( $postfilternumber < 100 )); then
	CRASHFILE=$TMP_DIR/ERROR_${EVX_RUNID}.txt

    echo "********* There are insufficent aligned reads in file: " $EVX_RUNID " to continue the aligned analysis *********"
    sleep 2
    echo "********* The sequencing files will now be scanned to see what has been seqeunced **********"
    sleep 2
    
    echo Pipeline exited due to insufficient aligned reads - Currently threshold is less than 100 >> $CRASHFILE
    echo Number of aligned reads: >> $CRASHFILE
    echo $postfilternumber >> $CRASHFILE
    echo "$EVX_RUNID,$REF_ID,FAIL,Insufficiant aligned reads" >> $LOGFILE
    echo "********* Input files will be returned to inputs folder *********"
    cp $EXP_DIR/*R* $HOMEPATH/ngs_pipeline/inputs

    #Run the whats in it? script to see what has been sequenced. Turned off for multi-site as it runs at the start
    #python unaligned_reads.py $TMP_DIR/ $SYNDNA $EVX_RUNID
    
    #Carry out a multiple alignment on output from whats in it
    #muscle -in $TMP_DIR/*_top15_reads.fasta -out $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp

    #Re-order the aligment so the most abundant variants are at the top (Using amended python script from Muscle author)
    #python stable.py $TMP_DIR/*_top15_reads.fasta $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp > $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta

    #Remove tempory alignment files
    #rm $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp

    #Make the multiple alignment plots using pymsaviz (https://pypi.org/project/pymsaviz/)
    #python unaligned_reads_plot.py $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta $EVX_RUNID $TMP_DIR

    #Create ERROR_Analysis sub folders named ERROR_Analysis_1 ERROR_Analysis_2 etc.....
    declare -i SYN_NUMBER=1
    SYN_DIR="$OUT_DIR/ERROR_Analysis_$SYN_NUMBER"

    #The until loop runs until the folder number is no longer present, this enables the same parent folder to be used for subsequent analysis
    #If it is present the number is increased until it is no longer present
    until [ ! -d $SYN_DIR ]; do
        if [ -d $SYN_DIR ]; then
            let SYN_NUMBER=SYN_NUMBER+1
            SYN_DIR="$OUT_DIR/ERROR_Analysis_$SYN_NUMBER"
        fi
    done

    #The folder to put the output into is made outside the until loop
    mkdir $SYN_DIR

    #Everything is now moved into this folder
    #mv $TMP_DIR/*_fastqc* $SYN_DIR
    mv $TMP_DIR/*.csv $SYN_DIR
    #mv $TMP_DIR/*.html $SYN_DIR
    #mv $TMP_DIR/*.png $SYN_DIR
    mv $TMP_DIR/*.sam $SYN_DIR
    mv $TMP_DIR/*.fasta $SYN_DIR
    mv $STATS $SYN_DIR
    mv $CRASHFILE $SYN_DIR
    rm -r $TMP_DIR

    echo "********* The pipeline run for $EVX_RUNID with $REF_ID will exit now, check the error files in the ERROR directory *********"

    conda deactivate
    echo "*********If running a multi analysis the pipeline will resume shortly**********"
    sleep 5
    return 5
fi	

#The main pipeline continues if the aligned reads are above 100.

#Calculate the auto threshold for the variant distribution script, currently set to calculate 1% of reads, previously was 
#calculated based upon the sequecing quality, this could be re-activated if required.
#Applying sequecing quality to threshold was too strict, sequencing errors are totally random accross all variants but sythesis errors occur
#in the same place in multiple variants due to PCR amplification.

echo "User threshold used for Variant Distribution:" >> $STATS

if [[ $THRESH = "AUTO" ]]; then
    echo "Pipeline in auto threshold mode"
    echo "Calculating Threshold for 1% of reads"
    THRESH=$(bc<<<"$postfilternumber*0.01")
    echo "1 % Threshold = $THRESH"
    echo "Threshold calculated for 1% of reads with an aligned read number of $postfilternumber" >> $STATS

    #Check not reauired below as sequencing quality not currently used
    #if [[ -z "$QUAL" ]]; then
    #    echo "Sequrncing quality not entered"
    #    echo "Please enter --quality to enable threshold calcuation"
    #    echo "Defaulting to 400 copy threshold"
    #    THRESH=400
    #fi
    sleep 2
fi

echo $THRESH >> $STATS

# Java
# jvarkit is a java package for bioinformatics, so java must be installed to run this script
# bioalcidaejdk is a java based version of awk, this pulls out the matches with some criteria
# and gives .csv with the matches.
echo "**************************************************************"
echo "Checking the number of matches at each point in reference with jvarkit for $EVX_RUNID with $REF_ID"
echo "**************************************************************"
java -jar ~/jvarkit/dist/jvarkit.jar bioalcidaejdk -e \
'stream().filter(R->!R.getReadUnmappedFlag() && R.getCigar()!=null).forEach(R->println(R.getReadName()+
"\t"+R.getContig()+"\t"+R.getStart()
+"\t"+R.getCigar().getCigarElements().stream().filter(CE->CE.getOperator().isAlignment()).mapToInt(CE->CE.getLength()).sum()));' \
$TMP_DIR/${EVX_RUNID}_sort.bam | awk '{print $4}' >> $TMP_DIR/${EVX_RUNID}_matches.csv 

# Awk: fetch reference sequence
SEQ=$(awk '{print $1}' $TMP_DIR/${REFASTA}*.fai) || pipeline_error "Error extracting Sam Index" #This must use awk to grab the reference sequences
# pysamstats: calculates summary of each read from alignment file in terms of mismatches, deletions, insertions etc
echo "**************************************************************"
echo "Carrying out variant calling with Pysamstats for $EVX_RUNID with $REF_ID"
echo "**************************************************************"
pysamstats --type variation --max-depth=1000000 --chromosome ${SEQ} \
-f $TMP_DIR/$REFASTA $TMP_DIR/${EVX_RUNID}_sort.bam > $TMP_DIR/${EVX_RUNID}_var_${SEQ}.csv || pipeline_error "Error with Pysamstats"

# concatenate json files into one .csv
#cat $TMP_DIR/*json | head -41 | tail -1 >> $TMP_DIR/${EVX_RUNID}_insert_size.csv || pipeline_error "Error concatenating json"

# using the python script
echo "*******************************************"
echo "Creating the plots for $EVX_RUNID with $REF_ID"
echo "*******************************************"
python ngs_python_wrapper.py $TMP_DIR/${EVX_RUNID}_var_${SEQ}.csv $TMP_DIR/${EVX_RUNID}_matches.csv $QC_DIR/${SAMPLE}_insert_size.csv \
$TMP_DIR/${EVX_RUNID}_filtered.sam $SYNDNA $TMP_DIR $EVX_RUNID $TMP_DIR/$REFASTA ${THRESH%%.*} || pipeline_error "Error with python plotting scripts"

###################################################################################################
# saving outputs
#Outputs are saved in folders named SYN_1, SYN_2, etc.
#This allows for multiple analysis against the same reference for the same pair of fastq files
#Removal of previous analysis is now no longer required
###################################################################################################

declare -i SYN_NUMBER=1
SYN_DIR="$OUT_DIR/Analysis_$SYN_NUMBER"

#The until loop runs until the folder number is no longer present
#If it is present the number is increased until it is no longer present
until [ ! -d $SYN_DIR ]; do
    if [ -d $SYN_DIR ]; then
        let SYN_NUMBER=SYN_NUMBER+1
        SYN_DIR="$OUT_DIR/Analysis_$SYN_NUMBER"
    fi
done

#adding pass/fails flags file
python flagger.py --syndna $SYNDNA --name $EVX_RUNID --input_dir $TMP_DIR

cp $TMP_DIR/${EVX_RUNID}_stats.csv $HOMEPATH/ngs_pipeline/outputs/$MULTIRUNID/Pooled-summary/${EVX_RUNID}_${SYNDNA}_stats.csv

#The folder to put the output into is made outside the until loop
mkdir $SYN_DIR

#Delete large unfiltered sam file
rm $TMP_DIR/$EVX_RUNID.sam

#Everything is now moved into this folder
#mv $TMP_DIR/*_fastqc* $SYN_DIR
mv $TMP_DIR/*.csv $SYN_DIR
mv $TMP_DIR/*.html $SYN_DIR
mv $TMP_DIR/*.png $SYN_DIR
mv $TMP_DIR/*.sam $SYN_DIR
mv $TMP_DIR/flags.txt $SYN_DIR
mv $STATS $SYN_DIR


##################################################################################################
# cleanup
#Errors are written to log file if present
##################################################################################################

if [ $PIPELINE_ERROR = "true" ]; then
    echo "Error running pipeline for $EVX_RUNID with $REF_ID, problem reported as $INFO, check $TMP_DIR for outputs"
    echo "$EVX_RUNID,$refID,FAIL,$INFO" >> $LOGFILE 
    echo "If running a multi analysis the pipeline will resume shortly"
    sleep 5
    else 
        echo "$EVX_RUNID,$refID,PASS" >> $LOGFILE 
        echo "********** Main NGS pipeline aligned analysis complete for $EVX_RUNID with $REF_ID **********"
        sleep 2
fi

##################################################################################################
#Initiate the whatsinit script to get the top 15 most abundant variant from unfiltered sequencing file
#Turned off for multi-site as it runs at the start
##################################################################################################

#echo "*******************************************"
#echo "Carrying out unaligned analysis"
#echo "*******************************************"

#python unaligned_reads.py $TMP_DIR/ $SYNDNA $EVX_RUNID

#echo "*******************************************"
#echo "Carrying out multiple alignment with Muscle for unaligned reads"
#echo "*******************************************"
#Multiple aligner Muscle https://www.drive5.com/muscle/
#muscle -in $TMP_DIR/*_top15_reads.fasta -out $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp

#Re-order the aligment so the most abundant variants are at the top (Using amended python script from Muscle author)
#python stable.py $TMP_DIR/*_top15_reads.fasta $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp > $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta

#Remove tempory alignment files
#rm $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta.tmp

#Make the multiple alignment plots using pymsaviz (https://pypi.org/project/pymsaviz/)
#python unaligned_reads_plot.py $TMP_DIR/${EVX_RUNID}_top15_reads.aligned.fasta $EVX_RUNID $TMP_DIR

#mkdir $SYN_DIR/Unaligned_Analysis
#mv $TMP_DIR/*aligned.fasta  $SYN_DIR/Unaligned_Analysis
#mv $TMP_DIR/*png $SYN_DIR/Unaligned_Analysis
#mv $TMP_DIR/*csv $SYN_DIR/Unaligned_Analysis

if [ ! $DEBUG == "true" ]; then
    rm -r $TMP_DIR
else
    echo "*******************************************"
    echo "Pipeline is in debugging mode, checking unaligned reads, and retaining all files in debug directory"
    echo "*******************************************"
    sleep 2
    mv $TMP_DIR $SYN_DIR
fi

conda deactivate

################################################################################################## 
# Due to running as source in case the same terminal is used twice
# Unset MULTIRUNID if pipeline is being run in single mode
# In multi mode this is handled by the multifastq bash script instead
################################################################################################## 

if [[ $MODE == "single" ]]; then
    unset MULTIRUNID
    unset THRESHMODE
    unset QUAL
    unset THRESHMODE
fi
