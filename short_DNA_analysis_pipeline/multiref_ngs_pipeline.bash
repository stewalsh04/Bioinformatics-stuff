function exit_with_bad_args {
    echo "Usage: source ngs_pipeline.bash --id <EV_ID> --ref <REF_FILE>"
    echo "Invalid arguments provided" >&1
    return # this stops the terminal closing when run as source
}

options=$(getopt -o '' -l id: -l ref: -- "$@") || exit_with_bad_args

# handling the commandline options
eval set -- "$options"
while true; do
    case "$1" in
        --id)
            shift
            EVX_RUNID="$1"
            ;;
        --ref)
            shift
            REF_ID="$1"
            ;;
        --)
            shift
            break
            ;;
    esac
    shift
done

echo $EVX_RUNID
echo $REF_ID

# check for required args
if [ -z "$EVX_RUNID" ] || [ -z "$REF_ID" ]; then
    echo "--id and --ref  require arguments"
    exit_with_bad_args
fi


ID=$EVX_RUNID
REFERENCEFILE=$REF_ID

REFERENCES=($(awk '{print $1}' $REFERENCEFILE))
SYNSEQS=($(awk '{print $2}' $REFERENCEFILE))

for ((i=0; i<${#REFERENCES[@]}; i++)); 
    do 
        source ngs_pipeline.bash --id $ID --ref "${REFERENCES[$i]}" --seq "${SYNSEQS[$i]}"; 
    done
