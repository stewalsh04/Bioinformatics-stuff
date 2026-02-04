#!/bin/sh -x
# handles wrong argument errors
function exit_with_bad_args {
    echo "Usage: . ./NGSPipeline-Install-WSL.sh [--update]"
    echo "Invalid arguments provided" >&2
    return # this stops the terminal closing when run as source
}

options=$(getopt -o '' -l update -- "$@") || exit_with_bad_args

# default arguments
UPDATE="false"

# handling the commandline options
eval set -- "$options"
while true; do
    case "$1" in
        --update)
            UPDATE="true"
        echo "Updating..."
            ;;
        --)
            shift
            break
            ;;
    esac
    shift
done

#Check if installation script is in update mode, if so attempt to update
if [[ $UPDATE == "true" ]]; then
    echo "Attempting to update pipeline"
    echo "Checking for current installation"

    #Check if conda is installed
    if command -v conda; then
        echo "Conda detected"
        sleep 2
    else 
        apt-get update
        apt-get install curl -y
        curl -O https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
        bash Anaconda3-2020.02-Linux-x86_64.sh
        source ~/.bashrc
        conda env create -f environment.yml
        conda activate NGS-Pipeline
    fi

    #Check if conda environment is installed
    if conda env list | grep "NGS-Pipeline"; then
        echo "Looks like the conda environment is installed"
        echo "updating conda environment"
        # update environment from environment.yml, uninstallign any no longer required dependencies
        conda env update --name NGS-Pipeline --file environment.yml --prune
        conda activate NGS-Pipeline
        sleep 2
    else
       conda env create -f environment.yml
       conda activate NGS-Pipeline
    fi
     
    #Check if bc is installed (For threshold calculation)
    if command -v bc &> /dev/null; then
        echo "bc already installed"   
    else
        sudo apt install bc
    fi

    #Check if required folders are available
    echo "Checking for required folders"
    
    if [ ! -d /mnt/c/ngs_fastq/logs ]; then
        echo "Directory: ngs_fastq not found, creating now"
        mkdir -p /mnt/c/ngs_fastq/logs
        cp ref.txt /mnt/c/ngs_fastq
    else 
        echo "ngs_fastq folder found"
    fi

    if [ ! -d /mnt/c/ngs_repository ]; then
        echo "Directory: ngs_repository not found, creating now"
        mkdir /mnt/c/ngs_repository
    else 
        echo "ngs_repository folder found"
    fi

    if [ ! -d /mnt/c/ngs_experiment_data ]; then
        echo "Directory: ngs_experiment_data not found, creating now"
        mkdir /mnt/c/ngs_experiment_data
    else 
        echo "ngs_experiment_data folder found"
    fi

    # copying files
    cp ref.txt /mnt/c/ngs_fastq/
    cp -r PhiX /mnt/c/ngs_repository

    # reinstalling jvarkit
    git submodule init
    git submodule update
    cd jvarkit
    ./gradlew jvarkit
    cd ..
    cp -r jvarkit ~/

    if [ -d ~/jvarkit ] && [ -f ~/jvarkit/gradlew ] && [ -f ~/jvarkit/dist/bioalcidaejdk.jar ] && [ -f ~/jvarkit/dist/jvarkit.jar ]; then
        echo "Bioalcidaejdk Found, Jvarkit found"
    else
        echo "oh no"
    fi

else
    if [ -d /mnt/c/ngs_experiment_data ] || [ -d /mnt/c/ngs_repository ] || [ -d /mnt/c/ngs_fastq ] || [ -d ~/jvarkit ]; then
        echo "Previous partial or full installation found, please rerun command with --update flag"

    else
        apt-get update
        apt-get install curl -y
        curl -O https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh 
        bash Anaconda3-2020.02-Linux-x86_64.sh
        rm Anaconda3-2020.02-Linux-x86_64.sh
        source ~/.bashrc
        conda config --set auto_activate_base false
        conda env create -f environment.yml
        conda activate NGS-Pipeline
        mkdir /mnt/c/ngs_repository
        mkdir /mnt/c/ngs_experiment_data
        mkdir -p /mnt/c/ngs_fastq/logs
        cp ref.txt /mnt/c/ngs_fastq/
        # cd jvarkit
        # ./gradlew jvarkit
        # cd ..
        # cp -r jvarkit ../
        cp -r PhiX /mnt/c/ngs_repository
        if ! command -v bc &> /dev/null; then
             sudo apt install bc
        fi
        echo "To complete installation restart terminal, navigate back to this directory and run . ./NGSPipeline.bash --update"
    fi
fi

#Install the multiple alignment plotter, this cannot be installed through conda! No idea why!
pip install pymsaviz

conda deactivate
