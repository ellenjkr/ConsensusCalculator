#!/bin/bash
# bash virtual_pcr.sh -f CCCAGTACAGCCCATTCACC -r AAAGCGGATAGGGCTCCTGTA -i FAT05083_pass_a4a97209_0.fastq.gz -m 7500 -l 9000


# ============================================================================================================
# The first part of the program extracts sequences between primers. 
# The implementation was extracted from this project in github https://github.com/Nucleomics-VIB/InSilico_PCR

version="1.1.3; 2022-02-08"

source /home/bioinfo/miniconda3/etc/profile.d/conda.sh
conda activate InSilico_PCR || \
  ( echo "# the conda environment 'InSilico_PCR' was not found on this machine" ;
    echo "# please read the top part of the script!" \
    && exit 1 )

########## Variables Initialization ##############

# REM: the workflow can be restarted by deleting all log files after a given step

# speed-up
thr=8
pigt=2
mem="1g"

# extract reads corresponding to the 16S PCR
forwardl="Forward"
reversel="Reverse"

# be stringent to avoid noisy reads
cut=0.85

barcodes_len=5

# split in 500k read bins and zip
lines=2000000

# quality phred scale of your data (for demo data it is 33)
qual=33

# expected amplicon size limits, set to 10 and 10000 by default
# adjust if you notice that the sequence extraction has unwanted tail(s)
readminlen=10
readmaxlen=100000

######## end of user editable region ############

# Arguments



while getopts ho:t:m:M:n:f:r:i: flag
do
    case "${flag}" in
        h) # display Help
            Help
            exit 1;;
        o) outfile=${OPTARG};;
        t) thr=${OPTARG};;
        m) readminlen=${OPTARG};;
        M) readmaxlen=${OPTARG};;
        n) primername=${OPTARG};;
        f) forwardp=${OPTARG};;
        r) reversep=${OPTARG};;
        i) infile=${OPTARG};;
    esac
done

shift $(( OPTIND - 1 ))

if [ -z "$primername" ]; then
        echo 'Missing -n (Primer Pair Name)'
        exit 1
elif [ -z "$forwardp" ]; then
        echo 'Missing -f (Forward Primer)'
        exit 1
elif [ -z "$reversep" ]; then
        echo 'Missing -r (Reverse Primer)'
        exit 1
elif [ -z "$infile" ]; then
        echo 'Missing -i (Sequences input file)'
        exit 1      

fi

# extract string for output names
name=$(basename ${infile%\.fasta})

##################
# prepare folders

# working default to local folder
data="$(pwd)"

# run logs
logs="run_logs_${name}"
mkdir -p ${logs}

# tmp output folder
tmpout="bbmap_out"
mkdir -p ${tmpout}

# parallel folders
WORKDIR=$PWD
mkdir -p $PWD/tmp
export TMPDIR=$PWD/tmp

# keep track of all
runlog=${logs}/runlog.txt
exec &> >(tee -i ${runlog})

######################################################
# find forward primer

echo "# searching for forward primer sequence: ${forwardp} in all files"
msa.sh -Xmx${mem} threads=10 qin=${qual} in=${infile} out=${tmpout}/forward.sam literal="${forwardp}" rcomp=t addr=t replicate=t cutoff="${cut}"

######################################################
# find reverse primer

echo "# searching for reverse primer sequence: ${reversep} in all files"
msa.sh -Xmx${mem} threads=10 qin=${qual} in=${infile} out=${tmpout}/reverse.sam literal="${reversep}" rcomp=t addr=t replicate=t cutoff="${cut}"

##########################################
# extract regions with BBMap cutprimers.sh

cutprimers.sh -Xmx${mem} qin=${qual} in=${infile} out=${tmpout}/${name}.fasta sam1=${tmpout}/forward.sam sam2=${tmpout}/reverse.sam include=f fixjunk=t

cp -f ${tmpout}/${name}.fasta ${infile}

conda deactivate

# ===============================================================================

rm -r $logs
rm -r $tmpout
rm -r $runlog
rm -r $PWD/tmp

