#!/bin/bash

helpFunction()
{
   echo -e "Usage: get_config.sh [-d b] [-r reference.fasta] [-m 16] [-s 16] [-p ONT] [-i input directory] [-o output directory] [-l reference headers list]\n[-q Minimum Mapping Quality] [-n Minimum Read Support for SV] [-g Genotype Ploidy]"
   echo -e "\t-d <string> Specify if Input is BAM ('b') or FASTQ ('f'). [default: b]"
   echo -e "\t-r <filename> Name of the Reference Nuclear Genome Fasta File (Include the Path if the fasta file is not in the same directory)"
   echo -e "\t-m <int> Number of threads to be used for Minimap2. [default: 16]"
   echo -e "\t-s <int> Number of threads for Sniffles. [default: 16]"
   echo -e "\t-p <string> Sequencing Platform (ONT/pb/HiFi). [default: ONT]"
   echo -e "\t-i <path> Absolute Path containing the FASTQ or BAM Files."
   echo -e "\t-o <path> Absolute Path for Pipeline Output and Intermediate Files."
   echo -e "\t-l <path> Absolute Path for the list containing chromosome reference headers."
   echo -e "\t-q <int> Minimum Map Quality for calling SVs. [default: 0]"
   echo -e "\t-n <int> Minimum Support for calling SVs. [default: 4]"
   echo -e "\t-g <int> Genotype Ploidy. [default: 2]"
   echo -e "\t-e <int> Blast E-value Cutoff. [default: 1e-3]"
   echo -e "\t-c <int> NUMT Coverage Cutoff. [default: 70]"
   echo -e "\t-a <int> Minimum Support for calling NUMTs from Supplementary Alignments. [default: 5]"
   exit 1 # Exit script after printing help
}

# Default values

data=b

threads_m=16

threads_s=16

platform=ONT

minq=1

mins=4

ploidy=2

evalue_cutoff=1e-3

coverage=70

support=5



while getopts "d:r:m:s:p:i:o:l:q:n:g:e:c:a:" opt
do
	case "$opt" in
	d )
		data="$OPTARG"
		;;
	r )
		ref="$OPTARG"
		;;
	m )
		threads_m="$OPTARG"
		;;
	s )
		threads_s="$OPTARG"
		;;
	p )
		platform="$OPTARG"
		;;
	i )
		sampledir="$OPTARG"
		;;
	o )
		workdir="$OPTARG"
		;;
	l )
		ref_list="$OPTARG"
		;;
	q )
		minq="$OPTARG"
		;;
	n )
		mins="$OPTARG"
		;;
	g )
		ploidy="$OPTARG"
		;;
	e )
		evalue_cutoff="$OPTARG"
		;;
	c )
		coverage="$OPTARG"
		;;
	a )
		support="$OPTARG"
		;;
	? )	helpFunction ;;
	esac
done

if [ -z "$sampledir" ] || [ -z "$workdir" ] || [ -z "$ref" ] || [ -z "$ref_list" ]
	then
	echo "Please Provide the Reference Fasta, Sample Directory and Working Directory"
	helpFunction
fi

workdir_path=$(realpath "$workdir")

sampledir_path=$(realpath "$sampledir")

ref_full=$(realpath "$ref")

ref_list_full=$(realpath "$ref_list")

if [ ! -d "$workdir" ]; then
	mkdir "$workdir"
fi

cat > $workdir_path/snake_config.yml <<EOF
# Config file for the pipeline
# General configuration
bam_or_fastq: $data  # Specify if the input is BAM ('b') or FASTQ ('f')
reference_nuclear: $ref_full  # Path to the nuclear reference genome
headers_nuclear: $ref_list_full
threads_minimap2: $threads_m  # Number of threads for Minimap2
read_type: $platform  # Sequencing platform (e.g., for Oxford Nanopore reads mention ONT; for PacBio reads mention 'pb'; for PacBio HiFi reads mention 'HiFi')
sample_directory: $sampledir_path  # Directory containing sample files
working_directory: $workdir_path  # Directory for pipeline output and intermediate files

# Sniffles parameters
minimum_numt_length: 15
min_support: $mins
minimum_mapq: $minq
long_ins_length: 2500
genotype_ploidy: $ploidy
quiet: True
allow_overwrite: False
threads_sniffles: $threads_s

# Other parameters
evalue_cutoff: $evalue_cutoff
numt_coverage: $coverage
numt_supporting_reads: $support
EOF
