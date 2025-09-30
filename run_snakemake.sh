#!/bin/bash

helpFunction()
{
        echo "Usage: bash run_snakemake.sh [-w /path/to/working_directory]"
        echo -e "\t-w <path> Absolute Path for Working Directory."
        echo -e "\t-t <int> Maximum Number of Threads Allowed For the Pipeline. [Default: 48]"
        exit 1
}

while getopts "w:t:" opt
do
        case "$opt" in
        w )
                workdir="$OPTARG"
                ;;
        t )
                threads="$OPTARG"
                ;;
        ? )     helpFunction ;;
        esac
done

threads=48

if [ -z "$workdir" ]
        then
        echo "Please Provide the Path to the Working Directory"
        helpFunction
fi


path=$(conda info --base)

cwd=$(dirname "$(readlink -f "$0")")

workdir_full=$(realpath "$workdir")

source $path/bin/activate anomaly

cd "$workdir_full" || exit

threads_sniffles=$(grep "threads_sniffles" snake_config.yml | awk -F " " '{print $2}')

threads_minimap2=$(grep "threads_minimap2" snake_config.yml | awk -F " " '{print $2}')

max=$(( $threads_sniffles >= $threads_minimap2 ? $threads_sniffles : $threads_minimap2 ))

jobs=$(($threads/$max))

snakemake --conda-frontend conda --use-conda --rerun-incomplete --cores $threads --jobs $jobs

cd "$cwd" || exit