#!/bin/bash

# Collect all the csv files for each type of fit in a -t bin into one directory
# Expects the following directory structure created by neutralb1.batch.submit
#   -t bin
#     - mass_X.YY bins
#       - primary fit results (best.csv, data.csv, etc.)
#           - bootstrap fit results
#           - randomized fit results

usage() {
    echo "Usage: collect_all_fits.bash -d <input_directory> -o <output_directory> [-p] [-h]"
    echo "  -d: -t bin directory containing subdirectories with fit results"
    echo "  -o: Output directory to store collected csv files"
    echo "  -p: Preview mode (do not execute commands)"
    echo "  -h: Display this help message"    
}

preview_mode="false"

while getopts ":d:o:ph" opt; do
    case $opt in
        d) input_dir="$OPTARG"
        ;;
        o) output_dir="$OPTARG"
        ;;
        p) preview_mode="true"
        ;;
        h) 
            usage
            exit 0
        ;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
        ;;
        :) echo "Option -$OPTARG requires an argument." >&2
            exit 1
        ;;
    esac
done

if [ -z "$input_dir" ]; then
    echo "Input directory is required."
    usage
    exit 1
fi

if [ -z "$output_dir" ]; then
    output_dir="./"
fi

mass_dirs=($(find "$input_dir" -maxdepth 1 -type d -name "mass_*"))
if [ ${#mass_dirs[@]} -eq 0 ]; then
    echo "No directories with 'mass_' found in $input_dir."
    exit 1
fi

# use preview flag to printout what files will be combined if in preview mode
if [ "$preview_mode" == "true" ]; then
    echo "Preview mode enabled. The following files will be combined:"
    preview="-p"
else
    preview=""
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/data.csv -o ${output_dir}/data.csv $preview
python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best.csv -o ${output_dir}/best.csv $preview
python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best_acceptance_corrected.csv -o ${output_dir}/best_acceptance_corrected.csv $preview

# handle optional files
files=(${input_dir}/mass*/best_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best_projected_moments.csv -o ${output_dir}/moments.csv $preview
fi

files=(${input_dir}/mass*/rand/rand.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand.csv -o ${output_dir}/rand.csv $preview
fi

files=(${input_dir}/mass*/rand/rand_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand_acceptance_corrected.csv -o ${output_dir}/rand_acceptance_corrected.csv $preview
fi

files=(${input_dir}/mass*/rand/rand_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand_projected_moments.csv -o ${output_dir}/rand_moments.csv $preview
fi

files=(${input_dir}/mass*/bootstrap/bootstrap.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap.csv -o ${output_dir}/bootstrap.csv $preview
fi

files=(${input_dir}/mass*/bootstrap/bootstrap_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap_acceptance_corrected.csv -o ${output_dir}/bootstrap_acceptance_corrected.csv $preview
fi

files=(${input_dir}/mass*/bootstrap/bootstrap_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap_projected_moments.csv -o ${output_dir}/bootstrap_moments.csv $preview
fi

files=(${input_dir}/mass*/truth/best.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best.csv -o ${output_dir}/truth.csv $preview
fi

files=(${input_dir}/mass*/truth/best_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best_acceptance_corrected.csv -o ${output_dir}/truth_acceptance_corrected.csv $preview
fi

files=(${input_dir}/mass*/truth/best_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best_projected_moments.csv -o ${output_dir}/truth_projected_moments.csv $preview
fi