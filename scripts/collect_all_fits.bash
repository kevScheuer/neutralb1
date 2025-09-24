#!/bin/bash

# Collect all the csv files for each type of fit in a -t bin into one directory
# Expects the following directory structure created by neutralb1.batch.submit
#   -t bin
#     - mass_X.YY bins
#       - primary fit results (best.csv, data.csv, etc.)
#           - bootstrap fit results
#           - randomized fit results

usage() {
    echo "Usage: collect_all_fits.bash -i <input_directory> -o <output_directory> [-p] [-h]"
    echo "  -i: -t bin directory containing subdirectories with fit results"
    echo "  -o: Output directory to store collected csv files"
    echo "  -p: Preview mode (do not execute commands)"
    echo "  -f: Force overwrite of existing output files"
    echo "  -h: Display this help message"
}

preview_mode="false"

while getopts ":i:o:pfh" opt; do
    case $opt in
        i) input_dir="$OPTARG"
        ;;
        o) output_dir="$OPTARG"
        ;;
        p) preview_mode="true"
        ;;
        f) force_overwrite="true"
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

if [ "$force_overwrite" == "true" ]; then
    echo "Force overwrite enabled. Existing files in $output_dir will be overwritten."
else
    echo "Force overwrite not enabled. Existing files in $output_dir will be preserved."
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$force_overwrite" == "true" ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/data.csv -o ${output_dir}/data.csv $preview
else
    echo "Skipping data.csv collection"
fi
if [ "$force_overwrite" == "true" ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best.csv -o ${output_dir}/best.csv $preview
else
    echo "Skipping best.csv collection"
fi
if [ "$force_overwrite" == "true" ]; then
    python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best_acceptance_corrected.csv -o ${output_dir}/best_acceptance_corrected.csv $preview
else
    echo "Skipping best_acceptance_corrected.csv collection"
fi


shopt -s nullglob # avoid literal wildcard if no matches

# handle optional files
files=(${input_dir}/mass*/best_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/best_projected_moments.csv -o ${output_dir}/moments.csv $preview
    else
        echo "Skipping best_projected_moments.csv collection"
    fi
fi

files=(${input_dir}/mass*/rand/rand.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand.csv -o ${output_dir}/rand.csv $preview
    else
        echo "Skipping rand.csv collection"
    fi
fi

files=(${input_dir}/mass*/rand/rand_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand_acceptance_corrected.csv -o ${output_dir}/rand_acceptance_corrected.csv $preview
    else
        echo "Skipping rand_acceptance_corrected.csv collection"
    fi
fi

files=(${input_dir}/mass*/rand/rand_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/rand/rand_projected_moments.csv -o ${output_dir}/rand_moments.csv $preview
    else
        echo "Skipping rand_projected_moments.csv collection"
    fi
fi

files=(${input_dir}/mass*/bootstrap/bootstrap.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap.csv -o ${output_dir}/bootstrap.csv $preview
    else
        echo "Skipping bootstrap.csv collection"
    fi
fi

files=(${input_dir}/mass*/bootstrap/bootstrap_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap_acceptance_corrected.csv -o ${output_dir}/bootstrap_acceptance_corrected.csv $preview
    else
        echo "Skipping bootstrap_acceptance_corrected.csv collection"
    fi
fi

files=(${input_dir}/mass*/bootstrap/bootstrap_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/bootstrap/bootstrap_projected_moments.csv -o ${output_dir}/bootstrap_moments.csv $preview
    else
        echo "Skipping bootstrap_projected_moments.csv collection"
    fi
fi

files=(${input_dir}/mass*/truth/best.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best.csv -o ${output_dir}/truth.csv $preview
    else
        echo "Skipping truth.csv collection"
    fi
fi

files=(${input_dir}/mass*/truth/best_acceptance_corrected.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best_acceptance_corrected.csv -o ${output_dir}/truth_acceptance_corrected.csv $preview
    else
        echo "Skipping truth_acceptance_corrected.csv collection"
    fi
fi

files=(${input_dir}/mass*/truth/best_projected_moments.csv)
if [ ${#files[@]} -gt 0 ]; then
    if [ "$force_overwrite" == "true" ]; then
        python ${script_dir}/collect_csv.py -i ${input_dir}/mass*/truth/best_projected_moments.csv -o ${output_dir}/truth_projected_moments.csv $preview
    else
        echo "Skipping truth_projected_moments.csv collection"
    fi
fi