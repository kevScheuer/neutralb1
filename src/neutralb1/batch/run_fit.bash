#!/bin/bash
start_time=$(date +%s)
echo -e "\nhost: \t\t\t\t$HOSTNAME\n"

# cleanup directory
rm -f ./*.fit 
rm -f ./bestFitPars.txt
rm -f ./vecps_fitPars.txt
rm -f ./*.pdf
rm -f ./*.csv
rm -f ./vecps_plot.root
rm -f ./vecps_plot_*.root
rm -f ./*.ni
rm -f ./rand/*
rm -f ./bootstrap/*

# ==== GET ALL THE FIT PARAMETERS USING OPTARG ====
usage() {
    echo "Usage $0 [-n] [-r] [-b] [-s] [-t] [-h]"
    echo "Options:"        
    echo " -n       number of randomized fits to perform"   
    echo " -r       Reaction name of fit"        
    echo " -b       Number of bootstrap fits to perform for uncertainty estimation"
    echo " -s       Optional seed to use for randomized fits"
    echo " -t       Optional cfg file that mimics generated signal MC, to fit and obtain truth info. Using this will override some other options."    
    echo " -h       Display this help message"
}

# variables labelled with "my" to avoid conflicts and not have to unset globals
while getopts ":o:n:r:b:s:t:h:" opt; do
    case "${opt}" in
    n)
        echo -e "num rand fits: \t\t$OPTARG\n"
        my_num_rand_fits=$OPTARG
    ;;
    r)
        echo -e "reaction: \t\t\t$OPTARG\n"
        my_reaction=$OPTARG
    ;;
    b)
        echo -e "num bootstrap fits: $OPTARG\n"
        my_num_bootstrap_fits=$OPTARG
    ;;
    s)
        echo -e "seed: \t\t\t\t$OPTARG\n"
        my_seed="-rs $OPTARG"
    ;;
    t)
        echo -e "is truth fit?: \t\t$OPTARG\n"
        my_truth_bool=$OPTARG
    ;;
    h) # help message
        usage
        exit 0
    ;;
    \?) # invalid option
        echo "Invalid option: -$OPTARG"
        usage
        exit 1
    ;;
    esac
done
# ==== END GETTING ARGUMENTS ====

source setup_gluex.sh version.xml
# following neutralb1 paths are hard-coded for now
source /w/halld-scshelf2101/kscheuer/neutralb1/.venv/bin/activate
export PATH="/w/halld-scshelf2101/kscheuer/neutralb1/build/release/bin:$PATH"

# print some details to log
pwd
ls -al
echo -e "check if GPU Card is active for GPU fits:\n"
nvidia-smi

### RUN FIT ###
echo -e "\n\n==================================================\n
Beginning fits \n\n\n\n"

# if AmpTools has been built for MPI, these libraries show up, meaning fitMPI should run
use_mpi=false
if [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_GPU_MPI.a ] \
|| [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_MPI.a ]; then
    use_mpi=true
fi

run_start_time=$(date +%s)
if [ $my_truth_bool -eq 1 ]; then
    fit -c fit.cfg -m 1000000 -s "bestFitPars.txt"
elif [ "$use_mpi" = true ]; then
    echo -e "\nCheck that needed commands resolve:\n"
    which fitMPI
    which nvcc
    which mpirun
    which mpicxx

    mpirun fitMPI -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars.txt" $my_seed
else
    fit -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars.txt" $my_seed
fi

run_end_time=$(date +%s)
echo "Fitting took: $((run_end_time - run_start_time)) seconds"

if ! [ -f "$my_reaction.fit" ]; then
    echo -e "\n\nError: $my_reaction.fit not found, assuming fit failed. Exiting."
    exit 1
fi

### FIT DIAGNOSTICS AND RESULTS ###
mv "$my_reaction".fit best.fit
vecps_start_time=$(date +%s)
vecps_plotter best.fit
hadd -f vecps_plot.root vecps_plot_*.root
vecps_end_time=$(date +%s)
echo "vecps_plotter took: $((vecps_end_time - vecps_start_time)) seconds"

# TODO: this file should be replaced by a python script using a csv of the rand fits
# if [ -z "$my_truth_file" ]; then
#     cwd=$(pwd)
#     root -l -n -b -q "$my_code_dir/rand_fit_diagnostic.C(\"$cwd\")" 
# fi

# get angular distributions for each orientation, skip if only one file
vecps_files=(./vecps_plot_*.root)
if [ ${#vecps_files[@]} -gt 1 ]; then
    for file in "${vecps_files[@]}"; do
        [ -e "$file" ] || continue    
        angle_plotter -f "$file" --gluex-style
        pol_angle=$(basename "$file" | sed -E 's/vecps_plot_(.*)\.root/\1/')
        mv angles.pdf "angles_$pol_angle.pdf"
    done
fi

# get angular distributions for all orientations together
angle_plotter -f ./vecps_plot.root --gluex-style

# get total data csv file
hadd -f all_data.root *Amplitude*.root
uv run convert_to_csv -i $(pwd)/all_data.root -o $(pwd)/data.csv

# convert best fit results to a csv file and its projected moments
uv run convert_to_csv -i $(pwd)/best.fit -o $(pwd)/best.csv
uv run convert_to_csv -i $(pwd)/best.fit -o $(pwd)/best_moments.csv --moments

# process random fits
if [ $my_num_rand_fits -ne 0 ]; then
    mv -f "$my_reaction"_*.fit rand/
    mv -f bestFitPars_*.txt rand/
    mv -f "$my_reaction".ni rand/    
    uv run convert_to_csv -i $(pwd)/rand/"$my_reaction"_*.fit -o $(pwd)/rand/rand.csv
    uv run convert_to_csv -i $(pwd)/rand/"$my_reaction"_*.fit -o $(pwd)/rand/rand_moments.csv --moments
    echo -e "\n\n==================================================\n
    Randomized fits have completed and been moved to the rand subdirectory\n\n"
fi

ls -al

## BOOTSTRAP FITS ##
if [ $my_num_bootstrap_fits -ne 0 ]; then
    all_bootstrap_start_time=$(date +%s)
    echo -e "\n\n==================================================\nBeginning bootstrap fits\n\n\n\n"

    # truth files are already set up to start at the "best" values
    if [ $my_truth_bool -eq 1 ]; then
        cp -f fit.cfg bootstrap_fit.cfg
    else
        # otherwise create special bootstrap cfg and set it to start at the best fit values
        cp -f fit.cfg bootstrap_fit.cfg
        echo "include bestFitPars.txt" >> bootstrap_fit.cfg
    fi

    # replace data reader with bootstrap version for just the "data" file    
    sed -i -e '/data/s/ROOTDataReader \([^ ]\+\)/ROOTDataReaderBootstrap \1 0/' bootstrap_fit.cfg

    # perform requested number N of bootstrap fits, using seeds 1,2,..N
    for ((i=1;i<=$my_num_bootstrap_fits;i++)); do
        echo -e "\nBOOTSTRAP FIT: $i\n"      
        # replace the bootstrap seed in the cfg  
        sed -i -E "s/(ROOTDataReaderBootstrap [^ ]+) [0-9]+/\1 $i/" bootstrap_fit.cfg        
        # run fit
        bootstrap_start_time=$(date +%s)
        if [ "$use_mpi" = true ]; then 
            mpirun fitMPI -c bootstrap_fit.cfg -m 1000000
        else
            fit -c bootstrap_fit.cfg -m 1000000
        fi
        mv "$my_reaction".fit "$my_reaction"_"$i".fit
        bootstrap_end_time=$(date +%s)
        echo "Bootstrap fit $i took: $((bootstrap_end_time - bootstrap_start_time)) seconds"
    done    

    ls -al 
    # convert bootstrap fits to csv and move them to a subdirectory
    mv -f "$my_reaction"_*.fit bootstrap/
    mv -f "$my_reaction".ni bootstrap/
    uv run convert_to_csv -i $(pwd)/bootstrap/"$my_reaction"_*.fit -o $(pwd)/bootstrap/bootstrap.csv
    uv run convert_to_csv -i $(pwd)/bootstrap/"$my_reaction"_*.fit -o $(pwd)/bootstrap/bootstrap_moments.csv --moments
    all_bootstrap_end_time=$(date +%s)
    echo "All bootstrap fits took: $((all_bootstrap_end_time - all_bootstrap_start_time)) seconds"
    echo -e "\n\n==================================================\nBootstrap fits have completed and been moved to the bootstrap subdirectory\n\n"
    ls -al
fi

# rename truth file outputs
if [ $my_truth_bool -eq 1 ]; then
    mv -f best.fit best_truth.fit    
fi

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Execution time: $elapsed_time seconds"