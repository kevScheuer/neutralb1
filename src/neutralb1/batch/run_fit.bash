#!/bin/bash
start_time=$(date +%s)
echo -e "\nhost: \t\t\t\t$HOSTNAME\n"

# cleanup directory
rm -f ./*.fit 
rm -f ./bestFitPars*.txt
rm -f ./fitPars*.txt
rm -f ./vecps_fitPars*.txt
rm -f ./*.pdf
rm -f ./*.csv
rm -f ./vecps_plot*.root
rm -f ./*.ni
rm -f ./rand/*
rm -f ./bootstrap/*
rm -f ./distributions/*
rm -f out.log
rm -f err.log

# =============== GET ALL THE FIT PARAMETERS USING OPTARG ===============
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
# =============== END GETTING ARGUMENTS ===============

# =============== SETUP ENVIRONMENT ===============
source setup_gluex.sh version.xml
# following neutralb1 paths are hard-coded for now
source /w/halld-scshelf2101/kscheuer/neutralb1/.venv/bin/activate
export PATH="/w/halld-scshelf2101/kscheuer/neutralb1/build/release/bin:$PATH"

# print some details to log
pwd
ls -al
echo -e "check if GPU Card is active for GPU fits:\n"
nvidia-smi

# if AmpTools has been built for MPI, these libraries show up, meaning fitMPI should run
use_mpi=false
if [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_GPU_MPI.a ] \
|| [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_MPI.a ]; then
    use_mpi=true
fi

#=============== PERFORM PRIMARY FITS ===============
# pwa fits
echo -e "\n\n==================================================\n
Beginning Amplitude fits \n\n\n\n"
amplitude_start_time=$(date +%s)
if [ $my_truth_bool -eq 1 ]; then
    fit -c fit.cfg -m 1000000 -s "fitPars.txt"
elif [ "$use_mpi" = true ]; then
    echo -e "\nCheck that needed commands resolve:\n"
    which fitMPI
    which nvcc
    which mpirun
    which mpicxx

    mpirun --mca btl ^openib --display-map --display-allocation fitMPI -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "fitPars.txt" $my_seed
else
    fit -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "fitPars.txt" $my_seed
fi

amplitude_end_time=$(date +%s)
echo "Amplitude fitting took: $((amplitude_end_time - amplitude_start_time)) seconds"

# moment fits
# these can be run in the same directory because output is <reaction>_moment.fit
# echo -e "\n\n==================================================\n
# Beginning Moment fits \n\n\n\n"
# moment_start_time=$(date +%s)
# if [ "$use_mpi" = true ]; then
#     mpirun --mca btl ^openib --display-map --display-allocation fitMPI -c fit_moment.cfg -m 1000000 -r $my_num_rand_fits -s "fitPars_moment.txt" $my_seed
# else
#     fit -c fit_moment.cfg -m 1000000 -r $my_num_rand_fits -s "fitPars_moment.txt" $my_seed
# fi
# moment_end_time=$(date +%s)
# echo "Moment fitting took: $((moment_end_time - moment_start_time)) seconds"

# Quit job if fit failed
if ! [ -f "$my_reaction.fit" ]; then
    echo -e "\n\nError: $my_reaction.fit not found, assuming fit failed. Exiting."
    exit 1
fi
# if ! [ -f "${my_reaction}_moment.fit" ]; then
#     echo -e "\n\nError: ${my_reaction}_moment.fit not found, assuming fit failed. Exiting."
#     exit 1
# fi

# =============== VECPS PLOTTER AND ANGULAR DISTRIBUTIONS ===============
mv "$my_reaction".fit best.fit
mv fitPars.txt bestFitPars.txt
# mv "${my_reaction}_moment.fit" best_moment.fit
# mv fitPars_moment.txt bestFitPars_moment.txt

# This section will create vecps_plot.root files with tons of angular distribution
# histograms, and combine the ones of primary interest into pdfs
vecps_start_time=$(date +%s)

# amplitudes
vecps_plotter best.fit
hadd -f vecps_plot.root vecps_plot_*.root
vecps_files=(./vecps_plot_*.root)
for file in "${vecps_files[@]}"; do
    [ -e "$file" ] || continue    
    angle_plotter -f "$file" --gluex-style
    # name the pdf files according to the reaction
    pol_angle=$(basename "$file" | sed -E 's/vecps_plot_(.*)\.root/\1/')
    mv angles.pdf "angles_$pol_angle.pdf"
done
angle_plotter -f ./vecps_plot.root --gluex-style
mv vecps_plot*.root ./distributions/
mv *.pdf ./distributions/

# moments
# vecps_plotter best_moment.fit
# hadd -f vecps_plot_moment.root vecps_plot_*.root
# vecps_moment_files=(./vecps_plot_*.root)
# for file in ${vecps_moment_files[@]}; do
#     # avoid running this on the combined file
#     if [ "$(basename "$file")" = "vecps_plot_moment.root" ]; then
#         continue
#     fi
#     [ -e "$file" ] || continue
#     angle_plotter -f "$file" --gluex-style
#     pol_angle=$(basename "$file" | sed -E 's/vecps_plot_(.*)\.root/\1/')
#     mv angles.pdf "angles_moment_$pol_angle.pdf"
#     # now rename the root files to not interfere with the amplitude ones
#     mv "$file" "$(basename "$file" .root)_moment.root"
# done
# angle_plotter -f ./vecps_plot_moment.root --gluex-style
# mv angles.pdf angles_moment.pdf
# mv vecps_plot*.root ./distributions/
# mv *.pdf ./distributions/

vecps_end_time=$(date +%s)
echo "vecps_plotter took: $((vecps_end_time - vecps_start_time)) seconds"

# =============== CONVERT RESULTS TO CSV ===============
csv_start_time=$(date +%s)

# convert best fit results to csvs
# first conversion include data file creation, and MINUIT covariance/correlation matrices
uv run convert_to_csv -i $(pwd)/best.fit -o $(pwd)/best.csv --covariance --correlation --verbose
uv run convert_to_csv -i $(pwd)/best.fit -o $(pwd)/best_acceptance_corrected.csv -a --no-data --verbose
uv run convert_to_csv -i $(pwd)/best.fit -o $(pwd)/best_projected_moments.csv --moments --no-data --verbose
# uv run convert_to_csv -i $(pwd)/best_moment.fit -o $(pwd)/best_moment.csv

# process random fits into subdirectory
if [ $my_num_rand_fits -ne 0 ]; then
    mv -f "$my_reaction"_*.fit rand/
    mv -f fitPars_*.txt rand/    
    uv run convert_to_csv -i $(ls rand/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/rand/rand.csv --no-data --verbose
    uv run convert_to_csv -i $(ls rand/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/rand/rand_acceptance_corrected.csv -a --no-data --verbose
    uv run convert_to_csv -i $(ls rand/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/rand/rand_projected_moments.csv --moments --no-data --verbose
    # uv run convert_to_csv -i $(pwd)/rand/"$my_reaction"_moment*.fit -o $(pwd)/rand/rand_moment.csv    
    echo -e "\n\n==================================================\n
    Randomized fits have completed and been moved to the rand subdirectory\n\n"
fi

ls -al

csv_end_time=$(date +%s)
echo "producing csv files took: $((csv_end_time - csv_start_time)) seconds"

# =============== BOOTSTRAP FITS ===============
if [ $my_num_bootstrap_fits -ne 0 ]; then    
    echo -e "\n\n==================================================\nBeginning bootstrap fits\n\n\n\n"

    # truth files are already set up to start at the "best" values
    if [ $my_truth_bool -eq 1 ]; then
        cp -f fit.cfg bootstrap_fit.cfg
    else
        # otherwise create special bootstrap cfg and set it to start at the best fit values
        cp -f fit.cfg bootstrap_fit.cfg
        echo "include bestFitPars.txt" >> bootstrap_fit.cfg
    fi

    # # do the same for the moments, but since truth not setup for it, always run
    # cp -f fit_moment.cfg bootstrap_fit_moment.cfg
    # echo "include bestFitPars_moment.txt" >> bootstrap_fit_moment.cfg

    # replace data reader with bootstrap version for just the "data" file    
    sed -i -e '/data/s/FSROOTDataReader \([^ ]\+\) \([^ ]\+\) \([^ ]\+\)/FSRootDataReaderBootstrap \1 \2 \3 0/' bootstrap_fit.cfg
    # sed -i -e '/data/s/FSROOTDataReader \([^ ]\+\) \([^ ]\+\) \([^ ]\+\)/FSRootDataReaderBootstrap \1 \2 \3 0/' bootstrap_fit_moment.cfg
    # TODO: do we want to bootstrap background as well?

    all_amplitude_bootstrap_start_time=$(date +%s)
    # perform requested number N of bootstrap fits, using seeds 1,2,..N
    for ((i=1;i<=$my_num_bootstrap_fits;i++)); do
        echo -e "\nAMPLITUDE BOOTSTRAP FIT: $i\n"      
        # replace the bootstrap seed in the cfg  
        sed -i -E "s/(FSRootDataReaderBootstrap [^ ]+ [^ ]+ [^ ]+) [0-9]+/\1 $i/" bootstrap_fit.cfg        
        # run fit
        bootstrap_start_time=$(date +%s)
        if [ "$use_mpi" = true ]; then 
            mpirun --mca btl ^openib --display-map --display-allocation fitMPI -c bootstrap_fit.cfg -m 1000000
        else
            fit -c bootstrap_fit.cfg -m 1000000
        fi
        mv "$my_reaction".fit "$my_reaction"_"$i".fit
        bootstrap_end_time=$(date +%s)
        echo "Amplitude bootstrap fit $i took: $((bootstrap_end_time - bootstrap_start_time)) seconds"
    done    
    all_amplitude_bootstrap_end_time=$(date +%s)
    echo "All amplitude bootstrap fits took: $((all_amplitude_bootstrap_end_time - all_amplitude_bootstrap_start_time)) seconds"

    # run bootstrap fits for moments
    # all_moment_bootstrap_start_time=$(date +%s)
    # for ((i=1;i<=$my_num_bootstrap_fits;i++)); do
    #     echo -e "\nMOMENT BOOTSTRAP FIT: $i\n"              
    #     sed -i -E "s/(ROOTDataReaderBootstrap [^ ]+) [0-9]+/\1 $i/" bootstrap_fit_moment.cfg        
    #     bootstrap_start_time=$(date +%s)
    #     if [ "$use_mpi" = true ]; then 
    #         mpirun --mca btl ^openib --display-map --display-allocation fitMPI -c bootstrap_fit.cfg -m 1000000
    #     else
    #         fit -c bootstrap_fit.cfg -m 1000000
    #     fi
    #     mv "${my_reaction}"_moment.fit "${my_reaction}"_moment_"$i".fit
    #     bootstrap_end_time=$(date +%s)
    #     echo "Moment bootstrap fit $i took: $((bootstrap_end_time - bootstrap_start_time)) seconds"
    # done    
    # all_moment_bootstrap_end_time=$(date +%s)
    # echo "All moment bootstrap fits took: $((all_moment_bootstrap_end_time - all_moment_bootstrap_start_time)) seconds"

    ls -al 
    # process bootstrap fits into subdirectory
    mv -f "$my_reaction"_*.fit bootstrap/
    mv -f "$my_reaction"_.ni bootstrap/
    mv -f bootstrap_fit.cfg bootstrap/
    mv -f bootstrap_fit_moment.cfg bootstrap/
    uv run convert_to_csv -i $(ls bootstrap/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/bootstrap/bootstrap.csv --no-data --verbose
    uv run convert_to_csv -i $(ls bootstrap/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/bootstrap/bootstrap_acceptance_corrected.csv -a --no-data --verbose
    uv run convert_to_csv -i $(ls bootstrap/"$my_reaction"_*.fit | grep -v '_moment') -o $(pwd)/bootstrap/bootstrap_projected_moments.csv --moments --no-data --verbose
    # uv run convert_to_csv -i $(pwd)/bootstrap/"$my_reaction"_moment*.fit -o $(pwd)/bootstrap/bootstrap_moment.csv
    
    echo -e "\n\n==================================================\nBootstrap fits have completed and been moved to the bootstrap subdirectory\n\n"
    ls -al
fi

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Execution time: $elapsed_time seconds"

# symlink the log file here for easy access
current_dir=$(pwd)
common_path=${current_dir#*"$USER"/}

ln -s "/farm_out/$USER/$common_path/log/out.log" ./out.log
ln -s "/farm_out/$USER/$common_path/log/err.log" ./err.log