#!/bin/sh
start_time=$(date +%s)
echo -e "\nhost: \t\t\t\t$HOSTNAME\n"

# cleanup directory
rm ./*.fit 
rm ./bestFitPars.txt
rm ./vecps_fitPars.txt
rm ./*.pdf
rm ./*.csv
rm ./*.root
rm ./*.ni
rm ./rand/*.fit
rm ./rand/bestFitPars*.txt
rm ./rand/*.ni
rm ./rand/*.pdf
rm ./rand/*.csv
rm ./bootstrap/*.fit
rm ./bootstrap/*.ni
rm ./bootstrap/*.csv

# ==== GET ALL THE FIT PARAMETERS USING OPTARG ====
usage() {
    echo "Usage $0 [-o] [-r] [-n] [-d] [-D] [-p] [-P] [-S] [-C] [-R] [-b] [-s] [-t] [-h]"
    echo "Options:"
    echo " -o       polarization orientations with '-' delimeter" 
    echo " -r       GlueX run period"
    echo " -n       number of randomized fits to perform"
    echo " -d       version # of data trees"
    echo " -D       MC selector option for data trees"
    echo " -p       version # of phasespace trees"
    echo " -P       MC selector option for phasespace trees"    
    echo " -S       directory where data Source files are stored"
    echo " -C       Code directory where scripts are stored"
    echo " -R       Reaction name of fit"        
    echo " -b       Number of bootstrap fits to perform for uncertainty estimation"
    echo " -s       Optional seed to use for randomized fits"
    echo " -t       Optional cfg file that mimics generated signal MC, to fit and obtain truth info. Using this will override some other options."    
    echo " -h       Display this help message"
}

# variables labelled with "my" to avoid conflicts and not have to unset globals
while getopts ":o:r:n:d:D:p:P:S:C:R:b:s:t:h:" opt; do
    case "${opt}" in
    o)
        echo -e "orientations: \t\t$OPTARG\n"        
        # replace "-" with " " in orientations so it can be looped over
        my_orientations=${OPTARG//[-]/ }
    ;;
    r)
        echo -e "run period: \t\t$OPTARG\n"
        my_run_period=$OPTARG
    ;;
    n)
        echo -e "num rand fits: \t\t$OPTARG\n"
        my_num_rand_fits=$OPTARG
    ;;
    d)
        echo -e "data version: \t\t$OPTARG\n"
        my_data_version=$OPTARG
    ;;
    D)
        echo -e "data MC option: \t$OPTARG\n"
        my_data_option=$OPTARG
    ;;
    p)
        echo -e "phasespace version: $OPTARG\n"
        my_phasespace_version=$OPTARG
    ;;
    P)
        echo -e "phasespace MC option: $OPTARG\n"
        my_phasespace_option=$OPTARG
    ;;
    S)
        echo -e "source file dir: \t$OPTARG\n"
        my_source_file_dir=$OPTARG
    ;;
    C)
        echo -e "code dir: \t\t\t$OPTARG\n"
        my_code_dir=$OPTARG
    ;;
    R)
        echo -e "reaction: \t\t\t$OPTARG\n"
        my_reaction=$OPTARG
    ;;
    b)
        echo -e "num bootstrap fits: \t$OPTARG\n"
        my_num_bootstrap_fits=$OPTARG
    ;;
    s)
        echo -e "seed: \t\t\t\t$OPTARG\n"
        my_seed="-rs $OPTARG"
    ;;
    t)
        echo -e "truth file: \t\t$OPTARG\n"
        my_truth_file=$OPTARG
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

source $my_code_dir/setup_gluex.sh
source /home/kscheuer/.bashrc # loading conda env allows for running python scripts below
conda activate neutralb1

# print some details to log
echo -e "check if GPU Card is active for GPU fits:\n"
nvidia-smi
pwd

### LINK ROOT FILES ###
tree="AmpToolsInputTree_sum"
amp_string="anglesOmegaPiAmplitude"
ph_string="anglesOmegaPiPhaseSpace"

# link data files
for ont in ${my_orientations}; do
    ont_num="$(cut -d'_' -f2 <<<"$ont")" # get the angle integer e.g. PARA_90 -> 90

    ln -sf ${my_source_file_dir}/${tree}_${ont}_${my_run_period}_${my_data_version}${my_data_option}.root\
     ./${amp_string}_${ont_num}.root
done

# link generated phasespace
ln -sf ${my_source_file_dir}/${ph_string}Gen_${my_run_period}_${my_phasespace_version}.root\
 ./${ph_string}.root

# link accepted phasespace
if [[ $my_data_option == *"_mcthrown"* ]]; then
    acc_string="Gen" # when using thrown MC, this means no detector effects are applied
    my_phasespace_option="" # remove option since it won't apply in this case
else
    acc_string="Acc"
fi
ln -sf ${my_source_file_dir}/${ph_string}${acc_string}_${my_run_period}_${my_phasespace_version}${my_phasespace_option}.root\
 ./${ph_string}Acc.root

ls -al 

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
if ! [ -z "$my_truth_file" ]; then
    fit -c $my_truth_file -m 1000000 -s "bestFitPars"
elif [ "$use_mpi" = true ]; then
    echo -e "\nCheck that needed commands resolve:\n"
    which fitMPI
    which nvcc
    which mpirun
    which mpicxx

    mpirun fitMPI -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars" $my_seed
else
    fit -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars" $my_seed
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
vecps_end_time=$(date +%s)
echo "vecps_plotter took: $((vecps_end_time - vecps_start_time)) seconds"

if [ -z "$my_truth_file" ]; then
    cwd=$(pwd)
    root -l -n -b -q "$my_code_dir/rand_fit_diagnostic.C(\"$cwd\")" # TODO: replace with python script
fi

# produce diagnostic plots of angles and mass distributions
if [[ $my_data_version == *"data"* ]]; then
    my_label="Data"
elif [[ $my_data_option == *"_mcthrown"* ]]; then
    my_label="Thrown MC"
else
    my_label="Recon MC"
fi

$my_code_dir/angle_plotter vecps_plot.root $my_label $my_reaction ./ true

# get total data csv file
hadd all_data.root ${amp_string}_*.root
python $my_code_dir/convert_to_csv.py -i $(pwd)/all_data.root -o $(pwd)/data.csv

# convert best fit results to a csv file and moments
python $my_code_dir/convert_to_csv.py -i $(pwd)/best.fit -o $(pwd)/best.csv
breit_wigner_arg=""
if ! [ -z "$my_truth_file" ]; then
    breit_wigner_arg=" -b" # truth files use breit wigner parameters
fi
python $my_code_dir/project_moments.py -i $(pwd)/best.fit -o $(pwd)/best_moments.csv $breit_wigner_arg

# process random fits
if [ $my_num_rand_fits -ne 0 ] && [ -z "$my_truth_file" ]; then
    mv -f "$my_reaction"_*.fit rand/
    mv -f bestFitPars_*.txt rand/
    mv -f "$my_reaction".ni rand/
    mv -f rand_fit_diagnostic.pdf rand/
    python $my_code_dir/convert_to_csv.py -i $(pwd)/rand/"$my_reaction"_*.fit -o $(pwd)/rand/rand.csv    
    python $my_code_dir/project_moments.py -i $(pwd)/rand/"$my_reaction"_*.fit -o $(pwd)/rand/rand_moments.csv -w 10
    echo -e "\n\n==================================================\n
    Randomized fits have completed and been moved to the rand subdirectory\n\n"
fi

ls -al

## BOOTSTRAP FITS ##
if [ $my_num_bootstrap_fits -ne 0 ]; then
    all_bootstrap_start_time=$(date +%s)
    echo -e "\n\n==================================================\nBeginning bootstrap fits\n\n\n\n"

    # truth files are already set up to start at the "best" values
    if ! [ -z "$my_truth_file" ]; then
        cp -f $my_truth_file bootstrap_fit.cfg
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
    python $my_code_dir/convert_to_csv.py -i $(pwd)/bootstrap/"$my_reaction"_*.fit -o $(pwd)/bootstrap/bootstrap.csv
    python $my_code_dir/project_moments.py -i $(pwd)/bootstrap/"$my_reaction"_*.fit -o $(pwd)/bootstrap/bootstrap_moments.csv -w 10 $breit_wigner_arg
    all_bootstrap_end_time=$(date +%s)
    echo "All bootstrap fits took: $((all_bootstrap_end_time - all_bootstrap_start_time)) seconds"
    echo -e "\n\n==================================================\nBootstrap fits have completed and been moved to the bootstrap subdirectory\n\n"
    ls -al
fi

# rename truth file outputs
if [ ! -z "$my_truth_file" ]; then
    mv -f best.fit best_truth.fit
    mv -f bestFitPars bestFitPars.txt # don't know why this only happens here
fi

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Execution time: $elapsed_time seconds"