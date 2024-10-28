#!/bin/sh
echo -e "\nhost: \t\t\t\t$HOSTNAME\n"

# cleanup directory
rm ./*.fit 
rm ./bestFitPars.txt
rm ./vecps_fitPars.txt
rm ./rand/*.fit
rm ./rand/bestFitPars*.txt
rm ./rand/*.ni
rm ./rand/*.pdf
rm ./bootstrap/*.fit
rm ./bootstrap/*.ni
rm ./*.pdf
rm ./*.root
rm ./*.ni

# ==== GET ALL THE FIT PARAMETERS USING OPTARG ====
usage() {
    echo "Usage $0 [-o] [-r] [-n] [-d] [-D] [-p] [-P] [-s] [-C] [-R] [-b] [-t]"
    echo "Options:"
    echo " -o       polarization orientations with '-' delimeter" 
    echo " -r       GlueX run period"
    echo " -n       number of randomized fits to perform"
    echo " -d       version # of data trees"
    echo " -D       MC selector option for data trees"
    echo " -p       version # of phasespace trees"
    echo " -P       MC selector option for phasespace trees"    
    echo " -s       directory where data Source files are stored"
    echo " -C       Code directory where scripts are stored"
    echo " -R       Reaction name of fit"        
    echo " -b       Perform bootstrapping for uncertainty estimation"
    echo " -t       Use cfg that generated signal MC to fit and obtain truth info"
}

# variables labelled with "my" to avoid conflicts and not have to unset globals
while getopts ":o:r:n:d:D:p:P:s:C:R:b:t:h:" opt; do
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
    s)
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
        echo -e "bootstrapped: \t\t$OPTARG\n"
        my_bootstrap_bool=$OPTARG
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

if ! [ -z "$my_truth_file" ]; then
    fit -c $my_truth_file -m 1000000 -s "bestFitPars"
elif [ "$use_mpi" = true ]; then
    echo -e "\nCheck that needed commands resolve:\n"
    which fitMPI
    which nvcc
    which mpirun
    which mpicxx

    mpirun fitMPI -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars"
else
    fit -c fit.cfg -m 1000000 -r $my_num_rand_fits -s "bestFitPars"
fi

if ! [ -f "$my_reaction.fit" ]; then
    echo -e "\n\nError: $my_reaction.fit not found, assuming fit failed. Exiting."
    exit
fi

# run diagnostic plotters
mv "$my_reaction".fit best.fit
vecps_plotter best.fit

if [ -z "$my_truth_file" ]; then
    cwd=$(pwd)
    root -l -b -q "$my_code_dir/rand_fit_diagnostic.C(\"$cwd\")"
fi

# produce diagnostic plots of angles and mass distributions
if [[ $my_data_version == *"data"* ]]; then
    my_label="Data"
elif [[ $my_data_option == *"_mcthrown"* ]]; then
    my_label="Thrown MC"
else
    my_label="Recon MC"
fi
root -l -b -q "$my_code_dir/angle_plotter.C(\"vecps_plot.root\", \"\", \"$my_label\")"

ls -al

# move all the randomized fits to the rand directory
if [ $my_num_rand_fits -ne 0 ]; then
mv -f "$my_reaction"_*.fit rand/
mv -f bestFitPars_*.txt rand/
mv -f "$my_reaction".ni rand/
mv -f rand_fit_diagnostic.pdf rand/
echo -e "\n\n==================================================\n
Randomized fits have completed and been moved to the rand subdirectory\n\n"
ls -al
fi


# Perform bootstrapping if requested
if [ $my_bootstrap_bool -ne 0 ] && [ -z "$my_truth_file" ]; then
    echo -e "\n\n==================================================\nBeginning bootstrap fits\n\n\n\n"
    # create special bootstrap cfg and set it to start at the best fit values
    cp -f fit.cfg bootstrap_fit.cfg
    echo "include bestFitPars.txt" >> bootstrap_fit.cfg

    # replace data reader with bootstrap version for just the "data" file
    sed -i -e 's/data LOOPREAC ROOTDataReader LOOPDATA/data LOOPREAC ROOTDataReaderBootstrap LOOPDATA 0/' bootstrap_fit.cfg

    # perform requested number N of bootstrap fits, using seeds 1,2,..N
    for ((i=1;i<=$my_bootstrap_bool;i++)); do
        echo -e "\nBOOTSTRAP FIT: $i\n"        
        sed -i -e "s/ROOTDataReaderBootstrap LOOPDATA $((i-1))/ROOTDataReaderBootstrap LOOPDATA $i/" bootstrap_fit.cfg
        # run fit
        if [ "$use_mpi" = true ]; then 
            mpirun fitMPI -c bootstrap_fit.cfg -m 1000000
        else
            fit -c bootstrap_fit.cfg -m 1000000
        fi
        mv "$my_reaction".fit "$my_reaction"_"$i".fit
    done    

    ls -al 
    mv -f "$my_reaction"_*.fit bootstrap/
    mv -f "$my_reaction".ni bootstrap/
    echo -e "\n\n==================================================\nBootstrap fits have completed and been moved to the bootstrap subdirectory\n\n"
    ls -al
fi

# rename truth file outputs
if [ ! -z "$my_truth_file" ]; then
    mv -f best.fit best_truth.fit
    mv -f bestFitPars bestFitPars.txt # don't know why this only happens here
fi