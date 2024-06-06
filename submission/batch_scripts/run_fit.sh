#!/bin/sh
echo -e "\nhost: \t\t\t$HOSTNAME\n"

# cleanup local running directory
rm ./*

# get command line arguments with optarg 
usage() {
    echo "Usage $0 [-o] [-r] [-n] [-d] [-p] [-D] [-s] [-c] [-R] [-t]"
    echo "Options:"
    echo " -o       polarization orientations with '-' delimeter" 
    echo " -r       GlueX run period"
    echo " -n       number of randomized fits to perform"
    echo " -d       version # of data trees"
    echo " -p       version # of phasespace trees"
    echo " -D       Data output directory on \volatile"
    echo " -s       directory where data Source files are stored"
    echo " -c       directory where scripts are stored"
    echo " -R       reaction name of fit"
    echo " -t       (optional) truth file for thrown MC"
}

# variables labelled with "my" to avoid conflicts and not have to unset globals
while getopts ":ornvpdscxt:" opt; do
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
    p)
        echo -e "phasespace version: \t\t$OPTARG\n"
        my_phasespace_version=$OPTARG
    ;;
    D)
        echo -e "data out dir: \t\t$OPTARG\n"
        my_data_out_dir=$OPTARG
    ;;
    s)
        echo -e "source file dir: \t\t$OPTARG\n"
        my_source_file_dir=$OPTARG
    ;;
    c)
        echo -e "code dir: \t\t$OPTARG\n"
        my_code_dir=$OPTARG
    ;;
    R)
        echo -e "reaction: \t\t$OPTARG\n"
        my_reaction=$OPTARG
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

source $my_code_dir/setup_gluex.sh

# print some details to log
echo -e "check if GPU Card is active for GPU fits:\n"
nvidia-smi
pwd

# copy in needed files
rsync -aq $my_code_dir/fit.cfg ./
if [[ ! -z "$my_truth_file" ]]; then
    rsync -aq $my_code_dir/$my_truth_file ./
fi

# LINK FILES TO THIS RUNNING DIRECTORY AND THE OUTPUT DIRECTORY
tree="AmpToolsInputTree_sum"
amp_string="anglesOmegaPiAmplitude"
ph_string="anglesOmegaPiPhasespace"

# link data files
for ont in ${my_orientations}; do
    ont_num="$(cut -d'_' -f2 <<<"$ont")" # get the angle integer e.g. PARA_90 -> 90

    ln -sf ${my_source_file_dir}/${tree}_${ont}_${my_run_period}_${my_data_version}.root\
     ./${amp_string}_${ont_num}.root
    ln -sf ${my_source_file_dir}/${tree}_${ont}_${my_run_period}_${my_data_version}.root\
     ${my_data_out_dir}/${amp_string}_${ont_num}.root
done

# link generated phasespace
ln -sf ${my_source_file_dir}/${ph_string}Gen_${my_run_period}_${my_phasespace_version}.root\
 ./${ph_string}.root
ln -sf ${my_source_file_dir}/${ph_string}Gen_${my_run_period}_${my_phasespace_version}.root\
 ${my_data_out_dir}/${ph_string}.root

# link accepted phasespace
if [[ $my_data_version == *"_mcthrown"* ]]; then
    acc_string="Gen" # when using thrown MC, this means no detector effects are applied
else
    acc_string="Acc"
fi
ln -sf ${my_source_file_dir}/${ph_string}${acc_string}_${my_run_period}_${my_phasespace_version}.root\
 ./${ph_string}Acc.root
ln -sf ${my_source_file_dir}/${ph_string}${acc_string}_${my_run_period}_${my_phasespace_version}.root\
 ${my_data_out_dir}/${ph_string}Acc.root


ls -al 

### RUN FIT ###

# if AmpTools has been built for MPI, these libraries show up, meaning fitMPI should run
if [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_GPU_MPI.a ] \
|| [ -f $AMPTOOLS_HOME/AmpTools/lib/libAmpTools_MPI.a ]; then
    echo -e "\nCheck that needed commands resolve:\n"
    which fitMPI
    which nvcc
    which mpirun
    which mpicxx

    mpirun fitMPI -c fit.cfg -m 1000000 -r $my_num_rand_fits
else
    fit -c fit.cfg -m 1000000 -r $my_num_rand_fits
fi

if ! [ -f "$my_reaction.fit" ]; then
    echo "Error: $my_reaction.fit not found, assuming fit failed. Exiting."
    exit
fi

# run vecps_plotter and plot_plotter to get diagnostic plots
mv $my_reaction.fit best.fit
vecps_plotter best.fit

# produce diagnostic plots of angles and mass distributions
if [[ $my_data_version == *"data"* ]]; then
    data_type="Data"
elif [[ $my_data_version == *"_mc" ]]; then
    data_type="Recon MC"
elif [[ $my_data_version == *"_mcthrown" ]]; then
    data_type="Thrown MC"
fi
root -l -b -q "$my_code_dir/plot_plotter.C(\"vecps_plot.root\", \"\", \"$data_type\")"

ls -al

# cleanup data out directory 
rm $my_data_out_dir/*.fit

# move fit results to output directory (force overwrite)
cp -f fit.cfg $my_data_out_dir
mv -f best.fit $my_data_out_dir
mv -f omegapi_*.fit $my_data_out_dir
mv -f omegapi.ni $my_data_out_dir
mv -f vecps_* $my_data_out_dir
mv -f fit.pdf $my_data_out_dir