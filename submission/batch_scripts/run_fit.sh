#!/bin/sh
echo -e "\nhost: \t\t\t$HOSTNAME\n"

# cleanup local running directory
rm ./*

# variables labelled with "my" to avoid conflicts
export my_orientations=$1
export my_run_period=$2
export my_num_rand_fits=$3

export my_data_version=$4
export my_phasespace_version=$5
export my_source_file_dir=$6
export my_data_out_dir=$7
export my_code_dir=$8

export my_low_mass=$9
export my_high_mass=${10}
export my_low_t=${11}
export my_high_t=${12}

export my_reaction=${13}
export truth_file=${14}

# replace "-" with " " in orientations so it can be handed off to overwrite_template
my_orientations=${my_orientations//[-]/ }

echo -e "orientations: \t\t$my_orientations\n
run period: \t\t$my_run_period\n
num rand fits: \t\t$my_num_rand_fits\n
data version: \t\t$my_data_version\n
phasespace version: \t$my_phasespace_version\n
source file dir: \t$my_source_file_dir\n
data out dir: \t\t$my_data_out_dir\n
code dir: \t\t$my_code_dir\n
low mass: \t\t$my_low_mass\n
high mass: \t\t$my_high_mass\n
low t: \t\t\t$my_low_t\n
high t: \t\t$my_high_t\n
reaction: \t\t$my_reaction\n
truth file: \t\t$truth_file\n
"

source $my_code_dir/setup_gluex.sh
alias python="/w/halld-scshelf2101/kscheuer/miniforge3/envs/neutralb1/bin/python"

# print some details to log
echo -e "check if GPU Card is active for GPU fits:\n"
nvidia-smi
pwd

# copy script to working directory and write config file that will be template 
# for each binned fit
rsync -aq $my_code_dir/bin_template.cfg ./
rsync -aq $my_code_dir/overwrite_template.py ./
if [[ ! -z "$truth_file" ]]; then
    rsync -aq $my_code_dir/$truth_file ./
fi

# replace template file with bin info, and link in data/phasespace files
python overwrite_template.py\
 --data $my_data_version $my_source_file_dir\
 --phasespace $my_phasespace_version $my_source_file_dir\
 --data_out_dir $my_data_out_dir\
 -o $my_orientations -r $my_run_period\
 --mass_range $my_low_mass $my_high_mass --t_range $my_low_t $my_high_t

ls -al

### run fit ###

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