#!/bin/bash
# Perform the generation of signal MC, and fit with amplitudes and moments to it
# NOTE: this assumes you are on a GPU node, and are running from the same directory as
# the script

WORKSPACE_DIR="/home/kscheuer/work/neutralb1/"
source $WORKSPACE_DIR/setup_gluex.sh


# GENERATE MC
if [ -e "data.root" ]; then
    echo "Data file already exists, skipping generation."
else
    echo "Generating data..."
    gen_vec_ps -c amplitudes.cfg\
        -o data.root\
        -l 1.20 -u 1.22\
        -n 50000\
        -a 8.2 -b 8.8\
        -tmin 0.1 -tmax 0.2
    if [ -e "data.root" ]; then
        echo "Data generation successful."
    else
        echo "Data generation failed."
        exit 1
    fi
fi

# FIT WITH AMPLITUDES
fit -c amplitudes.cfg -m 10000000 -r 50 > amplitude_fit.log
if [ -e "omegapi.fit" ]; then
    echo "Fit successful."    
else
    echo "Fit failed."
    exit 1
fi
vecps_plotter omegapi.fit    
$WORKSPACE_DIR/submission/batch_scripts/angle_plotter vecps_plot.root "Thrown MC" "" ./ true
python $WORKSPACE_DIR/submission/batch_scripts/convert_to_csv.py -i omegapi.fit -o result.csv
python $WORKSPACE_DIR/submission/batch_scripts/project_moments.py\
    -i omegapi.fit -o projected_moments.csv -m 1.21

mkdir -p amplitude_results
mv omegapi_*.fit amplitude_results/
mv omegapi.fit amplitude_results/
mv omegapi.ni amplitude_results/
mv vecps_plot.root amplitude_results/
mv vecps_fitPars.txt amplitude_results/
mv result.csv amplitude_results/
mv projected_moments.csv amplitude_results/

# FIT WITH MOMENTS
fit -c moments.cfg -m 10000000 -r 50 > moment_fit.log
if [ -e "omegapi.fit" ]; then
    echo "Fit successful."    
else
    echo "Fit failed."
    exit 1
fi
vecps_plotter omegapi.fit
$WORKSPACE_DIR/submission/batch_script/angle_plotter vecps_plot.root "Thrown MC" "" ./ true
python $WORKSPACE_DIR/submission/batch_scripts/convert_to_csv.py -i omegapi.fit -o result.csv

mkdir -p moment_results
mv omegapi_*.fit moment_results/
mv omegapi.fit moment_results/
mv omegapi.ni moment_results/
mv vecps_plot.root moment_results/
mv vecps_fitPars.txt moment_results/
mv result.csv moment_results/