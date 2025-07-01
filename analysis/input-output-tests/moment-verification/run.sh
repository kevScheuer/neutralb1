#!/bin/bash
# Perform the generation of signal MC, and fit with amplitudes and moments to it
# NOTE: this assumes you are on a GPU node

WORKSPACE_DIR="/home/kscheuer/work/neutralb1/"
source $WORKSPACE_DIR/setup_gluex.sh


# GENERATE MC
gen_vec_ps -c vec_ps_moment.cfg\
    -o data.root\
    -l 1.20 -u 1.22\
    -n 10000\
    -a 8.2 -b 8.8\
    -tmin 0.1 -tmax 0.2
if [ -e "data.root" ]; then
    echo "Data generation successful."
else
    echo "Data generation failed."
    exit 1
fi

# FIT WITH AMPLITUDES
fit -c amplitudes.cfg -m 10000000 > output.log
if [ -e "omegapi.fit" ]; then
    echo "Fit successful."
else
    echo "Fit failed."
    exit 1
fi
vecps_plotter omegapi.fit    
$WORKSPACE_DIR/submission/batch_script/angle_plotter vecps_plot.root "Thrown MC" "omegapi" ./ true
python $WORKSPACE_DIR/submission/batch_scripts/convert_to_csv.py -i omegapi.fit -o result.csv
python $WORKSPACE_DIR/submission/batch_scripts/project_moments.py -i omegapi.fit -o projected_moments.csv

# TODO: Moment fits and extraction here