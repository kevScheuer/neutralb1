# Combines multiple ROOT trees for different orientations into one data file
#
# When a fit is done over multiple orientations, a subdirectory will contain those 
# individual orientation files. For get_bin_data.C to work properly, we would like 
# to apply it to the TOTAL file (all orientations combined) that AmpTools is effectively
# using. The simplest way is to iterate over those directories of interest and create 
# an "ALL_{data_type}.root" file to run get_bin_data.C on
#
# NOTE: File currently only works when subdirectories of parent are of the form 
# "mass_{#}-{#}/{data_files_to_combine}"

PREFIX="AmpToolsInputTree_sum"

parent_dir=$1
data_type=$2 # suffix indicating data type, usually "data", or "ver_03.1_mc"

for d in $parent_dir* ; do        
    files=""
    for filename in "$d/${PREFIX}_P*_*_$data_type.root"; do        
        files="${files} ${filename}"
    done        
    hadd "$d/ALL_${data_type}.root" $files
done

