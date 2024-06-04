# Current file structure is to place all mass bin dirs within a specific t-bin
# dir. We only want to extract the 1.2-1.25 GeV mass bin for cross section, so to make a 
# csv of the relevant fit files, we need to place all the t_dirs and their
# needed files in a single mass bin dir

# this file will be deprecated if fitsToCsv becomes python and handles directories
# better

fit_dir="/lustre19/expphy/volatile/halld/home/kscheuer/ampToolsFits/omegapi/GlueXI/PARA_0-PARA_135-PERP_45-PERP_90/data/ver03/1m_1p_iso/"
mass_bin="mass_1.200-1.250"

for t_dir in ${fit_dir}t*/; do

    # make the corresponding dir name in our single mass bin folder
    t_bin=$(basename ${t_dir})
    if [ "$t_bin" = "t_0.40-0.50" ]; then # TEMP: not part of cross section calc
        continue
    fi
    mkdir -p ${fit_dir}${mass_bin}/${t_bin}/

    # copy the best.fit and the fit.pdf files to the previously made t-bin dir
    for b1_bin in ${t_dir}${mass_bin}/; do                
        cp ${b1_bin}best.fit ${fit_dir}${mass_bin}/${t_bin}/
        cp ${b1_bin}fit.pdf ${fit_dir}${mass_bin}/${t_bin}/
    done
done
