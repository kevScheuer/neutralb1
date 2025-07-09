source /etc/profile.d/modules.csh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.csh /w/halld-scshelf2101/kscheuer/neutralb1/version.xml

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load mpi/openmpi-x86_64
setenv FSROOT /w/halld-scshelf2101/kscheuer/my_build/FSRoot