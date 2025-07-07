source /etc/profile.d/modules.csh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.csh /w/halld-scshelf2101/kscheuer/neutralb1/version.xml

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load cuda
module load mpi/openmpi-x86_64
setenv CUDA_INSTALL_PATH /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2
set path=($path /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/bin/)
set path=($path /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/lib64)
setenv FSROOT /w/halld-scshelf2101/kscheuer/my_build/FSRoot