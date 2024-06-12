source /etc/profile.d/modules.sh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.sh /w/halld-scshelf2101/kscheuer/neutralb1/submission/batch_scripts/version.xml

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load mpi/openmpi-x86_64
module load cuda
export CUDA_INSTALL_PATH=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2