source /etc/profile.d/modules.sh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.sh /w/halld-scshelf2101/kscheuer/ambiguities/submission/version.xml

# # load mpi tools (el7)
# source /etc/profile.d/modules.sh
# module use /apps/modulefiles
# module load mpi/openmpi3-x86_64
# export CUDA_INSTALL_PATH=/apps/cuda/11.4.2/
# export PATH="/apps/cuda/11.4.2/bin:$PATH"
# export LD_LIBRARY_PATH="/apps/cuda/11.4.2/lib64:$LD_LIBRARY_PATH"

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load /usr/share/modulefiles/mpi/openmpi-x86_64
module load cuda
# export CUDA_INSTALL_PATH=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2
# export PATH="/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/bin/:$PATH"
# export LD_LIBRARY_PATH="/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/lib64:$LD_LIBRARY_PATH"