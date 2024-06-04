source /etc/profile.d/modules.csh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.csh /work/halld/kscheuer/my_build/version.xml

# # load mpi tools
# source /etc/profile.d/modules.csh
# module use /apps/modulefiles
# module load mpi/openmpi3-x86_64

# # load nvcc path for fitMPI
# module add cuda
# setenv CUDA_INSTALL_PATH /apps/cuda/11.4.2/
# set path=($path /apps/cuda/11.4.2/bin)
# set path=($path /apps/cuda/11.4.2/lib64/)


module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles

module load cuda
module load mpi/openmpi-x86_64
setenv CUDA_INSTALL_PATH /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2
set path=($path /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/bin/)
set path=($path /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2/lib64)
