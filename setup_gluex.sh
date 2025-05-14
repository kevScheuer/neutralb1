source /etc/profile.d/modules.sh
module purge

source /group/halld/Software/build_scripts/gluex_env_jlab.sh /w/halld-scshelf2101/kscheuer/neutralb1/version.xml

# add paths for c++ / ROOT to search for libraries and include files
export LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${CPLUS_INCLUDE_PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${DYLD_LIBRARY_PATH}

module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load mpi/openmpi-x86_64
module load cuda
export CUDA_INSTALL_PATH=/cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/cuda/12.2.2
export FSROOT=/w/halld-scshelf2101/kscheuer/my_build/FSRoot