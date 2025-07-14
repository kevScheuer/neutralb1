source /etc/profile.d/modules.sh
module purge

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Get the project root directory (one level up from config/)
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Source GlueX environment using relative path to version.xml
source /group/halld/Software/build_scripts/gluex_env_jlab.sh "$PROJECT_ROOT/config/version.xml"

# load necessary modules
module use /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/modulefiles
module load mpi/openmpi-x86_64

export PATH="$PROJECT_ROOT/build/bin:$PATH"

# add paths for c++ / ROOT to search for libraries and include files
# useful for compiling and linking outside of the project / makefile method
export LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${LIBRARY_PATH}
export CPLUS_INCLUDE_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${CPLUS_INCLUDE_PATH}
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${LD_LIBRARY_PATH}
export DYLD_LIBRARY_PATH=${ROOTSYS}/lib:${AMPTOOLS}/lib:${HALLD_SIM_HOME}/src/libraries:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/include:${HALLD_SIM_HOME}/Linux_Alma9-x86_64-gcc11.5.0/lib:${DYLD_LIBRARY_PATH}
