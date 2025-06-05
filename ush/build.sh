#!/bin/bash

set -eux

# Get the root of the cloned GSI directory
readonly DIR_ROOT=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )/.." && pwd -P)

# User Options
BUILD_TYPE=${BUILD_TYPE:-"Release"}
CMAKE_OPTS=${CMAKE_OPTS:-}
COMPILER=${COMPILER:-"intel"}
BUILD_DIR=${BUILD_DIR:-"${DIR_ROOT}/build"}
INSTALL_PREFIX=${INSTALL_PREFIX:-"${DIR_ROOT}/install"}
GSI_MODE=${GSI_MODE:-"Regional"}  # By default build Regional GSI (for regression testing)
ENKF_MODE=${ENKF_MODE:-"GFS"}     # By default build Global EnKF  (for regression testing)
REGRESSION_TESTS=${REGRESSION_TESTS:-"YES"} # Build regression test suite

#==============================================================================#

# Detect machine (sets MACHINE_ID)
source $DIR_ROOT/ush/detect_machine.sh

# Load modules
set +x
source $DIR_ROOT/ush/module-setup.sh
module use $DIR_ROOT/modulefiles
module load "gsi_${MACHINE_ID}.${COMPILER}"
module list
set -x

# Using self-compiled CRTM example
CRTM_VERSION="2.3"
CRTM_PATH=/scratch2/GFDL/gfdlscr/Mingjing.Tong/CRTM
CRTMPATH="$DIR_ROOT/../../../../CRTM"
if [[ ${CRTM_VERSION} == "v3" ]]; then
   CRTMVER="CRTMv3"
elif [[ ${CRTM_VERSION} == "2.3" ]]; then
   CRTMVER="REL-2.3.0_emc"
fi
if [[ ${CRTM_VERSION} != "2.4" ]]; then
   export CRTM_LIB=${CRTMPATH}/${CRTMVER}/install/lib/libcrtm.a
   export CRTM_INC=${CRTMPATH}/${CRTMVER}/install/include
fi

# Set CONTROLPATH variable to user develop installation
CONTROLPATH="$DIR_ROOT/../develop/install/bin"
# Collect BUILD Options
CMAKE_OPTS+=" -DCMAKE_BUILD_TYPE=$BUILD_TYPE"

# Install destination for built executables, libraries, CMake Package config
CMAKE_OPTS+=" -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX"

# Configure for GSI and EnKF
CMAKE_OPTS+=" -DGSI_MODE=$GSI_MODE -DENKF_MODE=${ENKF_MODE}"

# Build regression test suite (on supported MACHINE_ID where CONTROLPATH exists)
[[ ${REGRESSION_TESTS} =~ [yYtT] ]] && CMAKE_OPTS+=" -DBUILD_REG_TESTING=ON -DCONTROLPATH=${CONTROLPATH:-}"

# Re-use or create a new BUILD_DIR (Default: create new BUILD_DIR)
[[ ${BUILD_CLEAN:-"YES"} =~ [yYtT] ]] && rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR && cd $BUILD_DIR

# Configure, build, install
cmake $CMAKE_OPTS $DIR_ROOT
make -j ${BUILD_JOBS:-8} VERBOSE=${BUILD_VERBOSE:-1}
make install

exit
