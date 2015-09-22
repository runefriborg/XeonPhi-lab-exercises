#!/bin/bash
# Launch the nbody application on a Xeon Phi coprocessor

NBODIES=262144	# number of bodies
VERBOSE=",verbose"

# micnativeloadex will upload the required libraries
# micnativeloadex will not use the last core!
export SINK_LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH
set -x
micnativeloadex ./nbody-orig.mic -a ${NBODIES} -e "KMP_AFFINITY=compact,granularity=fine${VERBOSE}" -d 0

# alternative command if NFS is working for $PWD and $MIC_LD_LIBRARY_PATH (usually /opt/intel/...)
#ssh mic0 "cd $PWD; LD_LIBRARY_PATH=$MIC_LD_LIBRARY_PATH KMP_AFFINITY=compact,granularity=fine${VERBOSE} ./nbody-orig.mic ${NBODIES}"
