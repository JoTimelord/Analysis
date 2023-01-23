#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

alias cmssw1021='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_2_13'
alias cmssw1029='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_2_9'
alias cmssw1068='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_6_8'
alias cmssw10620='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_6_20'
alias cmssw1134='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_11_3_4'
alias cmssw1220='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_12_2_0'
alias cmssw1221prev4='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_12_1_0_prev4_ROOT624'

cmssw1220
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd - > /dev/null

#source $HOME/Utils/rooutil/root.sh
#source $HOME/Utils/rooutil/thisrooutil.sh
export PATH=$DIR:$PATH

LIBPATHS="
$HOME/Utils/rapido/src
$HOME/Utils/NanoTools/NanoCORE
$CMSSW_BASE/lib/$SCRAM_ARCH
"

for LIBPATH in $LIBPATHS; do
    # Add library to LD_LIBRARY_PATH
    if [[ "$LD_LIBRARY_PATH" != *"$LIBPATH"* ]]; then
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBPATH
    fi
    # Add library to ROOT_INCLUDE_PATH
    if [[ "$ROOT_INCLUDE_PATH" != *"$LIBPATH"* ]]; then
        export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$LIBPATH
    fi
done

PATH=$(printf %s "$PATH" \
     | awk -vRS=: -vORS= '!a[$0]++ {if (NR>1) printf(":"); printf("%s", $0) }' )

