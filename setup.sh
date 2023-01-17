#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

alias cmssw1021='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_2_13'
alias cmssw1029='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_2_9'
alias cmssw1068='export SCRAM_ARCH=slc7_amd64_gcc700; export CMSSW_VERSION=CMSSW_10_6_8'
alias cmssw1134='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_11_3_4'
alias cmssw1220='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_12_2_0'
alias cmssw1221prev4='export SCRAM_ARCH=slc7_amd64_gcc900; export CMSSW_VERSION=CMSSW_12_1_0_prev4_ROOT624'

cmssw1220
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd - > /dev/null

echo 'Setup following ROOT'
which root

export PATH=$DIR:$PATH
export LD_LIBRARY_PATH=$HOME/Utils/NanoTools:$HOME/Utils/rooutil:$PWD/tools:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PWD/tools:$HOME/Utils

PATH=$(printf %s "$PATH" \
     | awk -vRS=: -vORS= '!a[$0]++ {if (NR>1) printf(":"); printf("%s", $0) }' )

