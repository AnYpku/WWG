#!/bin/bash
ulimit -s unlimited
set -e
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/y/yian/work/GZ/CMSSW_10_6_20/src/PhysicsTools/NanoAODTools/WWG/2017/condor
python WWG_postproc.py -f $1 -y 2017
