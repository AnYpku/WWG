#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import createJMECorrector
from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.muonScaleResProducer import *
from WWG_Module import *

import argparse
import re
import optparse

parser = argparse.ArgumentParser(description='baseline selection')
parser.add_argument('-f', dest='infile', default='', help='local file input')
parser.add_argument('-y', dest='year', default='2016', help='year of dataset')
parser.add_argument('-m', dest='mode', default='local', help='runmode local/condor')
parser.add_argument('-d', dest='isdata',action='store_true',default=False)
parser.add_argument('-e', dest='era',default='B')
args = parser.parse_args()

print "mode: ", args.mode
print "year: ", args.year
print "dataset_name: ", args.infile
print "data_era", args.era

from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *

if args.year=='2016':
   jmeCorrections_ak4_MC = createJMECorrector(isMC=True,dataYear="UL2016",runPeriod="B",jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)
   jmeCorrections_ak4_Data = createJMECorrector(isMC=False,dataYear="UL2016",runPeriod=args.era,jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)

if args.year=='2016pre':
   jmeCorrections_ak4_MC = createJMECorrector(isMC=True,dataYear="UL2016_preVFP",runPeriod="B",jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)
   jmeCorrections_ak4_Data = createJMECorrector(isMC=False,dataYear="UL2016_preVFP",runPeriod=args.era,jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)

if args.year=='2017':
   jmeCorrections_ak4_MC = createJMECorrector(isMC=True,dataYear="UL2017",runPeriod="B",jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)
   jmeCorrections_ak4_Data = createJMECorrector(isMC=False,dataYear="UL2017",runPeriod=args.era,jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=False,splitJER=False)

if args.year=='2018':
   jmeCorrections_ak4_MC = createJMECorrector(isMC=True,dataYear="UL2018",runPeriod="A",jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=True,splitJER=False)
   jmeCorrections_ak4_Data = createJMECorrector(isMC=False,dataYear="UL2018",runPeriod=args.era,jesUncert="Total",jetType="AK4PFchs",noGroom=False,metBranchName="MET",applySmearing=True,isFastSim=False,applyHEMfix=True,splitJER=False)

# classify input files
if args.infile:
    infilelist = [args.infile]
    jsoninput = None
    fwkjobreport = False

    if args.mode == 'condor':
        infilelist.append(search.getValidSite(args.file)+args.file) 
    else:
        infilelist = [args.infile]
else:
    from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis
    infilelist = inputFiles()
    jsoninput = runsAndLumis()
    fwkjobreport = True

if args.isdata:
       Modules = [countHistogramsModule(),jmeCorrections_ak4_Data(),WWG_Module()]
else:
       if args.year=='2016':
          Modules = [countHistogramsModule(),WWG_Module(),jmeCorrections_ak4_MC(),puWeight_UL2016()]
       if args.year=='2016pre':
          Modules = [countHistogramsModule(),WWG_Module(),jmeCorrections_ak4_MC(),puWeight_UL2016()]
       if args.year=='2017':
          Modules = [countHistogramsModule(),WWG_Module(),jmeCorrections_ak4_MC(),puWeight_UL2017()]
       if args.year=='2018':
          Modules = [countHistogramsModule(),WWG_Module(),jmeCorrections_ak4_MC(),puWeight_UL2018()]

if args.isdata and args.year=='2018' and args.era=='D' and 'MuonEG' in args.infile:

    print 'special treatment for MuonEG_Run2018D'
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Config as cms

    lumisToProcess = cms.untracked.VLuminosityBlockRange( LumiList.LumiList(filename="./Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt").getCMSSWString().split(',') )
    # print lumisToProcess

    runsAndLumis_special = {}

    for l in lumisToProcess:
        if "-" in l:
            start, stop = l.split("-")
            rstart, lstart = start.split(":")
            rstop, lstop = stop.split(":")
        else:
            rstart, lstart = l.split(":")
            rstop, lstop = l.split(":")
        if rstart != rstop:
            raise Exception(
                "Cannot convert '%s' to runs and lumis json format" % l)
        if rstart not in runsAndLumis_special:
            runsAndLumis_special[rstart] = []
        runsAndLumis_special[rstart].append([int(lstart), int(lstop)])

    jsoninput = runsAndLumis_special

p=PostProcessor(".",infilelist,
                branchsel="WWG_keep_and_drop.txt",
                modules = Modules,
                provenance=True,
                justcount=False,
                noOut=False,
                fwkJobReport=fwkjobreport, 
                jsonInput=jsoninput, 
                outputbranchsel = "WWG_output_branch.txt")
p.run()



