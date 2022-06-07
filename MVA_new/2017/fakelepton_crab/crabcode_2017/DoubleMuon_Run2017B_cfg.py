from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "DoubleMuon_Run2017B_mva"
config.General.transferLogs = False 
config.General.workArea = "crab2017"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.maxJobRuntimeMin = 2630
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../../scripts/haddnano.py","../WWG_fakelepton/WWG_postproc.py","../WWG_fakelepton/WWGfakelepton_Module.py","../WWG_fakelepton/WWG_keep_and_drop.txt","../WWG_fakelepton/WWG_output_branch.txt","../WWG_fakelepton/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=data","year=2017","era=B"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/DoubleMuon/Run2017B-UL2017_MiniAODv2_NanoAODv9-v1/NANOAOD" 
config.Data.inputDBS = "global"
# config.Data.splitting = "FileBased"
# config.Data.unitsPerJob = 1
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 50
config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt" 

config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = "DoubleMuon_Run2017B_mva" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 
