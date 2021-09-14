from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "MuonEG_Run2018A"
config.General.transferLogs = False 
config.General.workArea = "crab2018"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../scripts/haddnano.py","../WWG_selector/WWG_postproc.py","../WWG_selector/WWG_Module.py","../WWG_selector/WWG_keep_and_drop.txt","../WWG_selector/WWG_output_branch.txt","../WWG_selector/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=data","year=2018","era=A"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/MuonEG/Run2018A-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD" 
config.Data.inputDBS = "global"
# config.Data.splitting = "LumiBased"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt" 

config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = "MuonEG_Run2018A" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 
