from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_2016pre_mva"
config.General.transferLogs = False 
config.General.workArea = "crab2016pre"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../../scripts/haddnano.py","../WWG_selector/WWG_postproc.py","../WWG_selector/WWG_Module.py","../WWG_selector/WWG_keep_and_drop.txt","../WWG_selector/WWG_output_branch.txt","../WWG_selector/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=MC","year=2016pre","era=B"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM" 
config.Data.inputDBS = "global"
# config.Data.splitting = "FileBased"
# config.Data.unitsPerJob = 1
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
#config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = "ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_2016pre_mva" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
#config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 