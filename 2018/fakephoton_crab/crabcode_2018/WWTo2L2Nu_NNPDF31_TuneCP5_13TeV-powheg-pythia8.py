from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8"
config.General.transferLogs = False 
config.General.workArea = "crab2018"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../scripts/haddnano.py","../WWG_fakephoton/WWG_postproc.py","../WWG_fakephoton/WWGfakephoton_Module.py","../WWG_fakephoton/WWG_keep_and_drop.txt","../WWG_fakephoton/WWG_output_branch.txt","../WWG_fakephoton/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=MC","year=2018"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM" 
#config.Data.inputDBS = "phys03"
config.Data.inputDBS = "global"
# config.Data.splitting = "LumiBased"
config.Data.splitting = "FileBased"
#config.Data.splitting = "EventAwareLumiBased" 
#config.Data.splitting = "Automatic" 
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.allowNonValidInputDataset = True
config.Data.outputDatasetTag = "WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 
