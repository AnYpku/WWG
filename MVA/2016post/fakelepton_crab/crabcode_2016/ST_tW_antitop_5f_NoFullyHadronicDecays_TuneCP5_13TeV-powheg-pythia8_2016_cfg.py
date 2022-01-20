from WMCore.Configuration import Configuration 

config = Configuration()
config.section_("General")
config.General.requestName = "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_2016"
config.General.transferLogs = False 
config.General.workArea = "crab2016"

config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "PSet.py"
config.JobType.scriptExe = "./WWG_crab_script.sh" 
config.JobType.inputFiles = ["../../../scripts/haddnano.py","../WWG_fakelepton/WWG_postproc.py","../WWG_fakelepton/WWGfakelepton_Module.py","../WWG_fakelepton/WWG_keep_and_drop.txt","../WWG_fakelepton/WWG_output_branch.txt","../WWG_fakelepton/DAS_filesearch.py"] #hadd nano will not be needed once nano tools are in cmssw 
config.JobType.scriptArgs = ["isdata=MC","year=2016"] 
config.JobType.sendPythonFolder  = True
config.JobType.allowUndistributedCMSSW = True 

config.section_("Data")
config.Data.inputDataset = "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM" 
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
config.Data.outputDatasetTag = "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_2016" 

config.section_("Site")
config.Site.storageSite = "T3_CH_CERNBOX"
config.Site.whitelist = ["T2_US_MIT","T2_US_Wisconsin","T2_US_Purdue","T2_US_UCSD","T2_US_Caltech","T2_US_Nebraska"] 
