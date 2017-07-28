from WMCore.Configuration import Configuration
config = Configuration()
#############################
config.section_('General')
config.General.workArea = 'taskManagement'
config.General.requestName = 'RAPGAPwithEOS'
config.General.transferOutputs = True
config.General.transferLogs = True
#############################
config.section_('JobType')
config.JobType.psetName = 'readHepMC_cfg.py'
config.JobType.pluginName = 'Analysis'
############################
config.section_('Data')
#config.Data.inputDataset = '/HIForward/HIRun2015-02May2016-v1/AOD'
config.Data.userInputFiles = open('inputFilesRapgap_v2.txt').readlines()
#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/HI/Cert_262548-263757_PromptReco_HICollisions15_JSON_MuonPhys_v2.txt'
#config.Data.lumiMask = 'notFinishedLumis.json'

config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.unitsPerJob = 1
#config.Data.splitting = 'LumiBased'
#NJOBS = 1
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = True

#################################
config.section_('User')
#################################
config.section_('Site')
config.Site.storageSite = 'T2_US_Wisconsin'
