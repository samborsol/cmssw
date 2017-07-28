# rename file to crab.cfg before usage
[CRAB]

jobtype = cmssw
#scheduler = glidein
###       or let crab chose one server automatically for you 
#use_server = 1
use_server = 0
scheduler = remoteGlidein

[CMSSW]

dbs_url=http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_01/servlet/DBSServlet
datasetpath=/HIAllPhysics/appeltel-Flow_Skim_Light4_Run2010_HIAllPhysics-batch1-c50c3a75d9ab32e1e182d6f87aba0b48/USER
runselection=151076-152705
lumi_mask=Cert_150436-152957_HI7TeV_StreamExpress_Collisions10_JSON_v2.txt
#total_number_of_events=-1
#total_number_of_events=10
#events_per_job = 1000
#number_of_jobs = 40
total_number_of_lumis=-1
lumis_per_job = 25

output_file = ForwardTrees_skim_alltracks_2010.root
get_edm_output = 1

### The ParameterSet you want to use
pset=forwardanalyzer_2011data_cfg.py

[USER]
return_data = 0


### OUTPUT files INTO A SE
copy_data = 1
#return_data = 1
### if you want to copy data in a "official CMS site"
### you have to specify the name as written in
storage_element = T3_US_UMD
### the user_remote_dir will be created under the SE mountpoint
### in the case of publication this directory is not considered
### make sure the directory has proper permissions
### https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrabFaq#4_Stage_out_in_your_own_director
user_remote_dir = ForwardTrees/2010/PanicTime/
# Setting to fix the "missing environment" error:
check_user_remote_dir=0


##  Black and White Lists management:
## By Storage
se_black_list = T0,T1
#se_white_list =

## By ComputingElement
#ce_black_list =
#ce_white_list =
