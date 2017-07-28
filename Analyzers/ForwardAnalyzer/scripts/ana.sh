#! /bin/bash

###################################################################
#               Usage Notes                                       #
#                                                                 #
# Usage: ./ana.sh AnalyzerCFGName InputRootFile NameOfOutput.root #
#                                                                 #
# What it does: ana.sh reads in the input root file and tells the #
# cfg to analyzer that root file. It also tells the cfg the name  #
# you want for the output root file. If someone uses just this    #
# script then you have to specify the name AND ADD (.root) at the #
# end.                                                            #
###################################################################

upchome=/afs/cern.ch/work/j/jgomez2/private/CMSSW_5_3_8_HI/src/Analyzers/ForwardAnalyzer/
upcana=$1
input=$2
output=$3

cp $upchome$input ./
cd $upchome
eval `scram runtime -sh`
cd -

cmsRun $upchome$upcana $input $output
scp -o "StrictHostKeyChecking no" $output lxplus418:/tmp/jgomez2/
