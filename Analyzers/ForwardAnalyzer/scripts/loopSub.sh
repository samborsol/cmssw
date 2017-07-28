#! /bin/bash

####################################################################################
#                   Usage Notes                                                    #
#   ./loopSub.sh InputTextFile NameOfOutput(DO NOT add .root at the end) NameOfcfg #
#                                                                                  #
# What it does: It takes your input text file, which has a list of root files from #
# EOS that this reads and makes temporary text files. Those temporary files are    #
# read in by ana.sh. Also the cfg you specify is the cfg.py that ana.sh uses. ana  #
# will submit jobs to LSF, one for each root file that is in the input text file.  #
####################################################################################

i=0
anahome=/afs/cern.ch/work/j/jgomez2/private/CMSSW_5_3_8_HI/src/Analyzers/ForwardAnalyzer
inputlist=$1
outfile=$2
anacfg=$3

for a in `cat $inputlist`
do
  echo $a>$anahome/tmp$i.txt
  
	bsub -q 1nh $anahome/ana.sh $anacfg $anahome/tmp$i.txt $outfile$i.root
	(( i++ ))
done
