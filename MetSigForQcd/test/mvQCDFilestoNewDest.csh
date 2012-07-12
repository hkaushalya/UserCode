#! /bin/csh
#1. set the absolute path to new dir
#2. source this 

set newDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/FactorizationResults/FactMoreJetTest/3jet5dphi"
set condorDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/FactorizationResults/CondorJobs/MC/FactMoreJetTest/3jet5dphi"
set qcdDir = "$condorDir"

cp  -v $qcdDir/qcd1/Merged.root $newDir/QCD_HTbin1.root
cp  -v $qcdDir/qcd2/Merged.root $newDir/QCD_HTbin2.root
cp  -v $qcdDir/qcd3/Merged.root $newDir/QCD_HTbin3.root
cp  -v $qcdDir/qcd4/Merged.root $newDir/QCD_HTbin4.root
cp  -v $qcdDir/qcd5/Merged.root $newDir/QCD_HTbin5.root
cp  -v $qcdDir/qcd6/Merged.root $newDir/QCD_HTbin6.root
cp  -v $qcdDir/qcd7/Merged.root $newDir/QCD_HTbin7.root
cp  -v $qcdDir/qcd8/Merged.root $newDir/QCD_HTbin8.root
cp  -v $qcdDir/qcd9/Merged.root $newDir/QCD_HTbin9.root
