#! /bin/csh
mv -v qcd1/Merged.root QCD_HT300to470.root
mv -v qcd2/Merged.root QCD_HT470to600.root
mv -v qcd3/Merged.root QCD_HT600to800.root
mv -v qcd4/Merged.root QCD_HT800to1000.root
mv -v qcd5/Merged.root QCD_HT1000to1400.root
mv -v qcd6/Merged.root QCD_HT1400to1800.root
mv -v qcd7/Merged.root QCD_HT1800.root
#hadd qcd_all.root *.root
