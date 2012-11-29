#! /bin/csh

#to copy files from dcahe to FIU

set srcDir = "11/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac"
set desDir = "/share/store/users/samantha/52X_PATs/PT1800"

set fileList = `ssh cmslpc22.fnal.gov "dcls -1 /pnfs/cms/WAX/11/store/user/lpcsusyhad/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/"`

foreach file ( $fileList )
	echo "copying file $file"
	time /opt/gLite/glite-UI-3.2.8/d-cache/srm/bin/srmcp "srm://cmssrm.fnal.gov:8443/${srcDir}/${file}" "file://localhost/${desDir}/${file}"
end
