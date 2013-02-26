#! /bin/csh

set curDir = `pwd`

set datasetname1 = "QCD_Pt_300to470"
set datasetname2 = "QCD_Pt_470to600"
set datasetname3 = "QCD_Pt_600to800"
set datasetname4 = "QCD_Pt_800to1000"
set datasetname5 = "QCD_Pt_1000to1400"
set datasetname6 = "QCD_Pt_1400to1800"
set datasetname7 = "QCD_Pt_1800"

set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/NtupleMaking/01282013/smallerFiles"
#set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/TEMP"
set outDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/01282013/smallerFiles"
set dirlist = "$mainDir/${datasetname1} $mainDir/${datasetname2} $mainDir/${datasetname3} $mainDir/${datasetname4} $mainDir/${datasetname5} $mainDir/${datasetname6} $mainDir/${datasetname7}" 
set outdirlist = ( "$outDir/${datasetname1}" "$outDir/${datasetname2}" "$outDir/${datasetname3}" "$outDir/${datasetname4}" "$outDir/${datasetname5}" "$outDir/${datasetname6}" "$outDir/${datasetname7}") 
set i = 0

set fileList="none"

foreach dir ( $dirlist )

	echo "dir = $dir" 
	if ( ! -d $dir ) then
		mkdir -vp $dir
		@ i += 1
		cp -v fnalCondorPrep.sh $dir
		cp -v fnalCondorSubmit.csh $dir
		cp -v fnalMakeSplitInputs.pl $dir
		cp -v runLostLeptonTreeMaker_condor.py $dir

		if ( $i == 1 ) then
			cp QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname1"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 1 2 $fileList $i $rootdir 
			#./fnalCondorSubmit.csh 7 100 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 13 50 $fileList $i $rootdir 
			cd $curDir
			
			echo "$PWD"
			#exit
		else if ( $i == 2 ) then
			cp QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname2"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 5 100 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 10 50 $fileList $i $rootdir 
			cd $curDir
		else if ( $i == 3 ) then
			cp QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir

			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname3"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 5 100 $fileList $i $rootdir
			./fnalCondorSubmit.csh 10 50 $fileList $i $rootdir
			cd $curDir
		else if ( $i == 4 ) then
			cp QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname4"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 5 100 $fileList $i $rootdir
			./fnalCondorSubmit.csh 10 50 $fileList $i $rootdir
			cd $curDir
		else if ( $i == 5 ) then
			cp QCD_Pt_1000to1400_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname5"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 3 100 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 5 50 $fileList $i $rootdir 
			cd $curDir

		else if ( $i == 6 ) then
			cp QCD_Pt_1400to1800_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname6"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 3 100 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 5 50 $fileList $i $rootdir 
			cd $curDir

		else if ( $i == 7 ) then
			cp QCD_Pt_1800_TuneZ2star_8TeV_pythia6_Summer12_53X.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname7"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 2 100 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 3 50 $fileList $i $rootdir 
			cd $curDir

		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end
