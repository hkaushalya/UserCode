#! /bin/csh

set curDir = `pwd`

set datasetname1 = "MG_QCD_HT_250To500"
set datasetname2 = "MG_QCD_HT_500To1000"
set datasetname3 = "MG_QCD_HT_1000ToInf"

set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/NtupleMaking/01282013"
set outDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/01282013"
set dirlist = "$mainDir/${datasetname1} $mainDir/${datasetname2} $mainDir/${datasetname3}" 
set outdirlist = ( "$outDir/${datasetname1}" "$outDir/${datasetname2}" "$outDir/${datasetname3}" ) 
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
			cp QCD_HT_250To500_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname1"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 1 2 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 60 85 $fileList $i $rootdir 
			cd $curDir
			#exit
		else if ( $i == 2 ) then
			cp QCD_HT_500To1000_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname2"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 50 82 $fileList $i $rootdir 
			cd $curDir
		else if ( $i == 3 ) then
			cp QCD_HT_1000ToInf_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname3"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 32 70 $fileList $i $rootdir
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end
