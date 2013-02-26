#! /bin/csh

set curDir = `pwd`
set datasetname1 = "JetHT_Run2012B_13Jul2012_v1_lpc1"
set datasetname2 = "JetHT_Run2012C_24Aug2012_v2_lpc1"
set datasetname3 = "JetHT_Run2012C_PromptReco_v2_lpc1"
set datasetname4 = "JetHT_Run2012D_PromptReco_v1_lpc1"


set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/NtupleMaking/01282013"
#set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/TEMP"
set outDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/01282013"
#set dirlist = "$mainDir/${datasetname1} $mainDir/${datasetname2} $mainDir/${datasetname3} $mainDir/${datasetname4}"
#set outdirlist = ( "$outDir/${datasetname1}" "$outDir/${datasetname2}" "$outDir/${datasetname3}"  "$outDir/${datasetname4}")

set dirlist = "$mainDir/${datasetname4}"
set outdirlist = ("$outDir/${datasetname4}")
set i = 3

set fileList="none"

foreach dir ( $dirlist )

	echo "dir = $dir" 
	if ( ! -d $dir ) then
		
		@ i += 1
		if ( $i == 4 ) then
			#continue
			mkdir -vp $dir
			cp -v fnalCondorPrep.sh $dir
			cp -v fnalCondorSubmit.csh $dir
			cp -v fnalMakeSplitInputs.pl $dir
			cp -v runLostLeptonTreeMaker_condor.py $dir
		else 
			continue
		endif

		if ( $i == 1 ) then
			cp JetHT_Run2012B_13Jul2012_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname1"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 1 2 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 13 100 $fileList $i $rootdir 
			cd $curDir
			
			echo "$PWD"
			#exit
		else if ( $i == 2 ) then
			cp JetHT_Run2012C_24Aug2012_v2_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname2"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 3 50 $fileList $i $rootdir 
			cd $curDir
		else if ( $i == 3 ) then
			cp JetHT_Run2012C_PromptReco_v2_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname3"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 17 100 $fileList $i $rootdir
			cd $curDir

		else if ( $i == 4 ) then
			cp JetHT_Run2012D_PromptReco_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname4"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 19 100 $fileList $i $rootdir
			./fnalCondorSubmit.csh 34 55 $fileList 5 $rootdir
			cd $curDir

		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		#exit
	endif 

end
