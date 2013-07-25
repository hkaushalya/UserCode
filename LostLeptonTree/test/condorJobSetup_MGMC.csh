#! /bin/csh

set curDir = `pwd`

set datasetname1 = "MG_QCD_HT_250To500"
set datasetname2 = "MG_QCD_HT_500To1000"
set datasetname3 = "MG_QCD_HT_1000ToInf"

#/uscmst1b_scratch/lpc1 is READ-ONLY from March 5th, 2013
#set mainDir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_5_3_5/src/UserCode/LostLeptonTree/test/NtupleMaking/01282013"
#set outDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/01282013"
set mainDir = "/uscms_data/d3/samantha/NtupleMaking/"
set outDir = "/uscms_data/d3/samantha/NtupleMaking/20130711_PatJetsWithMetPhiCorr"

#set mainDir = "/eos/uscms/store/user/samantha/NtupleMaking"
#set outDir = "/eos/uscms/store/user/samantha/NtupleMaking/03102013"

set dirlist = "$mainDir/${datasetname1} $mainDir/${datasetname2} $mainDir/${datasetname3}" 
set outdirlist = ( "$outDir/${datasetname1}" "$outDir/${datasetname2}" "$outDir/${datasetname3}" ) 
set i = 0

set fileList="none"

foreach dir ( $dirlist )

	echo "dir = $dir" 
	if ( ! -d $dir ) then
		mkdir -vp $dir
		@ i += 1
		cp -r -v fnalCondorPrep.sh $dir
		cp -r -v fnalCondorSubmit.csh $dir
		cp -r -v fnalMakeSplitInputs.pl $dir
#		cp -r -v runLostLeptonTreeMaker_condor.py $dir
		cp -r -v runLostLeptonTreeMaker_pfjets_condor.py $dir
		cp -r -v DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root $dir

		if ( $i == 1 ) then
			cp -r QCD_HT_250To500_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname1"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			#./fnalCondorSubmit.csh 1 2 $fileList $i $rootdir 
			./fnalCondorSubmit.csh 170 30 $fileList $i $rootdir 
			cd $curDir
			exit
		else if ( $i == 2 ) then
			cp -r QCD_HT_500To1000_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname2"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 134 30 $fileList $i $rootdir 
			cd $curDir
		else if ( $i == 3 ) then
			cp -r QCD_HT_1000ToInf_MGPythia_v1_lpc1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			set rootdir = "$outDir/$datasetname3"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 74 30 $fileList $i $rootdir
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end
