#! /bin/csh

set curDir = `pwd`

set datasetname1 = "TT_CT10"
set datasetname2 = "WJets_HT-300To400_2"
set datasetname3 = "WJets_HT-400ToInf_2"
set datasetname4 = "ZJets_400_HT_inf_2"
set datasetname5 = "WJets_HT-300To400_kasmi"

#/uscmst1b_scratch/lpc1 is READ-ONLY from March 5th, 2013
#set mainDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/"
#set outDir = "/uscmst1b_scratch/lpc1/3DayLifetime/samantha/NewNtuples/03012013"
set mainDir = "/uscms_data/d3/samantha/NtupleMaking/"
set outDir = "/uscms_data/d3/samantha/NtupleMaking/03042013"
set dirlist = "$mainDir/${datasetname1} $mainDir/${datasetname2} $mainDir/${datasetname3} $mainDir/${datasetname4}" 
set outdirlist = ( "$outDir/${datasetname1}" "$outDir/${datasetname2}" "$outDir/${datasetname3}" "$outDir/${datasetname4}") 
set i = 0

set fileList="none"

foreach dir ( $dirlist )

	echo "dir = $dir" 
	@ i += 1
	if ( ! -d $dir ) then
		mkdir -vp $dir
		cp -v fnalCondorPrep.sh $dir
		cp -v fnalCondorSubmit.csh $dir
		cp -v fnalMakeSplitInputs.pl $dir
		cp -v runLostLeptonTreeMaker_condor.py $dir
		cp -v DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root $dir

		if ( $i == 1 ) then
			#cp TT_CT10_TuneZ2star_8TeV-powheg-tauola.list $dir
			#cd $dir
			#set fileList =  `ls -1 *.list`
			#echo "$PWD -> $fileList"
			#set rootdir = "$outDir/$datasetname1"
			#if ( ! -d $rootdir ) then
			#	mkdir -vp $rootdir
			#endif
			#./fnalCondorSubmit.csh 52 100 $fileList $i $rootdir
			#cd $curDir
		else if ( $i == 2 ) then
			cp WJetsToLNu_HT-300To400_8TeV-madgraph_Summer12.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname2"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 7 30 $fileList $i $rootdir 
			cd $curDir

		else if ( $i == 3 ) then
			cp WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname3"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 67 30 $fileList $i $rootdir 
			cd $curDir

		else if ( $i == 4 ) then
			cp ZJetsToNuNu_400_HT_inf_8TeV_madgraph_Summer12_v1.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname4"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 7 30 $fileList $i $rootdir 
			cd $curDir

		#kasmi's WJet samples
		else if ( $i == 5 ) then
			cp WJetsToLNu_HT-300To400_8TeV-madgraph_kasmi.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			set rootdir = "$outDir/$datasetname5"
			echo "$PWD -> $fileList"
			if ( ! -d $rootdir ) then
				mkdir -vp $rootdir
			endif
			./fnalCondorSubmit.csh 21 100 $fileList $i $rootdir 
			cd $curDir

		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		#exit
	endif 

end
