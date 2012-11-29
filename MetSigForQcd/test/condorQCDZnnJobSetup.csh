#! /bin/csh

set curDir  = `pwd`
set mainDir = "Results/CondorJobs/3JetsIncl_DphiMinAfterMScut"
set dirlist = "$mainDir/QCD/qcd1 $mainDir/QCD/qcd2 $mainDir/QCD/qcd3 $mainDir/QCD/qcd4 $mainDir/QCD/qcd5 $mainDir/QCD/qcd6 $mainDir/QCD/qcd7 $mainDir/QCD/qcd8 $mainDir/QCD/qcd9"
set i = 0
#set dirlist = "$mainDir/qcd1 $mainDir/qcd2 $mainDir/qcd3 $mainDir/qcd4"
#set i = 9 #to run only madgraph samples

if ( -d $mainDir ) then
	echo "Output base dir ${mainDir} already exists! proceed"
	echo -n "Do you want to continue? (Y/N) "
	set input = $<
		
	switch ($input)
		#case [yY][eE][sS]:
		case [yY]:
			echo "argv[1] is 'yes'"
			breaksw
		#case [nN][oO]:
		case [nN]:
			echo "argv[1] is 'no'"
			exit
		default:
		#echo "invalid arg: must be yes or no"
			echo "invalid arg: must be y/Y or n/N"
			exit
	endsw
endif

set fileList="none"

foreach dir ( $dirlist )

	if ( ! -d $dir ) then
		mkdir -vp $dir
		@ i += 1
		cp -v fnalCondorPrep.sh $dir
		cp -v fnalCondorSubmit.csh $dir
		cp -v fnalMakeSplitInputs.pl $dir
		cp -v runMetSig_MCcondor.py $dir

		if ( $i == 1 ) then
			#cp QCD_Pt_120to170_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp -v 19September2011_QCD_Pt_120to170_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 40 20 $fileList $i 
			./fnalCondorSubmit.csh 6 50 $fileList $i 
			cd $curDir
			echo "$PWD"
		else if ( $i == 2 ) then
			#cp QCD_Pt_170to300_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp -v 19September2011_QCD_Pt_170to300_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 14 30 $fileList $i 
			./fnalCondorSubmit.csh 6 50 $fileList $i 
			cd $curDir
		else if ( $i == 3 ) then
			#cp QCD_Pt_300to470_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_300to470_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 14 30 $fileList $i 
			./fnalCondorSubmit.csh 6 50 $fileList $i 
			cd $curDir
		else if ( $i == 4 ) then
			#cp QCD_Pt_470to600_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_470to600_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 9 30 $fileList $i 
			./fnalCondorSubmit.csh 3 60 $fileList $i 
			cd $curDir
		else if ( $i == 5 ) then
			#cp QCD_Pt_600to800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_600to800_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 10 30 $fileList $i 
			./fnalCondorSubmit.csh 3 60 $fileList $i 
			cd $curDir
		else if ( $i == 6 ) then
			#cp QCD_Pt_800to1000_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_800to1000_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 10 30 $fileList $i 
			./fnalCondorSubmit.csh 3 60 $fileList $i 
			cd $curDir
		else if ( $i == 7 ) then
			#cp QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_1000to1400_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 6 30 $fileList $i 
			./fnalCondorSubmit.csh 2 55 $fileList $i 
			cd $curDir
		else if ( $i == 8 ) then
			#cp QCD_Pt_1400to1800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_1400to1800_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 5 30 $fileList $i 
			./fnalCondorSubmit.csh 2 60 $fileList $i 
			cd $curDir
		else if ( $i == 9 ) then
			#cp QCD_Pt_1800_TuneZ2_7TeV_pythia6_Summer11_PU_S3_START42_V11_v2_AODSIM.list $dir
			cp 19September2011_QCD_Pt_1800_TuneZ2_pythia6_Summer11PUS3.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 2 25 $fileList $i 
			./fnalCondorSubmit.csh 1 31 $fileList $i 
			cd $curDir
		else if ( $i == 10 ) then
			#cp QCD_TuneZ2_HT_100To250_7TeV_madgraph_Summer11_PU_S4_START42_V11_v1_AODSIM.list $dir
			cp QCD_TuneZ2_HT-100To250_7TeV-madgraph_anwar.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 50 22 $fileList $i 
			./fnalCondorSubmit.csh 48 100 $fileList $i 
			cd $curDir
		else if ( $i == 11 ) then
			#cp QCD_TuneZ2_HT_250To500_7TeV_madgraph_Summer11_PU_S4_START42_V11_v3_AODSIM.list $dir
			cp QCD_TuneZ2_HT-250To500_7TeV-madgraph_anwar.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 50 23 $fileList $i 
			./fnalCondorSubmit.csh 42 100 $fileList $i 
			cd $curDir
		else if ( $i == 12 ) then
			#cp QCD_TuneZ2_HT_500To1000_7TeV_madgraph_Summer11_PU_S4_START42_V11_v1_AODSIM.list $dir
			cp QCD_TuneZ2_HT-500To1000_7TeV-madgraph_anwar.list $dir
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 40 21 $fileList $i 
			./fnalCondorSubmit.csh 29 100  $fileList $i 
			cd $curDir
		else if ( $i == 13 ) then
			#cp QCD_TuneZ2_HT_1000_7TeV_madgraph_Summer11_PU_S4_START42_V11_v1_AODSIM.list $dir
			cp QCD_TuneZ2_HT-1000_7TeV-madgraph_anwar.list $dir 
			cd $dir
			set fileList =  `ls -1 *.list`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 10 40 $fileList $i 
			./fnalCondorSubmit.csh 13 100 $fileList $i 
			cd $curDir


		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end

exit

#now submit Znn condor jobs
echo ">>>>>>>> Submitting Znn jobs "

set znnDir = "${mainDir}/Znn"
mkdir -vp $znnDir
cp -v /uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/ZJetsToNuNu_200_HT_inf_7TeV-madgraph_Summer11_PU_S4_START42.list $znnDir
cp -v fnalCondorPrep.sh $znnDir
cp -v fnalCondorSubmit.csh $znnDir
cp -v fnalMakeSplitInputs.pl $znnDir
cp -v runMetSig_MCcondor.py $znnDir
cd $znnDir
./fnalCondorSubmit.csh 8 100 ZJetsToNuNu_200_HT_inf_7TeV-madgraph_Summer11_PU_S4_START42.list 50
cd $curDir

