#! /bin/csh

 
foreach systVariation (0 1 10 2 20 3 4)

if ( $systVariation == 0 ) then
	set subdir = "Mean"
else 
	set subdir = "Syst${systVariation}"
endif

set curDir = `pwd`
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/12112012_chsJets_AllNjetHtMhtBins_10kSamling/"
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/12112012_chsJets_AllNjetHtMhtBins_5kSamling/RatioFineBinned_Syst${systVariation}/"
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/01102013_5kSamling_NoLumiWgt_2010SystVariation/${subdir}"
set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/01292013_test2/${subdir}"
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/01192013_MC_QCDforBGoverlays/${subdir}"
set dirlist = "$mainDir/qcd1 $mainDir/qcd2 $mainDir/qcd3 $mainDir/qcd4 $mainDir/qcd5 $mainDir/qcd6 $mainDir/qcd7"
set i = 0

set fileList="none"

foreach dir ( $dirlist )

	if ( ! -d $dir ) then
		mkdir -vp $dir
		@ i += 1
		cp -v fnalCondorPrep.sh $dir
		cp -v fnalCondorSubmit.csh $dir
		cp -v fnalMakeSplitInputs.pl $dir
		cp -v ../runsmear $dir
		if ( $i == 1 ) then
			cp ../input_files/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			echo "condorJobSetup:: systVariation = ${systVariation}"
			#./fnalCondorSubmit.csh 4 4 $fileList $systVariation 
			./fnalCondorSubmit.csh 15 3 $fileList $systVariation 
			cd $curDir
			#exit
		else if ( $i == 2 ) then
			cp ../input_files/QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 9 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 3 ) then
			cp ../input_files/QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 9 3  $fileList $systVariation 
			cd $curDir
		else if ( $i == 4 ) then
			cp ../input_files/QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 10 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 5 ) then
			cp ../input_files/QCD_Pt_1000to1400_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 12 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 6 ) then
			cp ../input_files/QCD_Pt_1400to1800_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 12 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 7 ) then
			cp ../input_files/QCD_Pt_1800_TuneZ2star_8TeV_pythia6_Summer12_wFilters_farm_v01282013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 6 3 $fileList $systVariation 
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end

end
