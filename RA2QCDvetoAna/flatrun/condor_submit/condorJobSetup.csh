#! /bin/csh

set systVariation = 0
set 
#if ( $systVariation == 0 )

set curDir = `pwd`
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/12112012_chsJets_AllNjetHtMhtBins_10kSamling/"
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/12112012_chsJets_AllNjetHtMhtBins_5kSamling/RatioFineBinned_Syst${systVariation}/"
set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/12112012_chsJets_AllNjetHtMhtBins_5kSamling/RatioFineBinned/"
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
			cp ../QCD1_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 4 4 $fileList ${i} 
			./fnalCondorSubmit.csh 15 3 $fileList $systVariation 
			cd $curDir
			#exit
		else if ( $i == 2 ) then
			cp ../QCD2_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 9 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 3 ) then
			cp ../QCD3_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 9 3  $fileList $systVariation 
			cd $curDir
		else if ( $i == 4 ) then
			cp ../QCD4_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 10 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 5 ) then
			cp ../QCD5_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 12 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 6 ) then
			cp ../QCD6_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 12 3 $fileList $systVariation 
			cd $curDir
		else if ( $i == 7 ) then
			cp ../QCD7_wFilters_farm.files $dir
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
