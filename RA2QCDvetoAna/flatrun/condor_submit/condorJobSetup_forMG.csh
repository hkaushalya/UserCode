#! /bin/csh

set systVariation = 0

if ( $systVariation == 0 ) then
	set subdir = "Mean"
else 
	set subdir = "Syst${systVariation}"
endif

set curDir = `pwd`
set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/01222013_MG_5kSampling_newHtBins/${subdir}"
set dirlist = "$mainDir/qcd1 $mainDir/qcd2 $mainDir/qcd3"
set i = 0

foreach dir ( $dirlist )

	if ( ! -d $dir ) then
		mkdir -vp $dir
		@ i += 1
		cp -v fnalCondorPrep.sh $dir
		cp -v fnalCondorSubmit.csh $dir
		cp -v fnalMakeSplitInputs.pl $dir
		cp -v ../runsmear $dir
		if ( $i == 1 ) then
			cp ../input_files/QCD_HT_250To500_MGPythia_v1_lpc1_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 4 4 $fileList ${i} 
			#./fnalCondorSubmit.csh 10 8 $fileList $systVariation 
			./fnalCondorSubmit.csh 15 5 $fileList $systVariation 

			cd $curDir
			#exit
		else if ( $i == 2 ) then
			cp ../input_files/QCD_HT_500To1000_MGPythia_v1_lpc1_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 9 9 $fileList $systVariation 
			./fnalCondorSubmit.csh 15 6 $fileList $systVariation 
			cd $curDir
		else if ( $i == 3 ) then
			cp ../input_files/QCD_HT_1000ToInf_MGPythia_v1_lpc1_wFilters_farm.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 8 4  $fileList $systVariation 
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end
