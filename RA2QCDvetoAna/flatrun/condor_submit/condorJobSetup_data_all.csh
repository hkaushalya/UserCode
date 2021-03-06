#! /bin/csh

set curDir = `pwd`
#set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/02132013_dataAll_yields_HT350TrigpresCheck/"
set mainDir = "/share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit/03132013_newJECs_dataAll_noTrigPrescale"
set dirlist = "$mainDir"
set i = 0
set systVariation = 0 
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
			#cp ../input_files/HTMHT_AND_JetHT_ABCD_Periods_local_v01282013.files $dir
			cp ../input_files/HTMHT_AND_JetHT_ABCD_Periods_farm_v03012013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 4 4 $fileList $systVariation 
			#./fnalCondorSubmit.csh 21 7 $fileList $systVariation 
			./fnalCondorSubmit.csh 40 10 $fileList $systVariation 
			cd $curDir
		endif
		sleep 3
	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end

end
