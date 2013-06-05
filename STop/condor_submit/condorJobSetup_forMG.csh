#! /bin/csh

foreach systVariation (0 1 2 3 4 5 6)
#foreach systVariation (0)

if ( $systVariation == 0 ) then
	set subdir = "Mean"
else 
	set subdir = "Syst${systVariation}"
endif

set curDir = `pwd`
set mainDir = "${curDir}/05222013_MG_100/${subdir}"
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
			cp ../input_files/MG_QCD_HT_250To500_wFilters_farm_v03012013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 4 4 $fileList ${i} 
			./fnalCondorSubmit.csh 84 2 $fileList $systVariation 

			cd $curDir
			#exit
		else if ( $i == 2 ) then
			cp ../input_files/MG_QCD_HT_500To1000_wFilters_farm_v03012013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 67 2 $fileList $systVariation 
			cd $curDir
		else if ( $i == 3 ) then
			cp ../input_files/MG_QCD_HT_1000ToInf_wFilters_farm_v03012013.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 37 2  $fileList $systVariation 
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end

end
