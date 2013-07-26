#! /bin/csh

set systVariation = 0

set curDir = `pwd`
set mainDir = "$curDir/04232013_stopnjetcuttest"

if (! -d $mainDir ) then
	mkdir $mainDir
else
	echo "$mainDir  already exists!"
	exit
endif

set wjet300to400dir = "${mainDir}/WJets300to400"
set wjet400toInfdir = "${mainDir}/WJets400toInf"
set ttdir           = "${mainDir}/TTC10"
set znndir          = "${mainDir}/Znn"

#foreach i (1 2 3 4)
foreach i (3)

		if ( $i == 1 ) then
			mkdir -vp $wjet300to400dir
			cp -v fnalCondorPrep.sh $wjet300to400dir
			cp -v fnalCondorSubmit.csh $wjet300to400dir
			cp -v fnalMakeSplitInputs.pl $wjet300to400dir
			cp -v ../runsmear $wjet300to400dir
			cp ../input_files/WJets_HT-300To400_wFilters_farm_v03012013.files $wjet300to400dir
			cd $wjet300to400dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh 4 4 $fileList ${i} 
			./fnalCondorSubmit.csh 2 4 $fileList $systVariation 
			cd $curDir
			#exit
		else if ( $i == 2 ) then
			mkdir -vp $wjet400toInfdir
			cp -v fnalCondorPrep.sh $wjet400toInfdir
			cp -v fnalCondorSubmit.csh $wjet400toInfdir
			cp -v fnalMakeSplitInputs.pl $wjet400toInfdir
			cp -v ../runsmear $wjet400toInfdir
			cp ../input_files/WJets_HT-400ToInf_wFilters_farm_v03012013.files $wjet400toInfdir
			cd $wjet400toInfdir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 14 5 $fileList $systVariation 
			cd $curDir
		else if ( $i == 3 ) then
			mkdir -vp $ttdir
			cp -v fnalCondorPrep.sh $ttdir
			cp -v fnalCondorSubmit.csh $ttdir
			cp -v fnalMakeSplitInputs.pl $ttdir
			cp -v ../runsmear $ttdir
			cp ../input_files/TT_CT10_wFilters_farm_v03012013.files $ttdir
			cd $ttdir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 11 5  $fileList $systVariation 
			cd $curDir
		else if ( $i == 4 ) then
			mkdir -vp $znndir
			cp -v fnalCondorPrep.sh $znndir
			cp -v fnalCondorSubmit.csh $znndir
			cp -v fnalMakeSplitInputs.pl $znndir
			cp -v ../runsmear $znndir
			cp ../input_files/ZJets_400_HT_inf_wFilters_farm_v03012013.files $znndir
			cd $znndir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 2 4  $fileList $systVariation 
			cd $curDir
		endif
end

cd $curDir
