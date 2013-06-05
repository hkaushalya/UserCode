#! /bin/csh

set systVariation = 0

set curDir = `pwd`
set mainDir = "$curDir/05222013_UseDefaultJetRes_100sampling"

if (! -d $mainDir ) then
	mkdir $mainDir
else
	echo "$mainDir  already exists!"
	exit
endif




#set ttdir           = "${mainDir}/TTC10"

#foreach i (1 5 13 29 61 63) #these are the bitmasks enabling each cut cumulatively
foreach i (0 383 255) #these are the bitmasks enabling each cut cumulatively

			set ttdir = "${mainDir}/$i"
			mkdir -vp $ttdir
			cp -v fnalCondorPrep.sh $ttdir
			cp -v fnalCondorSubmit.csh $ttdir
			cp -v fnalMakeSplitInputs.pl $ttdir
			cp -v ../runsmear $ttdir
			cp /share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/input_files/TT_CT10_wFilters_farm_v03012013.files $ttdir
			cd $ttdir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 26 2  $fileList $systVariation $i 
			cd $curDir
end

cd $curDir
