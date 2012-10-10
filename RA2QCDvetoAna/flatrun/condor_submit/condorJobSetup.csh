#! /bin/csh

set curDir = `pwd`
set mainDir = "/share/store/users/samantha/CMSSW_DEV/525/FlatSmearingCode/optimize/test2/submit"
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
			cp ../qcd1.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 20 45 $fileList $i 
			cd $curDir
			echo "$PWD"
		else if ( $i == 2 ) then
			cp ../qcd2.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 20 26 $fileList $i 
			cd $curDir
		else if ( $i == 3 ) then
			cp ../qcd3.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 20 27  $fileList $i 
			cd $curDir
		else if ( $i == 4 ) then
			cp ../qcd4.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 25 22 $fileList $i 
			cd $curDir
		else if ( $i == 5 ) then
			cp ../qcd5.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 10 28 $fileList $i 
			cd $curDir
		else if ( $i == 6 ) then
			cp ../qcd6.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 10 28 $fileList $i 
			cd $curDir
		else if ( $i == 7 ) then
			cp ../qcd7.files $dir
			cd $dir
			set fileList =  `ls -1 *.files`
			echo "$PWD -> $fileList"
			./fnalCondorSubmit.csh 10 14 $fileList $i 
			cd $curDir
		endif

		sleep 3

	else
		echo "A directory named $dir already exists. Skipping.."
		exit
	endif 

end
