#! /bin/csh
#rename the resubmitedd crab jobs output files

set curDir = `pwd`
#set dir = "/uscms_data/d3/samantha/SUSYRA2work/CMSSW_4_2_5/src/UserCode/RA2QCDvetoAna/test/FactorizationResults/11192011_Data/Data4/res"
set dir = "/uscms_data/d3/samantha/METRelValTesting_new/CMSSW_4_4_0/src/Validation/RecoMET/test/HiPU/ttbar/res/"
set jobs = (3 8 16 17 20 22 26 31 34 36 40 47 50 51 60 66 67 68 71 73 77 84 85 98 99 102 104)
echo "Total jobs to be moved =  $#jobs from directory"
echo "$dir"

foreach job ( $jobs )
	echo "job = $job"
	set files = `ls -1 | grep "_${job}\."`
	foreach file ( $files )
		set newname="${file}_bad"
		echo "$file -> $newname"
		mv -v $file $newname
	end
end
