#! /bin/csh
#to check condor job out
#uses trigger report to find the total number of events processed.
set curDir = `pwd`

set PYTEVENTS = ( 526533 2077149 3114689 2037669 2157438 1984336 959011 868617 95373 )

foreach dir ( $argv )
	if ( -d $dir ) then
		cd -v $dir
		set total = `grep "TrigReport Events total" *.stderr | awk '{f+=$5}END{print f}'`
		echo "\033[1;32m TrigReport Events total = $total \033[0m"
		#set totalQCDfromSmearing = `grep -m1 "QCDfromSmearing" *.stderr | grep "TrigReport" | awk '{f+=$4}END{print f}'`
		#echo "\033[1;32m TrigReport QCDfromSmearing = $totalQCDfromSmearing \033[0m"
		set totalQCDfromSmearing = `grep -m1 "factorization_ht500" *.stderr | grep "TrigReport" | awk '{f+=$4}END{print f}'`
		echo "\033[1;32m TrigReport factorization_ht500 = $totalQCDfromSmearing \033[0m"
		cd $curDir

	else
		echo "\033[1;31m $dir is not a directory. Skipping.. \033[0m"
	endif 
end

echo "PYTHIA total events in record: "
echo $PYTEVENTS
