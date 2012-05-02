#! /bin/csh
#to check condor job output
set curDir = `pwd`

foreach dir ( $argv )

	if ( -d $dir ) then
		cd -v $dir
		set nLogFiles = `ls -l *.log | wc -l`
		echo "log files = $nLogFiles"
		set nRootFiles = `ls -l *.root | wc -l`
		echo "root/log files = $nRootFiles/$nLogFiles"
		if ( $nLogFiles != $nRootFiles ) then
			echo "\033[1;31m WARN SRC tobe processed does not match to the given! \033[0m"

			echo -n "Do you want to continue? (Y/N) "
			set input = $<
			switch ($input)
				#case [yY][eE][sS]:
				case [yY]:
					echo "argv[1] is 'yes'"
					breaksw
				#case [nN][oO]:
				case [nN]:
					echo "argv[1] is 'no'"
					exit
				default:
					#echo "invalid arg: must be yes or no"
					echo "invalid arg: must be y/Y or n/N"
					exit
			endsw

		endif

		set fileList = `ls *.list`
		if ( ! -e $fileList ) then
			echo "File list not found!"
			exit
		endif
		set filesToProcess = `cat $fileList | wc -l`
		echo "SRC files to be processed = $filesToProcess ($fileList)"
		set SRCtobeProc =  `grep "root" *.stdout| grep "dcap" | wc -l`
		echo "SRC files given     = $SRCtobeProc"
		#grep "root" *.stdout| grep "dcap" | wc -l
		set SRCprocessed = `grep "Successfully opened" *.stderr | wc -l`
		echo "Processed SRC files = $SRCprocessed"
		 
		if ( $filesToProcess != $SRCtobeProc ) then
			echo "\033[1;31m WARNING>> SRC tobe processed does not match to the given\!\!\! \033[0m"
		else 
			echo "\033[1;32m >>SRC tobe processed match to the given. \033[0m"
		endif

		if ( $filesToProcess != $SRCprocessed ) then
			echo "\033[1;31m WARNNING>> SRC tobe processed does not match to the processed\!\!\! \033[0m"
			#figure out which jobs failed
			echo "Jobs that indicates possoble errors are : "
			grep -li "FileOpenError" *.stderr
			grep -li "Error" *.stderr
			exit
		else 
			echo "\033[1;32m >>SRC tobe processed match to the processed. \033[0m"
		endif

		#echo -n "Total event processed = " 
		#grep "TrigReport Events total" *.stderr | awk '{f+=$5}END{print f}'

		echo "Merging any Available root files."
		
		if ( -e "Merged.root" ) then
			echo "\033[1;31m File name Merged.root already exitst. Not merging! \033[0m"
		else 
			hadd Merged.root [0-9]*.root
#			hadd Merged_StdRA2FlatTree.root  StdRA2FlatTree*.root	
		endif

		cd $curDir

	else
		echo "\033[1;31m $dir is not a directory. Skipping.. \033[0m"
	endif 

end
