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
			exit
		else 

		endif

		#set fileList = `ls *.list`
		set fileList = `ls *.rtlst`
		#if ( ! -e $fileList ) then
		#	echo "File list not found!"
		#	exit
		#endif
		set filesToProcess = `cat $fileList | wc -l`
		echo "SRC files to be processed = $filesToProcess ($fileList)"
		set SRCtobeProc =  `grep "root" *.stdout| grep "dcap" | wc -l`

		echo -n ">>> Testing events processed "
		set outfiles = `ls -1 *.stdout`

		set sucess = 1
		echo $outfiles
		set tot_found = 0;
		set tot_processs = 0
		foreach out ( $outfiles )
			#echo "processing $out"
			set p1 = `cat $out | grep 'Entries found' | cut -d ' ' -f4` 
			set p2 = `cat $out | grep 'Entries found' | cut -d ' ' -f6` 
			#echo "$p1 / $p2 "
			@ tot_found += $p1
			@ tot_processs += $p2
			if ( $p1 != $p2 ) then
				echo "\033[31m entries found /processed mismatch: $p1 / $p2 in job file $out \033[0m"
				@ sucess = 0
				exit
			endif
		end
		if ( $sucess == 1 ) then
			echo "\033[32m All events present are processed. $tot_found/$tot_processs \033[0m"
		endif
	

		echo "Merging any Available root files."
		
		if ( -e "Merged.root" ) then
			echo "\033[1;31m File name Merged.root already exitst. Not merging! \033[0m"
		else 
			hadd Merged.root [0-9]*.root
		endif

		cd $curDir

	else
		echo "\033[1;31m $dir is not a directory. Skipping.. \033[0m"
	endif 

end
