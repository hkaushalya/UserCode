#! /bin/csh
#to check condor job output
set curDir = `pwd`

foreach dir ( $argv )

	if ( -d $dir ) then
		cd -v $dir
		set failedJobs = `grep -l "error" *.stderr`
		foreach job ( $failedJobs )
			echo "---------- Failed job# $job"
			set prefix = `echo $job | cut -d "." -f1`
			set jid    = `echo $job | cut -d "_" -f3`
			set jno    = `echo $job | cut -d "_" -f2`

			#echo "prefix/jid/jno = $prefix $jid $jno"
			set rmlist = "$prefix.* $jid.root"

			#echo "Files to remove $rmlist"
			rm -vf $rmlist
			set subnew = "run_$jno.sh"
			echo "Submitting new job $subnew"
			condor_submit $subnew
		end
		cd -v $curDir
	endif 
end
