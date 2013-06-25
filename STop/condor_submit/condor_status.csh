#!/bin/csh
#############################################################
#new utils to do manipulate condor jobs
#by looking into log files.
#This is needed to handle jobs run with very small
#files as they exceed the input limits to most unix and ROOT
#commands.
#Sam June 2013
#############################################################
set ext = "log"
set logs = `ls -1 ${PWD}/*.${ext}`
if ( $status != 0 ) then
	echo "No log files with extension $ext found in $PWD"
	exit
endif

set holdIdle = 0
set releaseHold = 0
set resubmitfailed = 0
set nJobsInIdle = 0
set nJobsRunning = 0
set nJobsResubmitted = 0
set nJobsOnHold = 0
set nJobsPutOnHold = 0
set nJobsDoneWithSuccess = 0

foreach log ( $logs )
	set stat = `grep "return value" $log | cut -d " " --fields=2,3`
	set exitstat = `grep "return value" $log | cut -d " " --fields=6 | cut -d ")" -f1`
	set base = `basename $log`
	echo -n "${base}:"
	set jid = `echo $base | cut -d "_" -f3`
	if ( "${stat}" == "Normal termination" ) then
			if ( "${exitstat}" != "0" ) then
				echo -n "${stat} exit status $exitstat"
				echo -n " !=0. Resubmit job."
				if ( $resubmitfailed == 1 ) then
					#remove old log files
					#ls -l *${jid}*
					rm -vf *${jid}*
					set runfilenum = `echo $base | cut -d "_" -f2`
					set runfile = "run_${runfilenum}.sh"
					#echo "\tSubmit file $runfile"
					#ls -l $runfile
					if (-e $runfile ) then
						condor_submit $runfile
						@ nJobsResubmitted++
					else
						echo -n " Not resubmitted!"
						if (-e $runfile) then
							echo -n "Initial script $runfile not found!" 
						endif
					endif
				endif
				#put a newline 
				echo ""
			else 
				echo "${stat} exit status $exitstat"
				@ nJobsDoneWithSuccess++
			endif
	else 
		set cstat = `condor_q $jid | grep $jid`
		if ( $status == 0) then
			echo -n "Job#$jid is in the farm"
			set runstat = `condor_q -run $jid | grep $jid`
			if ( $status == 0 ) then
				echo " RUNNING"
				@ nJobsRunning++
			else 
				set holdstat = `condor_q -hold $jid | grep $jid`
				if ( $status == 0 ) then
					echo " HOLD"
					@ nJobsOnHold++
				else 
					echo " IDLE"
					@ nJobsInIdle++
					if ($holdIdle == 1) then
						echo -n "\t putting on hold as requested. New status "
						condor_hold $jid
						set holdstat = `condor_q -hold $jid | grep $jid`
						if ( $status == 0 ) then
							echo " HOLD"
							@ nJobsPutOnHold++
						else 
							echo " Attempt to HOLD #$jid UNSUCCESSFUL"
						endif
					endif
				endif
			endif

		else 
			echo "Job#$jid is not in the farm: $cstat"
		endif
	endif
end

echo "------ Summary --- "
echo "Procesed logfile# $#logs"
echo "nJobsDoneWithSuccess = $nJobsDoneWithSuccess"
echo "nJobsRunning = $nJobsRunning"
echo "nJobsOnHold  = $nJobsOnHold"
echo "nJobsInIdle  = $nJobsInIdle"
echo "nJobsOnHold  = $nJobsOnHold"
if ($holdIdle == 1) then
	echo "nJobsPutOnHold = $nJobsPutOnHold"
endif
if ($resubmitfailed == 1) then
	echo "nJobsResubmitted = $nJobsResubmitted"
endif
