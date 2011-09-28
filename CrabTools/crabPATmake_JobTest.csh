#! /bin/tcsh


if ($#argv < 1) then
	echo "|===================================================|"
	echo "| Usage: $0 dir1 dir2"
	echo "| Will check CRAB PAT job log files for completeness."
	echo "|===================================================|"
	exit
endif

#echo "given dir $1 $2"

########## COLOR OUTPUT STUFF
set LF="\n"
set CR="\r"
set INVT="\033[7m"
set NORM="\033[0m"
set BOLD="\033[1m"
set BLINK="\033[5m"
#UNDR="\033[2m\033[4m"; EOL="\033[0K"; EOD="\033[0J"
set UNDR="\033[4m"
set EOL="\033[0K"
set EOD="\033[0J"
set SOD="\033[1;1f"
set CUR_UP="\033[1A"
set CUR_DN="\033[1B"
set CUR_LEFT="\033[1D"
set CUR_RIGHT="\033[1C"
#-- ANSI code
set SCR_HOME="\033[0;0H" #-- Home of the display

set BLACK_F="\033[30m"
set BLACK_B="\033[40m"
set RED_F="\033[31m"
set RED_B="\033[41m"
set GREEN_F="\033[32m"
set GREEN_B="\033[42m"
set YELLOW_F="\033[33m"
set YELLOW_B="\033[43m"
set BLUE_F="\033[34m"
set BLUE_B="\033[44m"
set MAGENTA_F="\033[35m"
set MAGENTA_B="\033[45m"
set CYAN_F="\033[36m"
set CYAN_B="\033[46m"
set WHITE_F="\033[37m"
set WHITE_B="\033[47m"
###### END COLOR OUTPUT STUFF

foreach file ($argv)
	set RESDIR="$file/res/"
	set LOGDIR="$file/log/"
	set SHAREDIR="$file/share"
	set crabcfgFile = "${SHAREDIR}/crab.cfg"
	set DATASET = `grep "datasetpath" ${crabcfgFile} | cut -d "=" -f2`	
	set dbsTempFile = "dbs.tmp"
	dbs listFiles --path=${DATASET} > ${dbsTempFile}
	set DBS_FILELIST = `cat ${dbsTempFile} | grep "root"`
	set DBS_NFILES = `cat ${dbsTempFile} | grep "Total files" | cut -d " " -f4`
	set DBS_NEVENTS = `dbsql "find sum(block.numevents), dataset where dataset like ${DATASET}" | grep ${DATASET} | cut -d " " -f1`

	if (-d $RESDIR) then
		
		echo "\n${RED_F}>>>>>>>>> Processing ${BLUE_F}${RESDIR} ${NORM}"
		echo "DATASET (cfg) = ${DATASET}" 
		echo "DATASET files/events = ${DBS_NFILES} / ${DBS_NEVENTS}" 
		echo "\nChecking CRAB JOB STATUS:"
		set tempFile = "tmp.txt"
		#crab -status -c ${file} > ${tempFile}
		#crab -status -c ${file} | grep 
		#grep "ExitCodes"

#ls -1 ${RESDIR}
		echo "=========== N output files "
		set NCRABJOBS_CREATED = `grep -m1 "jobs created to run " ${LOGDIR}/crab.log | cut -d " " -f 4`
		echo "CRAB Creacted Jobs = ${NCRABJOBS_CREATED}"
		echo -n "*.stderr = "
		ls -1 ${RESDIR}/*.stderr | wc -l
		echo -n "*.stdout = "
		ls -1 ${RESDIR}/*.stdout | wc -l
		echo "=========== ${BOLD}Any errors ${NORM}"
		cat ${RESDIR}/*.stderr
		echo "=========== ${BOLD}Opened/Closed data files ${NORM}"
		echo "Input files From DBS = ${DBS_NFILES}"
		set openedDataFiles = `grep "Successfully opened file" ${RESDIR}/*.stdout | wc -l`
		set closedDataFiles = `grep "Closed file" ${RESDIR}/*.stdout | wc -l`
		echo "Opened / Closed = ${openedDataFiles} / ${closedDataFiles}"
		if ( ${DBS_NFILES} != ${openedDataFiles} ) then
			echo "${RED_F} FILES OPENED DO NOT MATCHED WITH DBS! ${NORM}"
		else 
			echo "${GREEN_F} Files Opened Matched With DBS! ${NORM}"
			if ( ${DBS_NFILES} != ${closedDataFiles} ) then
				echo "${RED_F} FILES CLOSED DO NOT MATCHED WITH DBS! ${NORM}"
			else
				echo "${GREEN_F} Files Closed Matched With DBS! ${NORM}"
			endif
		endif
		echo "=========== ${BOLD}Output Copy Successes ${NORM}"
#grep "Copy succedeed with srm-lcg utils" ${RESDIR}/*.stdout | grep "echo" | wc -l
		set copySuccessFiles = `grep "Copy succedeed with srm-lcg utils" ${RESDIR}/*.stdout | grep "echo" | wc -l`
		echo "Expected # of output root files = ${MAGENTA_F}${NCRABJOBS_CREATED}${NORM}"
		if ( ${copySuccessFiles} != ${NCRABJOBS_CREATED} ) then
			echo "${RED_F} # of  output files successfully copied DO NOT MATCH with created CRAB Jobs! ${NORM}"
		else
			echo "${GREEN_F} # of output files successfully copied match with created CRAB Jobs. ${NORM}"
		endif
#echo "=========== Output Copy Fail"
#grep -i "fail"  ${RESDIR}/*.stdout

		echo "=========== ${BOLD}Trigger Report ${NORM}"
		set totEvents = `grep "TrigReport Events total" ${RESDIR}/*.stdout | awk '{f+=$5}END{print f}'`
		echo "Total Processed Events  = ${BOLD}${CYAN_F}${totEvents}${NORM}"
		if ( ${totEvents} != ${DBS_NEVENTS} ) then
			echo "${RED_F} EVENTS PROCESSED DO NOT MATCHED WITH DBS! ${NORM}"
		else
			echo "${GREEN_F} EVENTS PROCESSED MATCHED WITH DBS! ${NORM}"
		endif
		echo "=========== Output file copy location"
		set OUTLOCATION = `cat ${RESDIR}/*.stdout | grep -m 1 "newLfn" | cut -d " " -f8`

		set OUTFILE = `echo ${OUTLOCATION:t}`
		set OUTDIR = `echo ${OUTLOCATION:h}`

		set ABSDIR = "/pnfs/cms/WAX/11${OUTDIR}"

		echo "OUT DIR = ${BLUE_F} ${OUTDIR} ${NORM}"
		echo "ABS PATH= ${BLUE_F} ${ABSDIR} ${NORM}"

		set filesFound = `ls -1 ${ABSDIR}/*.root |wc -l`
		echo "Root files found in OUT DIR = ${MAGENTA_F}${filesFound}${NORM}"

		if ( ${copySuccessFiles} == ${filesFound} ) then
			echo "${GREEN_F}Output root files FOUND MATCH WITH EXPCTED. ${NORM}"
		else 
			echo "${RED_F}${BLINK}Output root files FOUND DO NOT MATCH WITH EXPCTED. ${NORM}"
			if ( ${filesFound} > ${copySuccessFiles} ) then
				echo "Possible duplicate output files found! Executing crabPAT_listDupRootFiles.csh"
				echo "Carefully answer YES [y] or NO [n] to remove duplciates"
				echo "Files shown are time ordered."
				./crabPAT_listDupRootFiles.csh ${ABSDIR}
			endif
		endif
		echo "${RED_F}<<<<<<<<< END Report for ${BLUE_F}${RESDIR} ${NORM}"
			
		rm -rf ${tempFile}
	else
		echo "Copying $RESDIR"
		echo "Skipping $RESDIR (is not a directory)"
	endif

	rm -rf ${dbsTempFile}
end


