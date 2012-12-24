#! /bin/csh
###########################################
# For testing Flat Ntuple making.
# Only for jobs run manually on the
# cmslpcxx nodes.
# Not designed for jobs run on condor.
###########################################

if ( $#argv != 2 ) then
	echo "Usage: $0 job_directory output_directory"
	exit
endif

foreach dir ( $argv )
	if ( ! -d $dir || ! -e $dir ) then
		echo "\033[1;31m $dir is not a directory or does not exist.\033[0m"
		exit
	endif 
end

set curDir = `pwd`
set jobDir = $1
set outDir = $2

echo "Test 1: Basic file count ..."
set nLogFiles = `ls -l ${jobDir}/*.log | wc -l`
#echo "log files = $nLogFiles"
set nRootFiles = `ls -l ${outDir}/*.root | wc -l`
echo -n "root/log files = $nRootFiles/$nLogFiles"
if ( $nLogFiles != $nRootFiles ) then
	echo "\033[1;31m WARN SRC tobe processed does not match to the given! \033[0m"
	exit
else 
	echo "\033[1;32m Success. \033[0m"
endif

set fileList = `ls -1 ${jobDir}/inputList_*.rtlst`
if ( $#fileList == 0 ) then
	echo "File list/s not found!"
	exit
else 
	echo "File lists found"
endif

set logFileExt = 'log'
set filesToProcess = `cat ${jobDir}/*.rtlst | wc -l`
set SRCopened = `grep "Successfully opened" ${jobDir}/*.$logFileExt | wc -l`
set SRCclosed = `grep "Closed file" ${jobDir}/*.$logFileExt | wc -l`
set rootFiles = `ls -1 ${outDir} | wc -l`
set exceptions = `cat ${jobDir}/*.$logFileExt| grep -i "Exception" | wc -l`

echo "Total SRC files to be processed = $filesToProcess"
echo "Total SRC files opened          = $SRCopened"
echo "Total SRC files closed          = $SRCclosed"
echo "Total ROOT files                = $rootFiles"
echo "Total Exceptions occured        = $exceptions"

if ( $filesToProcess != $SRCopened ) then
	echo "\033[1;31m WARNING>> SRC tobe processed does not match opened\!\!\! \033[0m"
else 
	echo "\033[1;32m >>SRC tobe processed match opened. \033[0m"
endif

if ( $filesToProcess != $SRCclosed ) then
	echo "\033[1;31m WARNNING>> SRC tobe processed does not match closed\!\!\! \033[0m"
else 
	echo "\033[1;32m >>SRC tobe processed match closed. \033[0m"
endif

#now figure out which jobs failed to process all the inputs files, failed to
#produce output file, or had exceptions and quit prematurely.

if ( $filesToProcess != $SRCopened || $filesToProcess != $SRCclosed || $exceptions != 0) then
	echo "Processed files mismatch with inputs..."
	echo "Jobs that has not processed all files! * indicate missing output root file."
	foreach file ( $fileList )
		set num = `ls ${jobDir}/${file} | cut -d "_" -f2 | cut -d "." -f1`
		set logFile = "${jobDir}/job_${num}.${logFileExt}"
		set rootFile = "${outDir}/Ntuple_${num}.root"
		set nFiles2Proc = `cat ${jobDir}/inputList_${num}.rtlst | wc -l`
		set nOpened = `grep "Successfully opened" ${jobDir}/${logFile} | wc -l`
		set nClosed = `grep "Closed file" ${jobDir}/${logFile} | wc -l`
		
		if ( $nFiles2Proc != $nOpened || $nFiles2Proc != $nClosed) then
			echo -n " $num"
				if ( ! -e $rootFile ) then
					echo -n "*"
				endif
		endif
	end
	echo "\n"
endif


