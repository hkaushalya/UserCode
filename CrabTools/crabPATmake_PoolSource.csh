#! /bin/tcsh


if ($#argv < 1) then
	echo "|===================================================|"
	echo "| Usage: $0 dir1"
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

	if (-d $RESDIR) then
		
		echo "\n${RED_F}>>>>>>>>> Processing ${BLUE_F}${RESDIR} ${NORM}"
		echo "${DATASET}"
		set OUTLOCATION = `cat ${RESDIR}/*.stdout | grep -m 1 "newLfn" | cut -d " " -f8`

		set OUTFILE = `echo ${OUTLOCATION:t}`
		set OUTDIR = `echo ${OUTLOCATION:h}`

		set ABSDIR = "/pnfs/cms/WAX/11${OUTDIR}"

		echo "OUT DIR = ${BLUE_F} ${OUTDIR} ${NORM}"
		#echo "ABS PATH= ${BLUE_F} ${ABSDIR} ${NORM}"

		#set filesFound = `ls -1 ${ABSDIR}/*.root |wc -l`
		#echo "Root files found in OUT DIR = ${MAGENTA_F}${filesFound}${NORM}"
		echo "dataset = ${DATASET}" 
		echo ${DATASET} | sed "s/\//_/g" | sed "s/\//_/" | sed "s/^_//" | sed "s/-/_/g" 
		set py_prefix = `echo ${DATASET} | sed "s/\//_/g" | sed "s/\//_/" | sed "s/^_//" | sed "s/-/_/g"`
		echo "prefix = ${py_prefix}"
		set py_file = "${py_prefix}_cfi.py"
		if ( -e ${py_file} ) then
			echo "File ${py_file} already exists! exiting!"
			exit
		else 
			echo "Creating python source file ${py_file}"
		endif
	
		########################
		# Write the stuff to the file.
		########################

		echo "#${DATASET}" > ${py_file}
		echo "import FWCore.ParameterSet.Config as cms" >> ${py_file}
		echo "readFiles = cms.untracked.vstring()" >> ${py_file}
		echo "secFiles = cms.untracked.vstring()" >> ${py_file}
		echo 'source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)' >> ${py_file}
		echo "readFiles.extend( [" >> ${py_file}
		
		set files = `ls -1 ${ABSDIR}/*.root`
		set i = 1
		foreach file ( $files )
			if ( $i == 1) then
				echo "'${OUTDIR}/${file}'" >> ${py_file}
			else
				echo ",'${OUTDIR}/${file}'" >> ${py_file}
			endif
			@ i+=1
		end


		echo "] )" >> ${py_file}

		echo "${RED_F}<<<<<<<<< END Report for ${BLUE_F}${RESDIR} ${NORM}"
end


