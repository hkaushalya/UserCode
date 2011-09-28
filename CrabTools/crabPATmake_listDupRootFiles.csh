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

ls -1 ${1} | grep "root" | cut -d "_" -f3 | sort -n | uniq -d
set dupIds = `ls -1 ${1} | grep "root" | cut -d "_" -f3 | sort -n | uniq -d `
set ndups = `ls -1 ${1} | grep "root" | cut -d "_" -f3 | sort -n | uniq -d | wc -l`
echo "Duplicate for $ndups files found!"

foreach id ( $dupIds )
	echo ">>>>>>>>>>>>>>>>>>>>>>>>>> id = ${id}"
	set files = `ls -1 ${1}"/QCD_PtBinned_"${id}_*`
	ls -ltrh ${1}"/QCD_PtBinned_"${id}_*
	foreach file ( $files )
		rm -i  $file
	end
end
