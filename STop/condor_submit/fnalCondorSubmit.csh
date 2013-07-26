#!/bin/csh

setenv PYTHONHOME `python-config --prefix` 
#echo "PYTHONHOME : $PYTHONHOME\n"

echo "fnalCondorSubmit: nARGS = $#argv"
set nargs = $#argv
set maxNum=$1
set perNum=$2
set fileList=$3
set systVariation=$4
#set cutmask=( $argv[4-${nargs}] )
set cutmask=$5

rm -vrf inputList* run_*
echo "fnalCondorSubmit::  Inputs ========" 
echo " maxNum  = ${maxNum}"
echo " perNum  = ${perNum}"
echo " systVariation = ${systVariation}"
echo " filelist= ${fileList}"
echo " cutmask = $cutmask"

#pwd
#ls -l
if (-e ${PWD}/$fileList ) then
	cat $fileList
else	
	echo "fnalCondorSubmit::  filelist= ${fileList} not found!"
	exit 0;
endif

perl fnalMakeSplitInputs.pl ${perNum} ${fileList}

set i=0

while ( ${i} < ${maxNum} )
bash fnalCondorPrep.sh ${i} ${systVariation} ${cutmask}
@ i = ${i} + "1"
end

sleep 1s

set i=0

while ( ${i} < ${maxNum} )
	condor_submit run_${i}.sh   
	@ i = ${i} + "1"
end
