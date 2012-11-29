#!/bin/csh

setenv PYTHONHOME `python-config --prefix` 
echo "PYTHONHOME : $PYTHONHOME\n"

set maxNum=$1
set perNum=$2
set fileList=$3
set dataset=$4

rm -vrf inputList* run_*
echo "fnalCondorSubmit::  Inputs ========" 
echo " maxNum  = ${maxNum}"
echo " perNum  = ${perNum}"
echo " dataset = ${dataset}"
echo " filelist= ${fileList}"

pwd
ls
if (-e $fileList ) then
	cat $fileList
else	
	echo "fnalCondorSubmit::  filelist= ${fileList} not found!"
	exit 0;
endif

perl fnalMakeSplitInputs.pl ${perNum} ${fileList}

set i=0

while ( ${i} < ${maxNum} )
bash fnalCondorPrep.sh ${i} ${dataset}
@ i = ${i} + "1"
end

sleep 1s

set i=0

while ( ${i} < ${maxNum} )
condor_submit run_${i}.sh   
@ i = ${i} + "1"
end
