#!/bin/csh

setenv PYTHONHOME `python-config --prefix` 
echo "PYTHONHOME : $PYTHONHOME\n"

set maxNum=$1
set perNum=$2

perl fnalMakeSplitInputs.pl ${perNum}

set i=0

while ( ${i} < ${maxNum} )
bash fnalCondorPrep.sh ${i}
@ i = ${i} + "1"
end

sleep 1s

set i=0

while ( ${i} < ${maxNum} )
condor_submit run_${i}.sh   
@ i = ${i} + "1"
end
