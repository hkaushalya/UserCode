#!/bin/bash

export configSpecStr=$1
export dataset=$2
export outputrootfilepath=$3

export filetail="_"$configSpecStr

export currDir=`pwd`

#export outputDir=${currDir}'/qcd1/'
export outputDir=${currDir}'/'
export logDir=$outputDir'log'

echo $configSpecStr
echo $filetail
echo $outputDir
echo $currDir

export inputFileName='inputList_'$configSpecStr'.rtlst'
export inputlist=$currDir'/'$inputFileName

echo $inputlist

export filename='run_'$configSpecStr'.sh'
if [ -e $filename ]; then
  rm $filename
fi


#----------
cat >> $filename <<EOF
#!/bin/bash


#---

universe = vanilla
Executable = /uscmst1/prod/sw/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_5/bin/slc5_amd64_gcc462/cmsRun 
Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = true
initialdir = ${outputDir}
Output = ${logDir}_${configSpecStr}_\$(Cluster)_\$(Process).stdout
Error = ${logDir}_${configSpecStr}_\$(Cluster)_\$(Process).stderr
Log = ${logDir}_${configSpecStr}_\$(Cluster)_\$(Process).log
Arguments = ${currDir}/runLostLeptonTreeMaker_condor.py $filetail $inputlist \$(Cluster) ${configSpecStr} ${dataset} ${outputrootfilepath}

Queue 

#exit
EOF
#----------

chmod 755 $filename

#condor_submit xxx.sh
#condor_q -submitter lhx
#condor_rm -name cmslpc15 jobID

exit

