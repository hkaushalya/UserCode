#!/bin/bash

echo "Executing fnalCondorPrep.sh"

export configSpecStr=$1
export dataset=$2
export outputrootfilepath=$3



export filetail="_"$configSpecStr
export currDir=`pwd`

#export outputDir=${currDir}'/'
export outputDir= ${_CONDOR_SCRATCH_DIR}
export logDir=$outputDir'log'

echo $configSpecStr
echo $filetail
echo $outputDir
echo $currDir

export inputFileName='inputList_'$configSpecStr'.rtlst'
export inputlist=$inputFileName

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
+LENGTH="SHORT"
getenv = true
initialdir = ${currDir}
Output = ${configSpecStr}_\$(Cluster)_\$(Process).stdout
Error = ${configSpecStr}_\$(Cluster)_\$(Process).stderr
Log = ${configSpecStr}_\$(Cluster)_\$(Process).log
transfer_input_files = $currDir/runLostLeptonTreeMaker_condor.py, $currDir/$inputlist, $currDir/DataPileupHistogram_RA2Summer12_190456-208686_ABCD.root 
Arguments = runLostLeptonTreeMaker_condor.py $filetail $inputlist \$(Cluster) ${configSpecStr} ${dataset} ${_CONDOR_SCRATCH_DIR}


Queue 

#exit
EOF
#----------

chmod 755 $filename

#condor_submit xxx.sh
#condor_q -submitter lhx
#condor_rm -name cmslpc15 jobID

exit

