#! /bin/csh
#to check condor job output
set curDir = `pwd`
setenv CONDOR_JOB_DIR /share/store/users/samantha/CMSSW_5_2_5/src/UserCode/RA2QCDvetoAna/flatrun/condor_submit
echo $CONDOR_JOB_DIR

if ! $?CONDOR_JOB_DIR then
	echo "Environment variable CONDOR_JOB_DIR is not set!"
	exit 1
endif

if ( ! -d $CONDOR_JOB_DIR ) then
	echo "Environment variable CONDOR_JOB_DIR does not point to a valid  directory!"
	exit 1
endif

#check if scripts exists
if  ( ! -e $CONDOR_JOB_DIR/renameQCDfiles.csh ) then
	echo "Script named renameQCDfiles.csh not found in $CONDOR_JOB_DIR!"
	exit 1
endif


foreach dir ( $argv )

	if ( -d $dir ) then
		set curdir = `pwd`
		cd -v $dir
		#go one level up before renaming
		cd ..
		$CONDOR_JOB_DIR/condorJobTest.csh $dir
		$CONDOR_JOB_DIR/renameQCDfiles.csh
		root -b -q $CONDOR_JOB_DIR/haddws.C+
		root -b -q  $CONDOR_JOB_DIR/mergeSearchBins.C+
		cd curdir
	else 
		echo "$dir is not a directory!. skipping.."
	endif
end
