#! /bin/csh

set systVariation = 0

set curDir = `pwd`
#set mainDir = "$curDir/06182013_QCDMG_NotUseGenMetForRecoUnclMet_1bs_IntMediatePointsStudy_Syst${systVariation}/Met50_MG1"
set mainDir = "$curDir/06292013_QCDMG_newHistSchema_IntMediatePointsStudy_n200_MetAbove200_Syst${systVariation}"
#set mainDir = "$curDir/06182013_QCDMG_NotUseGenMetForRecoUnclMet_1bs_IntMediatePointsStudy_Syst${systVariation}/"
#set mainDir = "$curDir/06182013_QCDMG_NotUseGenMetForRecoUnclMet_2bs_4k_Syst${systVariation}"
#set mainDir = "$curDir/06092013_QCDMG_NoCut_SmallFileTest2_Syst${systVariation}"
#set mainDir = "$curDir/06052013_QCDMG_Cut255_test${systVariation}"
#set mainDir = "$curDir/05302013_TTbar_CutTests_withAllVarsLikeUnclMET"

if (! -d $mainDir ) then
	mkdir -p  $mainDir
else
	echo "$mainDir  already exists!"
#	exit
endif

set MG1 = "../input_files/MG_QCD_HT_250To500_wFilters_farm_v03012013.files"
#set MG2 = "../input_files/MG_QCD_HT_500To1000_wFilters_farm_v03012013.files"
#set MG3 = "../input_files/MG_QCD_HT_1000ToInf_wFilters_farm_v03012013.files"
#set TT  = "../input_files/TT_CT10_wFilters_farm_v03012013.files"
#set MG1 = "../input_files/MG_QCD_HT_250To500_smallFiles_wFilters_farm_v03012013.files"
set MG2 = "../input_files/MG_QCD_HT_500To1000_smallFiles_wFilters_farm_v03012013.files"
set MG3 = "../input_files/MG_QCD_HT_1000ToInf_smallFiles_wFilters_farm_v03012013.files"

set DATASETS = "$MG3 $MG2 $MG1"
#set DATASETS = "$MG3 $MG2"
#set DATASETS = "$MG1"
#set DATASETS = "$MG2 $MG3 $TT"
#set DATASETS = "$TT"
#set BITMASKS = "0 1 2 4 8 16 32 64 128 256 3 5" #these are the bitmasks enabling each cut cumulatively 
#set BITMASKS = "0 1 2 3 32 " #these are the bitmasks enabling each cut cumulatively 
#set BITMASKS = "255 383" #these are the bitmasks enabling each cut cumulatively 
#set BITMASKS = "207 239 255 383 " #these are the bitmasks enabling each cut cumulatively 
#set BITMASKS = "135 265"
#set BITMASKS = "131 129 139 257 259 267"
set BITMASKS = "1:2:3:128:129:130:131:195:256:257:258:259:323"


#foreach bitmask ( $BITMASKS ) #these are the bitmasks enabling each cut cumulatively
	foreach d ( $DATASETS ) #DATASETS

			set subdir = `echo $d | cut -c 16-34`
			echo "subdir = $subdir"
			#set ttdir = "${mainDir}/$bitmask/$subdir"
			set ttdir = "${mainDir}/$subdir"
			if (! -d $ttdir ) then
				mkdir -vp $ttdir
			else 
				echo "$ttdir already exists! Exiting!"
				exit 1
			endif

			cp -v fnalCondorPrep.sh $ttdir
			cp -v fnalCondorSubmit.csh $ttdir
			cp -v fnalMakeSplitInputs.pl $ttdir
			cp -v ../runsmear $ttdir
			cp $d $ttdir
			cd $ttdir
			set fileList =  `ls -1 *.files`
			set nfiles = `cat $fileList | wc -l`
			#echo "$PWD -> $fileList"
			#./fnalCondorSubmit.csh $nfiles 1  $fileList $systVariation $bitmask 
			./fnalCondorSubmit.csh $nfiles 1  $fileList $systVariation $BITMASKS 
			cd $curDir
			sleep 2
	end
#end

cd $curDir
