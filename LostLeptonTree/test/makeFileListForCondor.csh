#! /bin/tcsh

#this is the create a complete list of files 
#in EOS OR PNFS directory to be run over by a job
#in CONDOR. This will create perl file which will
#the input to fnalMakeSplitInputs.pl 

#LOCATION 1=EOS else=PNFS
set LOCATION=2  #default set to pnfs

#set COMMONDIR="2011RA2/QCD_TuneZ2_7TeV_pythia6/QCD_Pt_1000to1400"

#set OUTFILE="file.list"
#set OUTFILE="HT_Run2012A_13Jul2012_v1_lpc1.list"
#set OUTFILE="HTMHT_Run2012B_13Jul2012_v1_lpc1.list"
#set OUTFILE="HT_Run2012Arecover_06Aug2012_v1_lpc1.list"
#set OUTFILE="HTMHT_Run2012C_PromptReco_v2_lpc1.list"
#set OUTFILE="HTMHT_Run2012C_24Aug2012_v1_lpc1.list"
#set OUTFILE="ZJetsToNuNu_400_HT_inf_8TeV_madgraph_Summer12_v1"
#set OUTFILE="QCD_HT_250To500_MGPythia_v1_lpc1.files"
#set OUTFILE="QCD_HT_500To1000_MGPythia_v1_lpc1.files"
#set OUTFILE="QCD_HT_1000ToInf_MGPythia_v1_lpc1.files"
#set OUTFILE="TT_CT10_TuneZ2star_8TeV-powheg-tauola.list"

#set OUTFILE="QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_1000to1400_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_1400to1800_TuneZ2star_8TeV_pythia6_Summer12_53X.list"
#set OUTFILE="QCD_Pt_1800_TuneZ2star_8TeV_pythia6_Summer12_53X.list"

#jet PD dataset. because some of the HT triggers were moved to JetPD
#set OUTFILE = "JetHT_Run2012B_13Jul2012_v1_lpc1.list"
#set OUTFILE = "JetHT_Run2012C_24Aug2012_v2_lpc1.list"
#set OUTFILE = "JetHT_Run2012C_PromptReco_v2_lpc1.list"

#set OUTFILE = "HTMHT_Run2012D_PromptReco_v1_lpc1.list"
#set OUTFILE = "JetHT_Run2012D_PromptReco_v1_lpc1.list"
set OUTFILE = "WJetsToLNu_HT-300To400_8TeV-madgraph_kasmi.list"

rm -rf ${OUTFILE}

if ( ${LOCATION} == 1 ) then
#	set SRCDIR="/eos/uscms/store/user/samantha/${COMMONDIR}"
else 
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HT_Run2012A_13Jul2012_v1_lpc1/seema/HT/HT_Run2012A-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012B_13Jul2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012B-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HT_Run2012Arecover_06Aug2012_v1_lpc1/seema/HT/HT_Run2012A-recover-06Aug2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/8e90e164b25bf8451f454d8c6a3799cd"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_PromptReco_v2_lpc1/seema/HTMHT/HTMHT_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/993e6ea01263feda15a347a7daafedb6"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012C_24Aug2012_v1_lpc1/seema/HTMHT/HTMHT_Run2012C-24Aug2012-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/5752ae8db176959272f91d699d678612"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_8TeV_madgraph_Summer12_v1/seema/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/ZJetsToNuNu_400_HT_inf_8TeV_madgraph_START53_V7A_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/QCD_HT_250To500_MGPythia_v1_lpc1/seema/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/QCD_HT_500To1000_MGPythia_v1_lpc1/seema/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd"
	#set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/QCD_HT_1000ToInf_MGPythia_v1_lpc1/seema/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd" 
#	set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/kasmi/TT_CT10_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v2_NOCUTS_12Oct2012V3/0c2ca37b7b411866265a5953adbfd963"

	set SRCDIR="/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/kasmi/WJetsToLNu_HT-300To400_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v2_NOCUTS_12Oct2012V3/0c2ca37b7b411866265a5953adbfd963"

	#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_1000to1400_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_1400to1800_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"
#	set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/samantha/QCD_Pt_1800_TuneZ2star_8TeV_pythia6_Summer12/samantha/QCD_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X_PU_S10_START53_V7A_53X_HEPtagger_OCT262012/0c2ca37b7b411866265a5953adbfd963"

#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012B_13Jul2012_v1_lpc1/seema/JetHT/JetHT_Run2012B-13Jul2012-v1_NOCUTS_HLTPFHTInc_12Oct2012V3/32ed558ccb6d6c1922414263922f006c"
#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012C_24Aug2012_v2_lpc1/seema/JetHT/JetHT_Run2012C-24Aug2012-v2_NOCUTS_HLTPFHTInc_12Oct2012V3/25f502adba324b5e7eb053c17d4f5df5"
#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012C_PromptReco_v2_lpc1/seema/JetHT/JetHT_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/993e6ea01263feda15a347a7daafedb6"

#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/HTMHT_Run2012D_PromptReco_v1_lpc1/seema/HTMHT/HTMHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9"
#set SRCDIR = "/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9"

endif


if ( ! -e $SRCDIR ) then 
	echo "SRCDIR not found ! $SRCDIR"
	exit
endif

set FILES=`dcls -1 ${SRCDIR}`


foreach file ( ${FILES} )
#	if ( ${LOCATION} == 1 ) then
#		echo "root://cmseos.fnal.gov/${SRCDIR}/${file}" >> ${OUTFILE}
#	else 
		#echo "dcap:/${SRCDIR}/${file}"
		echo "dcap:/${SRCDIR}/${file}" >> ${OUTFILE}
#	endif
end

cat ${OUTFILE}
echo ""
echo -n "FILES# "
cat ${OUTFILE} | wc -l 
