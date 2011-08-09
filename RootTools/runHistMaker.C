{
	gSystem->Load("HistMaker_cc.so");
	
	HistMaker histMaker;

/*	histMaker.AddFile("METFilter_METPD.root","METPD-Tagged",0,2,1);
	histMaker.AddFile("METFilter_JETPD_Bad.root","JETPD-Tagged",0,2,8);
	histMaker.AddFile("METFilter_TTbarMC.root","t#bar{t} MC-Tagged",0,2,9);
	histMaker.AddFile("METFilter_Tais_Bad.root","MET Scanner Found",0,2,2);
*/

	histMaker.AddFile("METFilter_METPD.root","METPD",0,2,1);
	histMaker.AddFile("METFilter_JETPD.root","JETPD",0,2,28);
	histMaker.AddFile("METFilter_TTbarMC.root","t#bar{t} MC",0,2,9);
	histMaker.AddFile("METFilter_Tais_Bad.root","MET Scanner Found",0,2,2);
	histMaker.AddFile("METFilter_JETPD_Bad.root","JETPD-Tagged",0,2,8);

/*	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_algo"," Tracks Associate with Vertices;Algo;;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_pixHits"," Tracks Associate with Vertices;Valid Pixel Hits;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_purity"," Tracks Associate with Vertices;purity;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_pt"," Tracks Associate with Vertices;pt;;",1,1,100);
	
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_algo"," Tracks NOT Associate with Vertices;Algo;;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_pixHits"," Tracks NOT Associate with Vertices;Valid Pixel Hits;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_purity"," Tracks NOT Associate with Vertices;purity;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_pt"," Tracks NOT Associate with Vertices;pt;;",1,1,100);
	

	histMaker.AddHist("anomtrk/nTrksAssoWithVtx2ntrks"," Ratio = N tracks associate with Vertices/N tracks;Ratio;Events;",1,1,100);
	histMaker.AddHist("anomtrk/nopixhits2ntrks"," Ratio = N tracks with no pixel hits/N tracks;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx2ntrks"," Ratio = N tracks NOT associate with Vertices/N tracks;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/AfterNtrkCut_ratio","ratio after N track cut;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/AfterNtrkCut_ratio2"," ratio2 after N track cut;Ratio;Events;",1,1,100);
	histMaker.AddHist("anomtrk/AfterNtrkCut_ratio3"," ratio3 after N track cut;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/AfterNtrkRatioCut_ratio2"," ratio2 after N track+ratio cut;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/AfterNtrkRatioCut_ratio3"," ratio3 after N track+ratio cut;Ratio;Events;",1,1,1000);
	histMaker.AddHist("anomtrk/AfterNtrkRatioRatio2Cut_ratio3"," ratio3 after N track+ratio+ratio2 cut;Ratio;Events;",1,1,1000);


	histMaker.AddHist("anomtrk/trkPurityFraction","Purity Fraction in Tracks Associate with Vertices;purity fraction;;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_wgt","Weights of the Tracks Associate with Vertices;weights;;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksWgt_lt1","# of tracks with weights <0.1 more than 10% of Tracks Associate with Vertices;# of tracks;;",4,1,100);
	histMaker.AddHist("anomtrk/nTrksWgt_lt2","# of tracks with weights <0.2 more than 20% of Tracks Associate with Vertices;# of tracks;;",4,1,100);

	histMaker.AddHist("anomtrk/nTrksWgt_lt1_pureTrkFrac","pure track fraction from tracks with weights <0.1 more than 10% of Tracks Associate with Vertices; pure track fraction;;",2,1,100);
	histMaker.AddHist("anomtrk/nTrksWgt_lt2_pureTrkFrac","pure track fraction from tracks with weights <0.2 more than 20% of Tracks Associate with Vertices;pure track fraction;;",2,1,100);

	histMaker.AddHist("anomtrk/nTrksWgt_lt1_algo","algo of tracks from # of tracks with weights <0.1 more than 10% of Tracks Associate with Vertices;algo;;",1,1,1000);
	histMaker.AddHist("anomtrk/nTrksWgt_lt2_algo","algo of tracks from # of tracks with weights <0.2 more than 20% of Tracks Associate with Vertices;algo;;",1,1,1000);

	histMaker.AddHist("anomtrk/nTrksWgt_lt1_validPixHits","Valid Pix Hits in tracks from # of tracks with weights <0.1 more than 10% of Tracks Associate with Vertices;valid pix hits;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksWgt_lt2_validPixHits","Valid Pix Hits in tracks from # of tracks with weights <0.2 more than 20% of Tracks Associate with Vertices;valid pix hits;;",1,1,100);

	histMaker.AddHist("anomtrk/nTrksWgt_lt1_valid","Valid Pix Hits in tracks from # of tracks with weights <0.1 more than 10% of Tracks Associate with Vertices;valid pix hits;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksWgt_lt2_validPixHits","Valid Pix Hits in tracks from # of tracks with weights <0.2 more than 20% of Tracks Associate with Vertices;valid pix hits;;",1,1,100);


	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_ptErr","#Delta p_{T}/p_{T} of the Tracks Associate with Vertices;;;",4,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_dxyErr","dxyErr of the Tracks Associate with Vertices;;;",2,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_d0Err"," d0Err of the Tracks Associate with Vertices;;;",2,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_exptInnerHits","Expected Inner Hits of the Tracks Associate with Vertices;;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_exptOuterHits","Expected Outer Hits of the Tracks Associate with Vertices;;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_ptVsptErr","#Delta p_{T}/p_{T} vs p_{T} of the Tracks Associate with Vertices;;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_validHitFraction","Valid Hit Fraction of the Tracks Associate with Vertices;;;",2,1,1000);

	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_ptErr","#Delta p_{T}/p_{T} of the Tracks Not Associate with Vertices;;;",4,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_dxyErr","dxyErr of the Tracks Not Associate with Vertices;;;",2,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_d0Err"," d0Err of the Tracks Not Associate with Vertices;;;",2,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_exptInnerHits","Expected Inner Hits of the Tracks Not Associate with Vertices;;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_exptOuterHits","Expected Outer Hits of the Tracks Not Associate with Vertices;;;",1,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_ptVsptErr","#Delta p_{T}/p_{T} vs p_{T} of the Tracks Not Associate with Vertices;;;",2,0,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_validHitFraction","Valid Hit Fraction of the Tracks Not Associate with Vertices;;;",2,1,1000);

	histMaker.AddHist("anomtrk/nTrksAssoWithVtx_ntrks","Number of Tracks Associate with Vertices;N-tracks;;",8,1,100);
	histMaker.AddHist("anomtrk/nTrksNotAssoWithVtx_ntrks","Number of Tracks Not Associate with Vertices;N-tracks;;",8,1,100);
	*/
/*
	histMaker.AddHist("anomtrk/my_PFMet","Tagged Events;PFMET [GeV];Events;",4,1,100);
	histMaker.AddHist("anomtrk/my_CaloMet","Tagged Events;CaloMET [GeV];Events;",4,1,100);
	histMaker.AddHist("anomtrk/my_CaloPFMet","Tagged Events;(CaloMET-PFMET)/(CaloMET+PFMET);Events;",4,1,1000);
	histMaker.AddHist("anomtrk/my_Calo2PFMet","Tagged Events;(CaloMET-PFMET)/PFMET;Events;",8,1,100);;
*/
//	histMaker.AddHist("anomtrk/all_PFMet","PFMET All Events;;;",4,1,100);
	//histMaker.AddHist("anomtrk/andrea_PFMet","PFMET from Andrea event selection;;;",4,1,100);

	histMaker.AddHist("anomtrk/all_PFMet","All Events;PFMET [GeV];Events;",4,1,100);
	histMaker.AddHist("anomtrk/all_CaloMet","All Events;CaloMET [GeV];Events;",4,1,100);
	histMaker.AddHist("anomtrk/all_CaloPFMet","All Events;(CaloMET-PFMET)/(CaloMET+PFMET);Events;",4,1,10000);
	histMaker.AddHist("anomtrk/all_Calo2PFMet","All Events;(CaloMET-PFMET)/PFMET;Events;",8,1,1000);;
	histMaker.PrintFormat("eps");
	
	histMaker.DrawAll();
/*	histMaker.Draw(0);
	histMaker.Draw(4);
	histMaker.Draw(1);
	histMaker.Draw(5);
	histMaker.Draw(2);
	histMaker.Draw(6);
	histMaker.Draw(3);
	histMaker.Draw(7);


	histMaker.Draw(8);
	histMaker.Draw(9);
	histMaker.Draw(10);
	histMaker.Draw(11);
	histMaker.Draw(12);
	histMaker.Draw(13);
	histMaker.Draw(14);
	histMaker.Draw(15);
	histMaker.Draw(16);
*/

}
