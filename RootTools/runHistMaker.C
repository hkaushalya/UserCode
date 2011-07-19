{
	gSystem->Load("HistMaker_cc.so");
	
	HistMaker histMaker;

	histMaker.AddFile("METFilter_METPDSkim.root","METPD",0,2,1);
	histMaker.AddFile("METFilter_METPDSkim_Bad.root","METPD Bad",0,2,6);
	histMaker.AddFile("METFilter_TTbarMC.root","t#bar{t} MC",0,2,9);

	histMaker.AddHist("anomtrk2/nTrksAssoWithVtx_algo"," Tracks Associate with Vertices;Algo;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksAssoWithVtx_pixHits"," Tracks Associate with Vertices;Valid Pixel Hits;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksAssoWithVtx_purity"," Tracks Associate with Vertices;purity;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksAssoWithVtx_pt"," Tracks Associate with Vertices;pt;;",1,100);
	
	histMaker.AddHist("anomtrk2/nTrksNotAssoWithVtx_algo"," Tracks NOT Associate with Vertices;Algo;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksNotAssoWithVtx_pixHits"," Tracks NOT Associate with Vertices;Valid Pixel Hits;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksNotAssoWithVtx_purity"," Tracks NOT Associate with Vertices;purity;;",1,100);
	histMaker.AddHist("anomtrk2/nTrksNotAssoWithVtx_pt"," Tracks NOT Associate with Vertices;pt;;",1,100);
	

	histMaker.AddHist("anomtrk2/nTrksAssoWithVtx2ntrks"," Ratio = N tracks associate with Vertices/N tracks;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/nopixhits2ntrks"," Ratio = N tracks with no pixel hits/N tracks;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/nTrksNotAssoWithVtx2ntrks"," Ratio = N tracks NOT associate with Vertices/N tracks;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkCut_ratio","ratio after N track cut;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkCut_ratio2"," ratio2 after N track cut;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkCut_ratio3"," ratio3 after N track cut;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkRatioCut_ratio2"," ratio2 after N track+ratio cut;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkRatioCut_ratio3"," ratio3 after N track+ratio cut;Ratio;Events;",1,100);
	histMaker.AddHist("anomtrk2/AfterNtrkRatioRatio2Cut_ratio3"," ratio3 after N track+ratio+ratio2 cut;Ratio;Events;",1,100);



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
