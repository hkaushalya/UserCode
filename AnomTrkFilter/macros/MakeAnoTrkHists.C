{
	TFile *f = new TFile("METFilter_METPD.root");
	const float rebin = 4;

	/*
	TH1* hist4_AtanVsMet = dynamic_cast<TH1*> (f->FindObjectAny("Atan2VsMet_iter5"));
	assert (hist4_AtanVsMet != NULL && "hist4_AtanVsMet not found!");
	new TCanvas();
	hist4_AtanVsMet->Draw();
	
	TH1* hist4_RatioIt4toIt01 = dynamic_cast<TH1*> (f->FindObjectAny("hRatioIt4toIt01"));
	assert(hist4_RatioIt4toIt01 != NULL && "hist4_RatioIt4toIt01");
	new TCanvas();
	gPad->SetLogy();
	hist4_RatioIt4toIt01->Draw();

	TH1* hist5_RatioIt5toIt01 = dynamic_cast<TH1*> (f->FindObjectAny("hRatioIt5toIt01"));
	assert(hist5_RatioIt5toIt01 != NULL && "hist5_RatioIt5toIt01");
	new TCanvas();
	gPad->SetLogy();
	hist5_RatioIt5toIt01->Draw();
	*/


	TH1* histPFmet = dynamic_cast<TH1*> (f->FindObjectAny("PfMet_iter4"));
	assert (histPFmet != NULL && "histPFmet not found!");
	TH1* histPFht = dynamic_cast<TH1*> (f->FindObjectAny("PfHt_iter4"));
	assert (histPFht != NULL && "histPFht not found!");


	TH1* histPFmet_it4lt7_pure = dynamic_cast<TH1*> (f->FindObjectAny("PfMet_highPureTrks_atanLt7_iter4"));
	assert (histPFmet_it4lt7_pure != NULL && "histPFmet_it4lt7_pure not found!");
	TH1* histPFmet_it4gt7_pure = dynamic_cast<TH1*> (f->FindObjectAny("PfMet_highPureTrks_atanGt7_iter4"));
	assert (histPFmet_it4gt7_pure != NULL && "histPFmet_it4gt7_pure not found!");


	TH1* histPFht_it4lt7_pure = dynamic_cast<TH1*> (f->FindObjectAny("PfHt_highPureTrks_atanLt7_iter4"));
	assert (histPFht_it4lt7_pure != NULL && "histPFht_it4lt7_pure not found!");
	TH1* histPFmet_it5lt7_pure = dynamic_cast<TH1*> (f->FindObjectAny("PfMet_highPureTrks_atanLt7_iter5"));
	assert (histPFmet_it5lt7_pure != NULL && "histPFmet_it5lt7_pure not found!");
	TH1* histPFht_it5lt7_pure = dynamic_cast<TH1*> (f->FindObjectAny("PfHt_highPureTrks_atanLt7_iter5"));
	assert (histPFht_it5lt7_pure != NULL && "histPFht_it5lt7_pure not found!");

	histPFmet->Rebin(rebin);
	histPFmet_it4lt7_pure->Rebin(rebin);
	histPFmet_it5lt7_pure->Rebin(rebin);
	histPFht->Rebin(rebin);
	histPFht_it4lt7_pure->Rebin(rebin);
	histPFht_it5lt7_pure->Rebin(rebin);

	histPFmet_it4lt7_pure->SetLineColor(kRed);
	histPFmet_it5lt7_pure->SetLineColor(kBlue);
	histPFht_it4lt7_pure->SetLineColor(kRed);
	histPFht_it5lt7_pure->SetLineColor(kBlue);


	histPFmet->SetLineWidth(2);
	histPFmet_it4lt7_pure->SetLineWidth(2);
	histPFmet_it5lt7_pure->SetLineWidth(2);
	histPFht->SetLineWidth(2);
	histPFht_it4lt7_pure->SetLineWidth(2);
	histPFht_it5lt7_pure->SetLineWidth(2);

	TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
//	leg->AddEntry(histPFmet,"PFMET");
//	leg->AddEntry(histPFmet_it4lt7_pure,"PFMET (atan2(it 4/it0+1) < 0.7)");
//	leg->AddEntry(histPFmet_it5lt7_pure,"PFMET (atan2(it 5/it0+1) < 0.7)");


	histPFmet->SetTitle("Subtracted PFMET (Rejected #slash{E}_{T} Shape); PFMET ;Events");

	new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	//histPFmet->Add(histPFmet_it4lt7_pure,-1.0);
	TH1* histPFmet_copy1 = dynamic_cast<TH1*> (histPFmet->Clone("histPDFmet_copy1"));
	TH1* histPFmet_copy2 = dynamic_cast<TH1*> (histPFmet->Clone("histPDFmet_copy2"));
	histPFmet_copy1->SetLineColor(kBlue);
	histPFmet_copy2->SetLineColor(kRed);
	histPFmet_copy1->Add(histPFmet_it4lt7_pure, -1);
	histPFmet_copy2->Add(histPFmet_it5lt7_pure, -1);
	leg->AddEntry(histPFmet,"PFMET");
	leg->AddEntry(histPFmet_copy1,"PFMET (atan2(it 4/it0+1) > 0.7)");
	leg->AddEntry(histPFmet_copy2,"PFMET (atan2(it 5/it0+1) > 0.7)");
	histPFmet->DrawClone();
	histPFmet_copy1->DrawClone("same");
	histPFmet_copy2->DrawClone("same");
	leg->Draw();

	std::cout << histPFmet->GetName() << ":"<< histPFmet->GetEntries() << ", " << histPFmet->Integral() << std::endl;
	std::cout << histPFmet_it4lt7_pure->GetName() << ":"<< histPFmet_it4lt7_pure->GetEntries() << ", " << histPFmet_it4lt7_pure->Integral() << std::endl;
	std::cout << histPFmet_it4gt7_pure->GetName() << ":"<< histPFmet_it4gt7_pure->GetEntries() << ", " << histPFmet_it4gt7_pure->Integral() << std::endl;

/*return;
	//check the diff in what is removed
	histPFmet_it4lt7_pure->Divide(histPFmet);
	histPFmet_it5lt7_pure->Divide(histPFmet);
	new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	histPFmet_it4lt7_pure->SetMarkerStyle(3);
	histPFmet_it4lt7_pure->SetMarkerStyle(2);
	histPFmet_it4lt7_pure->DrawClone("P");
	histPFmet_it5lt7_pure->DrawClone("Psame");
	leg->Draw();
	



	TLegend *leg1 = new TLegend(0.5,0.7,0.9,0.9);
	leg1->AddEntry(histPFht,"PFHT");
	leg1->AddEntry(histPFht_it4lt7_pure,"PFHT (atan2(it 4/it0+1) < 0.7)");
	leg1->AddEntry(histPFht_it5lt7_pure,"PFHT (atan2(it 5/it0+1) < 0.7)");


	histPFht->SetTitle(";PFHT;Events");
*/
/*	new TCanvas();
	gPad->SetLogy();
	gStyle->SetOptStat(0);
	histPFht->Draw();
	histPFht_it4lt7_pure->Draw("same");
	histPFht_it5lt7_pure->Draw("same");
	leg1->Draw();
*/


/*
	new TCanvas();
	gPad->SetLogy();
	TH1* histAtanRatioIt4to01_pure = dynamic_cast<TH1*> (f->FindObjectAny("hRatioIt4toIt01pure"));
	assert (histAtanRatioIt4to01_pure != NULL && "hRatioIt4toIt01 not found!");
	TH1* histAtanRatioIt5to01_pure = dynamic_cast<TH1*> (f->FindObjectAny("hRatioIt5toIt01pure"));
	assert (histAtanRatioIt5to01_pure != NULL && "hRatioIt5toIt01 not found!");

	histAtanRatioIt4to01_pure->SetTitle("atan2 ratio from high purity tracks;atan ratio;Events");
	histAtanRatioIt4to01_pure->SetLineWidth(4);
	histAtanRatioIt5to01_pure->SetTitle("atan2 ratio from high purity tracks;atan ratio;Events");
	histAtanRatioIt5to01_pure->SetLineWidth(2);
	histAtanRatioIt5to01_pure->SetLineColor(4);
	histAtanRatioIt4to01_pure->Draw();
	histAtanRatioIt5to01_pure->Draw();


	TLegend *leg2 = new TLegend(0.5,0.7,0.9,0.9);
	leg2->AddEntry(histAtanRatioIt4to01_pure,"atan2(Iter4 / Iter0+1)");
	leg2->AddEntry(histAtanRatioIt5to01_pure,"atan2(Iter5 / Iter0+1)");
	leg2->Draw();
	*/
}
