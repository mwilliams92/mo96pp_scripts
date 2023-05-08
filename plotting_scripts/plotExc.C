{

	TChain* Chain_Coinc = new TChain ("t1");
	TChain* Chain_Singles = new TChain ("t1");

	int firstRun = 60;
	int lastRun = 78;
	char fname[64];

	for (int i=firstRun; i<=lastRun; i++) {

		if (i==62) continue;

		sprintf(fname,"../root_outputs/Run00%d_Coinc.root", i);
		Chain_Coinc->Add(fname);
	
		sprintf(fname,"../root_outputs/Run00%d_Singles.root", i);
		Chain_Singles->Add(fname);
	}

	TCut Proton_Cut = "(Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5) || (Range_BB10>2300 && Range_BB10<2900 && Energy>Range_BB10)";
	TCut Proton_Cut_Barrel = "Range_BB10>2300 && Range_BB10<2900";
	TCut Proton_Cut_QQQ5 = "Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5 && Energy<12000";
	TCut Time_Cut = "GretinaTDC>1500 && GretinaTDC<3000";
	TCut HighAngle = "Angle>65";
	TCut LowAngle = "Angle<65";
	TCut Weird_Events = "Excitation>-0.2 && Excitation<0.1 && Egamma>500";
	TCut Oxygen_Gamma = "Egamma>6120 && Egamma<6140";
	TCut Oxygen_Gamma2 = "Egamma>6800 && Egamma<7200";
	TCut Oxygen_Gamma3 = "Egamma>6800 && Egamma<7000";
	TCut Egam_778 = "Egamma>774 && Egamma<782";
	TCut Egam_204 = "Egamma>200 && Egamma<207";

	TFile *OutFile = new TFile("Excitation.root", "recreate");

	TCanvas *c1 = new TCanvas("c1");
	c1->SetLogy();
/*
	TH1D *hExc = new TH1D("hExc","Excitation Energy Barrel",2000,-5,15);
	Chain_Singles->Draw("Excitation>>hExc",Proton_Cut_Barrel);
	hExc->SetLineColor(1);
*/
	TH1D *hExc_Gated = new TH1D("hExc_Gated","Excitation Energy Barrel (Gated)",2000,-5,15);
	Chain_Coinc->Draw("Excitation>>hExc_Gated",Proton_Cut_Barrel + Egam_778);
	hExc_Gated->SetLineColor(2);

	TH1D *hExc_Gated2 = new TH1D("hExc_Gated2","Excitation Energy Barrel (Gated)",2000,-5,15);
	Chain_Coinc->Draw("Excitation>>hExc_Gated2",Proton_Cut_Barrel + Egam_204,"same");
	hExc_Gated2->SetLineColor(4);

	TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   	//legend->AddEntry(hExc,"Singles","l");
	legend->AddEntry(hExc_Gated,"^{96}Mo 2_{1}^{+} -> 0_{g.s.}^{+}","l");
	legend->AddEntry(hExc_Gated2,"^{95}Mo 3/2_{1}^{+} -> 5/2_{g.s.}^{+}","l");
   	legend->Draw();

	c1->Write();
	//hExc->Write();
	hExc_Gated->Write();
	hExc_Gated2->Write();

}
