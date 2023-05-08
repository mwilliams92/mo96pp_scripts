{

	TChain* Chain = new TChain ("t1");
  
	Chain->Add("./Run0060_Coinc.root");
	Chain->Add("./Run0061_Coinc.root");
	Chain->Add("./Run0062_Coinc.root");
/*	Chain->Add("./Run0063_Sorted.root");
	Chain->Add("./Run0064_Sorted.root");
	Chain->Add("./Run0065_Sorted.root");
	Chain->Add("./Run0066_Sorted.root");
	Chain->Add("./Run0067_Sorted.root");
	Chain->Add("./Run0068_Sorted.root");
	Chain->Add("./Run0069_Sorted.root");
	Chain->Add("./Run0070_Sorted.root");
	Chain->Add("./Run0071_Sorted.root");
	Chain->Add("./Run0072_Sorted.root");
	Chain->Add("./Run0073_Sorted.root");
	Chain->Add("./Run0074_Sorted.root");
	Chain->Add("./Run0075_Sorted.root");
	Chain->Add("./Run0076_Sorted.root");
	Chain->Add("./Run0077_Sorted.root");
	Chain->Add("./Run0078_Sorted.root");
*/
	TCut Proton_Cut = "(Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5) || (Range_BB10>2300 && Range_BB10<2900 && Energy>Range_BB10)";
	TCut Proton_Cut_Barrel = "Range_BB10>2300 && Range_BB10<2900 && Energy>3000";
	TCut Proton_Cut_QQQ5 = "Range_QQQ5>2700 && Range_QQQ5<4200 && Energy>Range_QQQ5 && Energy<12000";
	TCut Time_Cut = "GretinaTDC>1500 && GretinaTDC<3000";
	TCut HighAngle = "Angle>65";
	TCut LowAngle = "Angle<65";
	TCut Weird_Events = "Excitation>-0.2 && Excitation<0.1 && Egamma>500";

	TCanvas *c1 = new TCanvas();
	TH2D * hExc_Egam;  // Get the Histogram out
  	t1->GetObject("hExc_Egam",hExc_Egam);
  	hExc_Egam->Draw();

/*
	TCanvas *c1 = new TCanvas();
	TH2D *hExc_Egam = new TH2D("hExc_Egam","Excitation vs Gamma Energy",10000,0,10000,2500,-5,20);
	Chain->Draw("Excitation:Egamma>>hExc_Egam","","colz");
	
	TCanvas *c2 = new TCanvas();
	TH1D *hExc = new TH1D("hExc","Excitation Energy Barrel",500,-5,15);
	Chain->Draw("Excitation>>hExc",Proton_Cut_Barrel);
	TH1D *hExc_Gated = new TH1D("hExc_Gated","Excitation Energy Barrel (Gated)",500,-5,15);
	Chain->Draw("Excitation>>hExc_Gated",Proton_Cut_Barrel + Time_Cut,"same");
	hExc_Gated->SetLineColor(2);

	TH1D *hExc_Gated2 = new TH1D("hExc_Gated2","Excitation Energy Barrel (Gated)",500,-5,15);
	Chain->Draw("Excitation>>hExc_Gated2",Proton_Cut_Barrel + HighAngle + Time_Cut,"same");
	hExc_Gated2->SetLineColor(3);

	TH1D *hExc_Gated3 = new TH1D("hExc_Gated3","Excitation Energy Barrel (Gated)",500,-5,15);
	Chain->Draw("Excitation>>hExc_Gated3",Proton_Cut_Barrel + LowAngle + Time_Cut,"same");
	hExc_Gated3->SetLineColor(4);

	TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   	legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   	legend->AddEntry(hExc,"Gated on Protons","l");
   	legend->AddEntry(hExc_Gated,"Gated on Protons and GretinaTDC>0","l");
	legend->AddEntry(hExc_Gated2,"Gated on Protons and Angle>65 deg","l");
	legend->AddEntry(hExc_Gated3,"Gated on Protons and Angle<65 deg","l");
   	legend->Draw();

	TCanvas *c3 = new TCanvas();
	TH2D *hKin = new TH2D("hKin","Energy vs Angle",180,0,180,2000,0,20000);
	Chain->Draw("Energy:Angle>>hKin",Proton_Cut_Barrel,"colz");

	TCanvas *c4 = new TCanvas();
	TH2D *hKin2 = new TH2D("hKin2","Excitation vs Angle",180,0,180,2000,0,20);
	Chain->Draw("Excitation:Angle>>hKin2",Proton_Cut_Barrel,"colz");

	//TCanvas *c5 = new TCanvas();
	//TH2D *hKin3 = new TH2D("hKin3","Excitation vs Angle",180,0,180,2000,0,20);
	//Chain->Draw("Excitation_Cor:Angle>>hKin3",Proton_Cut_Barrel,"colz");

/*
	TCanvas *c3 = new TCanvas();
	TH1D *hExc2 = new TH1D("hExc2","Excitation Energy QQQ5",500,-5,15);
	Chain->Draw("Excitation>>hExc2",Proton_Cut_QQQ5);
	TH1D *hExc_Gated2 = new TH1D("hExc_Gated2","Excitation Energy QQQ5 (Gated)",500,-5,15);
	Chain->Draw("Excitation>>hExc_Gated2",Proton_Cut_QQQ5 + Time_Cut,"same");
	hExc_Gated2->SetLineColor(2);

	TCanvas *c5 = new TCanvas();
	TH2D *hBarrel_PID = new TH2D("hBarrel_PID","Barrel PID",2000,0,20000,2000,0,20000);
	Chain->Draw("Range_BB10:Energy>>hBarrel_PID","Range_BB10>0 && Energy>0","colz");

	TCanvas *c5 = new TCanvas();
	TH2D *hQQQ5_PID = new TH2D("hQQQ5_PID","QQQ5 PID",2000,0,20000,2000,0,20000);
	Chain->Draw("Range_QQQ5:Energy>>hQQQ5_PID","Range_QQQ5>0 && Energy>0 && GretinaTDC>0","colz");
*/
}
