{  
  gStyle->SetLineWidth(2);
   TFile *fdata = TFile::Open("../root_outputs/SumRun_59-93.root", "READ");
   //TFile *fdata = TFile::Open("../root_outputs/Run0060_sorted.root", "READ");
   //TFile *fdata = TFile::Open("../root_outputs/SumRun_11-29.root", "READ");
   TTree *tsingles = (TTree*)fdata->Get("t1");
   TTree *tcoinc = (TTree*)fdata->Get("t2");

   tsingles->SetAlias("RangeMeV","Range*0.001");
   tsingles->SetAlias("EtotMeV","Etot*0.001");

   TCut barrel_protons = "Range>2400 && Range<3000 && EndCap==0 && Upstream==0";
   TCut barrel = "EndCap==0 && Upstream==0";

   TCut endcap = "EndCap==1 && Upstream==0";
   TCut endcap_left = "EndCap==1 && Upstream==0 && Left==1";
   TCut endcap_right = "EndCap==1 && Upstream==0 && Left==0";
   TCut endcap_protons = "Range>3000 && Range<4200 && EndCap==1 && Upstream==0";
   TCut endcap_protons_right = "Range>3000 && Range<4200 && EndCap==1 && Upstream==0 && Left==0";
   TCut endcap_protons_left = "Range>3000 && Range<4200 && EndCap==1 && Upstream==0 && Left==1";

   TCut all_protons = "(Range>2400 && Range<3000 && EndCap==0) || (Range>3000 && Range<4200 && EndCap==1)";

   TCut EndCap_Left = "EndCap==1 && Upstream==0 && Left==1";

   TCut angle_cutoff = "ThetaLab>65 && ThetaLab<88";
   TCut angle_cutoff2 = "ThetaLab>32 && ThetaLab<41";
   TCut angle_lol = "ThetaLab>77 && ThetaLab<87";
   TCut avoid_contaminant = "!(ThetaLab>45 && ThetaLab<65)";

   TCut CarbonState = "Egam>4200 && Egam<4600";
   TCut Egam778keV = "Egam>770 && Egam<784";

   TCut DividedTrigger = "tdcSilicon_Div>390 && tdcSilicon_Div<420";

   //TCut DeltaTS_Cut = "deltaTS>20 && deltaTS<100";
   TCut DeltaTS_Cut = "deltaTS>100 && deltaTS<350";

   TCutG *cutg = new TCutG("cutg",9);
   cutg->SetVarX("Etot");
   cutg->SetVarY("Range");
   cutg->SetPoint(0,8.55731,2.28507);
   cutg->SetPoint(1,6.36197,2.23542);
   cutg->SetPoint(2,4.02028,2.29926);
   cutg->SetPoint(3,2.64035,2.40564);
   cutg->SetPoint(4,3.22578,2.80282);
   cutg->SetPoint(5,4.68933,2.92339);
   cutg->SetPoint(6,6.40379,2.92339);
   cutg->SetPoint(7,8.85002,2.89502);
   cutg->SetPoint(8,8.55731,2.28507);

  TCutG *cutg2 = new TCutG("cutg2",10);
  cutg2->SetVarX("Etot");
  cutg2->SetVarY("DeltaE");
  cutg2->SetPoint(0,3051.576,2096.82);
  cutg2->SetPoint(1,4770.774,1531.498);
  cutg2->SetPoint(2,5988.539,1286.719);
  cutg2->SetPoint(3,7313.754,1111.878);
  cutg2->SetPoint(4,9641.834,948.692);
  cutg2->SetPoint(5,8853.868,668.9451);
  cutg2->SetPoint(6,6955.587,768.0221);
  cutg2->SetPoint(7,4161.891,1135.19);
  cutg2->SetPoint(8,2621.776,1752.964);
  cutg2->SetPoint(9,3051.576,2096.82);

  TCut graphical_cut = "cutg";
  TCut graphical_cut2 = "cutg2";

  TCutG *dE_E_EndCap = new TCutG("dE_E_EndCap",13);
  dE_E_EndCap->SetVarX("Etot");
  dE_E_EndCap->SetVarY("DeltaE");
  dE_E_EndCap->SetPoint(0,9479.09,974.729);
  dE_E_EndCap->SetPoint(1,7722.83,1175.63);
  dE_E_EndCap->SetPoint(2,6307.59,1411.17);
  dE_E_EndCap->SetPoint(3,5216.33,1785.27);
  dE_E_EndCap->SetPoint(4,4210.31,2256.36);
  dE_E_EndCap->SetPoint(5,3425.96,2942.2);
  dE_E_EndCap->SetPoint(6,4125.06,3344.01);
  dE_E_EndCap->SetPoint(7,4858.25,2644.31);
  dE_E_EndCap->SetPoint(8,5693.76,2117.8);
  dE_E_EndCap->SetPoint(9,6904.38,1736.78);
  dE_E_EndCap->SetPoint(10,8507.18,1466.6);
  dE_E_EndCap->SetPoint(11,9513.2,1300.33);
  dE_E_EndCap->SetPoint(12,9479.09,974.729);

  TCut Cut_dE_E_EndCap = "dE_E_EndCap";

  TCutG *range_EndCap = new TCutG("range_EndCap",5);
  range_EndCap->SetVarX("Etot");
  range_EndCap->SetVarY("Range");
  range_EndCap->SetPoint(0,10365.8,3369.03);
  range_EndCap->SetPoint(1,10444.3,3817.74);
  range_EndCap->SetPoint(2,4169.21,3855.59);
  range_EndCap->SetPoint(3,4012.34,3401.47);
  range_EndCap->SetPoint(4,10365.8,3369.03);
  TCut Cut_range_EndCap = "range_EndCap";
/*
  /// Particle-ID
  TCanvas *cPID_Barrel = new TCanvas("cPID_Barrel");
  TH2F *hPID_Barrel = new TH2F("hPID_Barrel","Barrel Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("DeltaE:Etot>>hPID_Barrel","EndCap==0 && Upstream==0","colz");
  TH2F *hPID_Barrel_Coinc = new TH2F("hPID_Barrel_Coinc","Barrel Particle ID",2500,0,25000,2000,0,20000);
  //tcoinc->Draw("DeltaE:Etot>>hPID_Barrel_Coinc","EndCap==0 && Upstream==0 && Egam>770 && Egam<784","same");
  //hPID_Barrel_Coinc->SetMarkerStyle(20); hPID_Barrel_Coinc->SetMarkerColor(2); hPID_Barrel_Coinc->SetMarkerSize(0.1);
  hPID_Barrel->GetXaxis()->SetTitle("Total Energy (keV)");
  hPID_Barrel->GetYaxis()->SetTitle("#Delta E Energy (keV)");
  //cutg2->SetLineColor(1); cutg2->SetLineWidth(2);
  //cutg2->Draw("same");

  //gStyle->SetLineWidth(1.5);
  TCanvas *cPID_Barrel2 = new TCanvas("cPID_Barrel2"); cPID_Barrel2->SetLogz(); cPID_Barrel2->SetTicky(); cPID_Barrel2->SetTickx();
  TH2F *hPID_Barrel2 = new TH2F("hPID_Barrel2","Particle ID",2000,0,20,2000,0,20);
  tsingles->Draw("RangeMeV:EtotMeV>>hPID_Barrel2","EndCap==0 && Upstream==0","colz");
  //TH2F *hPID_Barrel2_Coinc = new TH2F("hPID_Barrel2_Coinc","Barrel Particle ID",2500,0,25000,2000,0,20000);
  //tcoinc->Draw("Range:Etot>>hPID_Barrel2_Coinc","EndCap==0 && Upstream==0 && Egam>770 && Egam<784","same");
  //hPID_Barrel2_Coinc->SetMarkerStyle(20); hPID_Barrel2_Coinc->SetMarkerColor(2); hPID_Barrel2_Coinc->SetMarkerSize(0.1);
  hPID_Barrel2->GetXaxis()->SetTitle("Total Energy (MeV)"); hPID_Barrel2->GetXaxis()->SetRangeUser(0,16);
  hPID_Barrel2->GetYaxis()->SetTitle("#DeltaE Range (MeV)"); hPID_Barrel2->GetYaxis()->SetRangeUser(0,5);
  hPID_Barrel2->GetXaxis()->SetTitleSize(0.05); hPID_Barrel2->GetYaxis()->SetTitleSize(0.05);
  hPID_Barrel2->GetXaxis()->CenterTitle(); hPID_Barrel2->GetYaxis()->CenterTitle();
  cutg->SetLineColor(1); cutg->SetLineWidth(2); cutg->SetLineColor(2);
  //cutg->Draw("same");

  TCanvas *cPID_EndCap = new TCanvas("cPID_EndCap");
  TH2F *hPID_EndCap = new TH2F("hPID_EndCap","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("DeltaE:Etot>>hPID_EndCap","EndCap==1 && Upstream==0 && Left==1 && ThetaLab>32 && ThetaLab<41","colz");
  TH2F *hPID_EndCap_Coinc = new TH2F("hPID_EndCap_Coinc","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tcoinc->Draw("DeltaE:Etot>>hPID_EndCap_Coinc","EndCap==1 && Upstream==0 && Left==1 && ThetaLab>32 && ThetaLab<41 && Egam>770 && Egam<784","same");
  hPID_EndCap_Coinc->SetMarkerStyle(20); hPID_EndCap_Coinc->SetMarkerColor(2); hPID_EndCap_Coinc->SetMarkerSize(0.1);
  hPID_EndCap->GetXaxis()->SetTitle("#DeltaE + E1 + E2 Energy (keV)");
  hPID_EndCap->GetYaxis()->SetTitle("#DeltaE (keV)");
  dE_E_EndCap->SetLineColor(1); dE_E_EndCap->SetLineWidth(2);
  dE_E_EndCap->Draw("same");

  TCanvas *cPID_EndCap2 = new TCanvas("cPID_EndCap2");
  TH2F *hPID_EndCap2 = new TH2F("hPID_EndCap2","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("Range:Etot>>hPID_EndCap2","EndCap==1 && Upstream==0 && ThetaLab>32 && ThetaLab<41 && Left==1","colz");
  TH2F *hPID_EndCap2_Coinc = new TH2F("hPID_EndCap2_Coinc","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tcoinc->Draw("Range:Etot>>hPID_EndCap2_Coinc","EndCap==1 && Upstream==0 && Left==1 && ThetaLab>32 && ThetaLab<41 && Egam>770 && Egam<784","same");
  hPID_EndCap2_Coinc->SetMarkerStyle(20); hPID_EndCap2_Coinc->SetMarkerColor(2); hPID_EndCap2_Coinc->SetMarkerSize(0.1);
  hPID_EndCap2->GetXaxis()->SetTitle("#DeltaE + E1 + E2 Energy (keV)");
  hPID_EndCap2->GetYaxis()->SetTitle("Range (keV)");
  range_EndCap->SetLineColor(1); range_EndCap->SetLineWidth(2);
  range_EndCap->Draw("same");

  /// Kinematics
  TCanvas *cKin = new TCanvas("cKin");
  TH2F *hKin = new TH2F("hKin","Kinematics",180,0,180,2500,0,25);
  tsingles->Draw("Etot/1e3:ThetaLab>>hKin","Left==1 || Left==-1","colz");
  hKin->GetXaxis()->SetTitle("Lab Angle (degrees)");
  hKin->GetYaxis()->SetTitle("Lab Energy (MeV)");
  hKin->GetXaxis()->CenterTitle(); hKin->GetYaxis()->CenterTitle(); hKin->GetXaxis()->SetTitleSize(0.05); hKin->GetYaxis()->SetTitleSize(0.05);
*/
  TCanvas *cExc_v_Ang = new TCanvas("cExc_v_Ang");
  TH2F *hExc_v_Ang = new TH2F("hExc_v_Ang","Exciation v Angle",180,0,180,2000,-5,15);
  tsingles->Draw("Exc:ThetaLab>>hExc_v_Ang",barrel_protons,"colz");
  hExc_v_Ang->GetXaxis()->SetTitle("Lab Angle (degrees)");  hExc_v_Ang->GetXaxis()->SetRangeUser(50,90);
  hExc_v_Ang->GetYaxis()->SetTitle("Excitation Energy (MeV)"); hExc_v_Ang->GetYaxis()->SetRangeUser(0.3,12.5);
  hExc_v_Ang->GetXaxis()->SetTitleSize(0.05);  hExc_v_Ang->GetYaxis()->SetTitleSize(0.05);
  hExc_v_Ang->GetXaxis()->CenterTitle();  hExc_v_Ang->GetYaxis()->CenterTitle();
  //TH2F *hExc_v_Ang_Coinc = new TH2F("hExc_v_Ang_Coinc","Exciation v Angle Coincidences",180,0,180,2000,-5,15);
  //tcoinc->Draw("Exc:ThetaLab>>hExc_v_Ang_Coinc",barrel_protons + Egam778keV,"colz");
/*
  TCanvas *cExc_EndCap = new TCanvas("cExc_EndCap");
  TH1F *hExc = new TH1F("hExc","Exciation EndCap",2000,-5,15);
  tsingles->Draw("Exc>>hExc",endcap_right);
  //TH1F *hExc_Coinc = new TH1F("hExc_Coinc","Exciation v Angle Coincidences",2000,-5,15);
  //tcoinc->Draw("Exc>>hExc_Coinc",endcap_left + Egam778keV,"same");
  TH1F *hExc_Coinc_ECR = new TH1F("hExc_Coinc_ECR","Exciation v Angle Coincidences",2000,-5,15);
  tcoinc->Draw("Exc>>hExc_Coinc_ECR",endcap_right + Egam778keV,"same");
  TH1F *hExc_Coinc_ECR2 = new TH1F("hExc_Coinc_ECR2","Exciation v Angle Coincidences",2000,-5,15);
  tcoinc->Draw("Exc>>hExc_Coinc_ECR2",endcap_right + Egam778keV + DeltaTS_Cut,"same");
  
  
  //hExc_Coinc->Scale(1/(0.082*0.52));
  //hExc_Coinc->SetLineColor(2);
  hExc_Coinc_ECR->SetLineColor(2);
  hExc_Coinc_ECR2->SetLineColor(6);

  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(hExc,"Singles","l");
  //legend->AddEntry(hExc_Coinc,"Coincidences Left","l");
  legend->AddEntry(hExc_Coinc_ECR,"Coincidences Right","l");
  legend->AddEntry(hExc_Coinc_ECR2,"Coincidences Right (deltaTS BG)","l");
  legend->Draw();

  TCanvas *cExc_Barrel = new TCanvas("cExc_Barrel");
  TH1F *hExc2 = new TH1F("hExc2","Exciation Barrel",2000,-5,15);
  tsingles->Draw("Exc>>hExc2",barrel + DividedTrigger);
  TH1F *hExc_Coinc2 = new TH1F("hExc_Coinc2","Exciation v Angle Coincidences",2000,-5,15);
  tcoinc->Draw("Exc>>hExc_Coinc2",barrel + Egam778keV,"same");
  //hExc_Coinc->Scale(1/(0.082*0.52));
  hExc_Coinc2->SetLineColor(2);

  TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
  legend2->AddEntry(hExc2,"Singles","l");
  legend2->AddEntry(hExc_Coinc2,"Coincidences","l");
  legend2->Draw();

  TCanvas *cEgam = new TCanvas("cEgam");
  TH2F *hEgam = new TH2F("hEgam","Egam",2000,0,200000000000000,10000,0,10000);
  tcoinc->Draw("Egam:deltaTS>>hEgam","","colz");

  TCanvas *cExc_v_Ang_O16 = new TCanvas("cExc_v_Ang_O16");
  TH2F *hExc_v_Ang_O16 = new TH2F("hExc_v_Ang_O16","Exciation v Angle Oxygen-16",180,0,180,2000,-5,15);
  tsingles->Draw("Exc_O16:ThetaLab>>hExc_v_Ang_O16",graphical_cut + barrel_protons,"colz");

  TCanvas *cExc_v_Ang_Coinc = new TCanvas("cExc_v_Ang_Coinc");
  TH2F *hExc_v_Ang_Coinc = new TH2F("hExc_v_Ang_Coinc","Exciation v Angle Coincidences",180,0,180,2000,-5,15);
  tcoinc->Draw("Exc:ThetaLab>>hExc_v_Ang_Coinc",graphical_cut + barrel_protons + Egam778keV,"colz");

  TCanvas *cExc_v_Ang_O16_Coinc = new TCanvas("cExc_v_Ang_O16_Coinc");
  TH2F *hExc_v_Ang_O16_Coinc = new TH2F("hExc_v_Ang_O16_Coinc","Exciation v Angle Ox
  ygen-16 Coincidences",180,0,180,2000,-5,15);
  tcoinc->Draw("Exc_O16:ThetaLab>>hExc_v_Ang_O16_Coinc",graphical_cut + barrel_protons + Egam778keV,"colz");

// Exciation Energy all protons
  TCanvas *c2 = new TCanvas("c2");
  TH1F *hExc_Singles = new TH1F("hExc_Singles","Exciation Energy",2000,-5,15);
  tsingles->Draw("Exc>>hExc_Singles",Cut_dE_E_EndCap + Cut_range_EndCap + EndCap_Left + angle_cutoff2);
  TH1F *hExc_Coinc = new TH1F("hExc_Coinc","Exciation Energy",2000,-5,15);
  tcoinc->Draw("Exc>>hExc_Coinc",Cut_dE_E_EndCap + Cut_range_EndCap + EndCap_Left + angle_cutoff2 + Egam778keV,"goff");
  hExc_Coinc->Scale(1/(0.083*0.55));
  hExc_Coinc->SetLineColor(2);
  hExc_Coinc->Draw("HIST+same");
*/
/// Particle - Gamm Matrix
TCanvas *cEgam_v_Exc = new TCanvas("cEgam_v_Exc"); cEgam_v_Exc->SetTicky(); cEgam_v_Exc->SetTickx(); cEgam_v_Exc->SetLogz(); 
TH2F *hEgam_v_Exc = new TH2F("hEgam_v_Exc","Egam vs Exc",10000,0,10000,2000,-5,15);
tcoinc->Draw("Exc:Egam>>hEgam_v_Exc",barrel_protons,"colz");
hEgam_v_Exc->GetXaxis()->SetTitle("#gamma-ray Energy (keV)"); hEgam_v_Exc->GetXaxis()->SetRangeUser(650,1000);
hEgam_v_Exc->GetYaxis()->SetTitle("Excitation Energy (MeV)"); hEgam_v_Exc->GetYaxis()->SetRangeUser(0.3,12.5);
hEgam_v_Exc->GetXaxis()->SetTitleSize(0.05); hEgam_v_Exc->GetYaxis()->SetTitleSize(0.05);
hEgam_v_Exc->GetXaxis()->CenterTitle(); hEgam_v_Exc->GetYaxis()->CenterTitle();
TLine *l1 = new TLine(650,9.15,1000,9.15); l1->SetLineColor(2); l1->SetLineWidth(2); l1->SetLineStyle(9); l1->Draw();

/*
    // Crystal vs Egam
    TCanvas *cCrystal = new TCanvas("cCrystal");
    TH2F *hCrys_v_Egam = new TH2F("hCrys_v_Egam","Crys v Egam",10000,0,10000,50,0,50);
    tcoinc->Draw("Crys:Egam>>hCrys_v_Egam",barrel_protons + angle_cutoff + avoid_contaminant,"colz");

    TCanvas *c3 = new TCanvas("c3");
    TH1F *hExc_C12_DS_Q5 = new TH1F("hExc_C12_DS_Q5","Downstream QQQ5 C12 Exciation Energy",2000,-5,15);
    tcoinc->Draw("Exc_C12>>hExc_C12_DS_Q5",endcap_protons);
    TH1F *hExc_C12_DS_SX3 = new TH1F("hExc_C12_DS_SX3","Downstream SX3 C12 Exciation Energy",2000,-5,15);
    tcoinc->Draw("Exc_C12>>hExc_C12_DS_SX3",barrel_protons,"same");
    hExc_C12_DS_SX3->SetLineColor(2);

    TCanvas *c4 = new TCanvas("c4");
    hExc_C12_DS_Q5->Draw();
    hExc_C12_US_Q5->Draw("same");
    hExc_C12_US_Q5->SetLineColor(2);
*/
}