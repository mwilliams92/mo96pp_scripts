{  
  TFile *fdata = TFile::Open("../root_outputs/SumRun_13-44.root", "READ");
  //TFile *fdata = TFile::Open("../root_outputs/Run0036_sorted.root", "READ");
  //TFile *fdata = TFile::Open("../root_outputs/SumRun_11-29.root", "READ");
  TTree *tsingles = (TTree*)fdata->Get("t1");
  TTree *tcoinc = (TTree*)fdata->Get("t2");

  //tsingles->SetAlias("Etot_Cor","Etot*1.0167");

  TCut barrel_protons = "Range>2300 && Range<2900 && EndCap==0 && Upstream==0";
  //TCut endcap_protons_right = "Range>2800 && Range<4200 && EndCap==1 && Upstream==0 && Left==0";
  TCut endcap_protons_left = "EndCap==1 && Upstream==0 && Left==1";
  //TCut all_protons = "(Range>2400 && Range<3000 && EndCap==0) || (Range>2800 && Range<4200 && EndCap==1)";
  //TCut EndCap_Left = "EndCap==1 && Upstream==0 && Left==1";

  TCut angle_cutoff = "ThetaLab>65 && ThetaLab<88";
  TCut angle_cutoff2 = "ThetaLab>32 && ThetaLab<41";
  TCut angle_lol = "ThetaLab>77 && ThetaLab<87";
  TCut avoid_contaminant = "!(ThetaLab>45 && ThetaLab<65)";

  TCut CarbonState = "Egam>4200 && Egam<4600";
  TCut Egam778keV = "Egam>770 && Egam<784";

  TCutG *cutg = new TCutG("cutg",9);
  cutg->SetVarX("Etot");
  cutg->SetVarY("Range");
  cutg->SetPoint(0,8557.31,2285.07);
  cutg->SetPoint(1,6361.97,2235.42);
  cutg->SetPoint(2,4020.28,2299.26);
  cutg->SetPoint(3,2640.35,2405.64);
  cutg->SetPoint(4,3225.78,2802.82);
  cutg->SetPoint(5,4689.33,2923.39);
  cutg->SetPoint(6,6403.79,2923.39);
  cutg->SetPoint(7,8850.02,2895.02);
  cutg->SetPoint(8,8557.31,2285.07);

  TCutG *dE_E_Barrel = new TCutG("dE_E_Barrel",10);
  dE_E_Barrel->SetVarX("Etot");
  dE_E_Barrel->SetVarY("DeltaE");
  dE_E_Barrel->SetPoint(0,2904.93,2069.55);
  dE_E_Barrel->SetPoint(1,4611.02,1510.98);
  dE_E_Barrel->SetPoint(2,5848.77,1264.83);
  dE_E_Barrel->SetPoint(3,7153.42,1084.95);
  dE_E_Barrel->SetPoint(4,9495.11,924.008);
  dE_E_Barrel->SetPoint(5,8692.25,649.457);
  dE_E_Barrel->SetPoint(6,6818.9,744.13);
  dE_E_Barrel->SetPoint(7,4008.87,1113.35);
  dE_E_Barrel->SetPoint(8,2470.04,1728.73);
  dE_E_Barrel->SetPoint(9,2904.93,2069.55);

  TCut graphical_cut = "cutg";
  TCut gCut_dE_E_Barrel = "dE_E_Barrel";

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

  TCutG *range_EndCap = new TCutG("range_EndCap",13);
  range_EndCap->SetVarX("Etot");
  range_EndCap->SetVarY("Range");
  range_EndCap->SetPoint(0,3123.48,2988.38);
  range_EndCap->SetPoint(1,4148.7,3878.17);
  range_EndCap->SetPoint(2,11211.3,3817.96);
  range_EndCap->SetPoint(3,15061.6,4226.06);
  range_EndCap->SetPoint(4,19025.7,4379.93);
  range_EndCap->SetPoint(5,21509,4687.68);
  range_EndCap->SetPoint(6,21486.3,3677.46);
  range_EndCap->SetPoint(7,15585.6,2740.85);
  range_EndCap->SetPoint(8,14674.3,3215.85);
  range_EndCap->SetPoint(9,11780.9,3015.14);
  range_EndCap->SetPoint(10,9570.96,2908.1);
  range_EndCap->SetPoint(11,4991.65,2914.79);
  range_EndCap->SetPoint(12,3123.48,2988.38);
  TCut Cut_range_EndCap = "range_EndCap";
/*
  /// Particle-ID
  TCanvas *cPID_Barrel = new TCanvas("cPID_Barrel");
  TH2F *hPID_Barrel = new TH2F("hPID_Barrel","Barrel Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("DeltaE:Etot>>hPID_Barrel","EndCap==0 && Upstream==0","colz");
  TH2F *hPID_Barrel_Coinc = new TH2F("hPID_Barrel_Coinc","Barrel Particle ID",2500,0,25000,2000,0,20000);
  tcoinc->Draw("DeltaE:Etot>>hPID_Barrel_Coinc","EndCap==0 && Upstream==0 && Egam>770 && Egam<784","same");
  hPID_Barrel_Coinc->SetMarkerStyle(20); hPID_Barrel_Coinc->SetMarkerColor(2); hPID_Barrel_Coinc->SetMarkerSize(0.1);
  hPID_Barrel->GetXaxis()->SetTitle("BB10 + SX3 Energy (keV)");
  hPID_Barrel->GetYaxis()->SetTitle("BB10 Energy (keV)");
  dE_E_Barrel->SetLineColor(1); dE_E_Barrel->SetLineWidth(2);
  dE_E_Barrel->Draw("same");

  TCanvas *cPID_Barrel2 = new TCanvas("cPID_Barrel2");
  TH2F *hPID_Barrel2 = new TH2F("hPID_Barrel2","Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("Range:Etot>>hPID_Barrel2","EndCap==0 && Upstream==0","colz");
  TH2F *hPID_Barrel2_Coinc = new TH2F("hPID_Barrel2_Coinc","Barrel Particle ID",2500,0,25000,2000,0,20000);
  tcoinc->Draw("Range:Etot>>hPID_Barrel2_Coinc","EndCap==0 && Upstream==0 && Egam>770 && Egam<784","same");
  hPID_Barrel2_Coinc->SetMarkerStyle(20); hPID_Barrel2_Coinc->SetMarkerColor(2); hPID_Barrel2_Coinc->SetMarkerSize(0.1);
  hPID_Barrel2->GetXaxis()->SetTitle("BB10 + SX3 Energy (keV)");
  hPID_Barrel2->GetYaxis()->SetTitle("BB10 Range (keV)");
  cutg->SetLineColor(1); cutg->SetLineWidth(2);
  cutg->Draw("same");

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
*/
  TCanvas *cPID_EndCap2 = new TCanvas("cPID_EndCap2");
  TH2F *hPID_EndCap2 = new TH2F("hPID_EndCap2","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tsingles->Draw("Range:Etot>>hPID_EndCap2","EndCap==1 && Upstream==0 && Left==1 && ThetaLab>27 && ThetaLab<41","colz");
  TH2F *hPID_EndCap2_Coinc = new TH2F("hPID_EndCap2_Coinc","EndCap Particle ID",2500,0,25000,2000,0,20000);
  tcoinc->Draw("Range:Etot>>hPID_EndCap2_Coinc","EndCap==1 && Upstream==0 && Left==1 && ThetaLab>27 && ThetaLab<41 && Egam>770 && Egam<784","same");
  hPID_EndCap2_Coinc->SetMarkerStyle(20); hPID_EndCap2_Coinc->SetMarkerColor(2); hPID_EndCap2_Coinc->SetMarkerSize(0.1);
  hPID_EndCap2->GetXaxis()->SetTitle("#DeltaE + E1 + E2 Energy (keV)");
  hPID_EndCap2->GetYaxis()->SetTitle("Range (keV)");
  range_EndCap->SetLineColor(1); range_EndCap->SetLineWidth(2);
  range_EndCap->Draw("same");

  /// Kinematics
  TCanvas *cKin = new TCanvas("cKin");
  TH2F *hKin = new TH2F("hKin","Kinematics",180,0,180,2500,0,25000);
  tsingles->Draw("Etot:ThetaLab>>hKin","Left==1","colz");
/*
  TCanvas *cExc_v_Ang = new TCanvas("cExc_v_Ang");
  TH2F *hExc_v_Ang = new TH2F("hExc_v_Ang","Exciation v Angle",180,0,180,2000,0,20);
  tsingles->Draw("Exc:ThetaLab>>hExc_v_Ang",barrel_protons,"colz");
  TCanvas *cExc_v_Ang2 = new TCanvas("cExc_v_Ang2");
  TH2F *hExc_v_Ang_Coinc = new TH2F("hExc_v_Ang_Coinc","Exciation v Angle Coincidences",180,0,180,2000,0,20);
  tcoinc->Draw("Exc:ThetaLab>>hExc_v_Ang_Coinc",barrel_protons + Egam778keV,"colz");
*/
  TCanvas *cExc = new TCanvas("cExc");
  TH1F *hExc = new TH1F("hExc","Exciation",500,0,20);
  tsingles->Draw("Exc>>hExc",endcap_protons_left + Cut_range_EndCap,"goff");
  TH1F *hExc_Coinc = new TH1F("hExc_Coinc","Exciation v Angle Coincidences",500,0,20);
  tcoinc->Draw("Exc>>hExc_Coinc",endcap_protons_left + Cut_range_EndCap + Egam778keV,"goff");
  hExc->Scale(0.082); // 0.52
  hExc->Draw();
  hExc_Coinc->Draw("HIST+same");
  hExc_Coinc->SetLineColor(2);
/*
  TCanvas *cExc_v_Ang_O16 = new TCanvas("cExc_v_Ang_O16");
  TH2F *hExc_v_Ang_O16 = new TH2F("hExc_v_Ang_O16","Exciation v Angle Oxygen-16",180,0,180,2000,-5,15);
  tsingles->Draw("Exc_O16:ThetaLab>>hExc_v_Ang_O16",graphical_cut + barrel_protons,"colz");

  TCanvas *cExc_v_Ang_Coinc = new TCanvas("cExc_v_Ang_Coinc");
  TH2F *hExc_v_Ang_Coinc = new TH2F("hExc_v_Ang_Coinc","Exciation v Angle Coincidences",180,0,180,2000,-5,15);
  tcoinc->Draw("Exc:ThetaLab>>hExc_v_Ang_Coinc",graphical_cut + barrel_protons + Egam778keV,"colz");

  TCanvas *cExc_v_Ang_O16_Coinc = new TCanvas("cExc_v_Ang_O16_Coinc");
  TH2F *hExc_v_Ang_O16_Coinc = new TH2F("hExc_v_Ang_O16_Coinc","Exciation v Angle Oxygen-16 Coincidences",180,0,180,2000,-5,15);
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

/// Particle - Gamm Matrix
TCanvas *cEgam_v_Exc = new TCanvas("cEgam_v_Exc");
TH2F *hEgam_v_Exc = new TH2F("hEgam_v_Exc","Egam vs Exc",10000,0,10000,2000,-5,15);
tcoinc->Draw("Exc:Egam>>hEgam_v_Exc",graphical_cut + barrel_protons + angle_cutoff,"colz");

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