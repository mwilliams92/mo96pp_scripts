{  
   //TFile *fdata = TFile::Open("../root_outputs/SumRun_59-93.root", "READ");
   TFile *fdata = TFile::Open("../root_outputs/Run0060_sorted.root", "READ");
   //TFile *fdata = TFile::Open("../root_outputs/SumRun_11-29.root", "READ");
   TTree *tsingles = (TTree*)fdata->Get("t1");
   TTree *tcoinc = (TTree*)fdata->Get("t2");

   tsingles->SetAlias("RangeMeV","Range*0.001");
   tsingles->SetAlias("EtotMeV","Etot*0.001");

   TCut barrel_protons = "Range>2400 && Range<3000 && EndCap==0 && Upstream==0";
   TCut barrel = "EndCap==0 && Upstream==0";
   TCut endcap = "EndCap==1 && Upstream==0";

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
   TCanvas *cDeltaTS = new TCanvas("cDeltaTS");
   TH1F *hTS = new TH1F("hTS","Time Stamp Difference",2000,-1000,1000);
   tcoinc->Draw("deltaTS>>hTS");

   TH1F *hTS_Gated = new TH1F("hTS_Gated","Time Stamp Difference",2000,-1000,1000);
   tcoinc->Draw("deltaTS>>hTS_Gated","tdcSilicon_GRETINA>810 && tdcSilicon_GRETINA<840","same");
   hTS_Gated->SetLineColor(2);

   TH1F *hTS_Gated2 = new TH1F("hTS_Gated2","Time Stamp Difference",2000,-1000,1000);
   tcoinc->Draw("deltaTS>>hTS_Gated2","tdcSilicon_Div>390 && tdcSilicon_Div<420","same");
   hTS_Gated2->SetLineColor(6);

   TCanvas *cTDC = new TCanvas("cTDC");
   TH1F *hTDC = new TH1F("hTDC","Time Stamp Difference",4096,0,4096);
   tcoinc->Draw("tdcSilicon_GRETINA>>hTDC");
   TH1F *hTDC_Div = new TH1F("hTDC_Div","Time Stamp Difference",4096,0,4096);
   tcoinc->Draw("tdcSilicon_Div>>hTDC_Div","","same");
   hTDC_Div->SetLineColor(2);

   TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(hTDC,"tdcSiliocn_GRETINA","l");
   legend->AddEntry(hTDC_Div,"tdcSilicon_Div","l");
   legend->Draw();

   TCanvas *cTDC2 = new TCanvas("cTDC2");
   TH1F *hTDC_s = new TH1F("hTDC_s","Time Stamp Difference",4096,0,4096);
   tsingles->Draw("tdcSilicon_GRETINA>>hTDC_s");
   TH1F *hTDC_Div_s = new TH1F("hTDC_Div_s","Time Stamp Difference",4096,0,4096);
   tsingles->Draw("tdcSilicon_Div>>hTDC_Div_s","","same");
   hTDC_Div_s->SetLineColor(2);

   TLegend *legend_s = new TLegend(0.1,0.7,0.48,0.9);
   legend_s->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend_s->AddEntry(hTDC_s,"tdcSiliocn_GRETINA Singles","l");
   legend_s->AddEntry(hTDC_Div_s,"tdcSilicon_Div Singles","l");
   legend_s->Draw();
*/
   TCanvas *cDeltaTS_tdcGRETINA = new TCanvas("cDeltaTS_tdcGRETINA");
   TH2F *hDeltaTS_tdcGRETINA = new TH2F("hDeltaTS_tdcGRETINA_Silicon","Time Stamp Difference vs tdcGRETINA",4097,-1,4096,2000,-1000,1000);
   tcoinc->Draw("deltaTS:tdcGRETINA>>hDeltaTS_tdcGRETINA_Silicon","EndCap==1 && Upstream==0","colz");
   hDeltaTS_tdcGRETINA->GetXaxis()->SetTitle("TDC GRETINA");
   hDeltaTS_tdcGRETINA->GetYaxis()->SetTitle("ORRUBA - GRETINA Timestamp");

}