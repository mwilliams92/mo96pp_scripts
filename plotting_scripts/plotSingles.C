{
    // In File (all runs)
    //TFile *infile = TFile::Open("../root_outputs/Sum_Singles.root");
    //TTree *data = (TTree*)infile->Get("t1");

    // In File (one run)
    TFile *infile = TFile::Open("../root_outputs/Run0060_Singles.root");
    TTree *data = (TTree*)infile->Get("t1");

/*
    // For a subset of runs:
    TChain* data = new TChain ("t1");
    int firstRun = 60;
	int lastRun = 62;
	char fname[64];

    for (int i=firstRun; i<=lastRun; i++) {
	
		sprintf(fname,"../root_outputs/Run00%d_Singles.root", i);
		data->Add(fname);

	}
*/
    // Out File
    TFile *outfile = new TFile("SinglesPlots.root","RECREATE");
    // Directories:
    TDirectory *dirPID = outfile->mkdir("ParticleID");
    TDirectory *dirKin = outfile->mkdir("Kinematics");
    TDirectory *dirExc1D = outfile->mkdir("Excitation_1D");
    TDirectory *dirExc2D = outfile->mkdir("Excitation_2D");
    TDirectory *dirTiming = outfile->mkdir("Timing");

    // Cuts
    TCut ValidHit = "Energy>0";
    TCut Protons = "Range_BB10>2200 && Range_BB10<3000";
    TCut NoPunchThru = "Energy<8000";
    TCut GRETINA_TDC = "tdcGRETINA>1500 && tdcGRETINA<3000";
    TCut Silicon_GRETINA_TDC = "tdcSilicon_GRETINA>750 && tdcSilicon_GRETINA<850";

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Particle ID Plots

    dirPID->cd();
    TCanvas *cPID_Barrel = new TCanvas("cPID_Barrel");
    TH2D *hPID_Barrel = new TH2D("hPID_Barrel","Particle ID in Barrel",2000,0,20000,2000,0,20000);
    data->Draw("Range_BB10:Energy>>hPID_Barrel",ValidHit,"colz");
    cPID_Barrel->Write();
    hPID_Barrel->Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Kinematic 2D Plots

    dirKin->cd();
    TCanvas *cKin = new TCanvas("cKin");
    TH2D *hKin = new TH2D("hKin","Energy v Lab Angle",180,0,180,2000,0,20000);
    data->Draw("Energy:Angle>>hKin",ValidHit,"colz");
    hKin->Write();
    cKin->Write();

    TCanvas *cKinGated = new TCanvas("cKinGated");
    TH2D *hKinGated = new TH2D("hKinGated","Energy v Lab Angle (Proton Gated)",180,0,180,2000,0,20000);
    data->Draw("Energy:Angle>>hKinGated",ValidHit + Protons,"colz");
    hKinGated->Write();
    cKinGated->Write();

    TCanvas *cKinGated2 = new TCanvas("cKinGated2");
    TH2D *hKinGated2 = new TH2D("hKinGated2","Energy v Lab Angle (Proton Gated + GRETINA TDC Hit)",180,0,180,2000,0,20000);
    data->Draw("Energy:Angle>>hKinGated2",ValidHit + Protons + GRETINA_TDC,"colz");
    hKinGated2->Write();
    cKinGated2->Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Excitation 1D Plots

    dirExc1D->cd();
    TCanvas *cExc = new TCanvas("cExc");
    TH1D *hExc = new TH1D("hExc","Proton Gated Excitation Energy",2000,-5,15);
    data->Draw("Excitation>>hExc",ValidHit + Protons);
    hExc->SetLineColor(1);
    TH1D *hExc_Gated = new TH1D("hExc_Gated","Excitation Energy Gated",2000,-5,15);
    data->Draw("Excitation>>hExc_Gated",ValidHit + Protons + GRETINA_TDC,"same");
    hExc_Gated->SetLineColor(2);
    TH1D *hExc_Gated2 = new TH1D("hExc_Gated2","Excitation Energy Gated",2000,-5,15);
    data->Draw("Excitation>>hExc_Gated2",ValidHit + Protons + Silicon_GRETINA_TDC,"same");
    hExc_Gated2->SetLineColor(4);

    TLegend *legend0 = new TLegend(0.1,0.7,0.48,0.9);
   	legend0->AddEntry(hExc,"No Time Gate","l");
	legend0->AddEntry(hExc_Gated,"1500 < TDC_GRETINA < 3000","l");
    legend0->AddEntry(hExc_Gated2,"750 < TDC_Silicon_GRETINA < 850","l");
   	legend0->Draw();

    hExc->Write();
    hExc_Gated->Write();
    hExc_Gated2->Write();
    cExc->Write();

    TCanvas *cExc16O = new TCanvas("cExc16O");
    TH1D *hExc16O = new TH1D("hExc16O","Proton Gated Excitation Energy 16O",2000,-5,15);
    data->Draw("Excitation_O16>>hExc16O",ValidHit + Protons);
    hExc16O->SetLineColor(1);
    TH1D *hExc_Gated16O = new TH1D("hExc_Gated16O","Excitation Energy Time Gated",2000,-5,15);
    data->Draw("Excitation_O16>>hExc_Gated16O",ValidHit + Protons + GRETINA_TDC,"same");
    hExc_Gated16O->SetLineColor(2);
    //TH1D *hExc_Gated2 = new TH1D("hExc_Gated2","Excitation Energy Gated",2000,-5,15);
    //data->Draw("Excitation>>hExc_Gated2",ValidHit + Protons + Silicon_GRETINA_TDC,"same");
    //hExc_Gated2->SetLineColor(4);

    TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9);
   	legend3->AddEntry(hExc16O,"No Time Gate","l");
	legend3->AddEntry(hExc_Gated16O,"1500 < TDC_GRETINA < 3000","l");
   	legend3->Draw();

    hExc16O->Write();
    hExc_Gated16O->Write();
    cExc16O->Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Excitation 2D Plots

    dirExc2D->cd();
    TCanvas *cExcAng = new TCanvas("cExcAng");
    TH2D *hExcAng = new TH2D("hExcAng","Excitation Energy v Lab Angle",180,0,180,2000,-5,15);
    data->Draw("Excitation:Angle>>hExcAng",ValidHit + Protons,"colz");
    hExcAng->Write();
    cExcAng->Write();

    TCanvas *cExcAngCoinc = new TCanvas("cExcAngCoinc");
    TH2D *hExcAngCoinc = new TH2D("hExcAngCoinc","Excitation Energy v Lab Angle GRETINA Gated",180,0,180,2000,-5,15);
    data->Draw("Excitation:Angle>>hExcAngCoinc",ValidHit + Protons + GRETINA_TDC,"colz");
    hExcAngCoinc->Write();
    cExcAngCoinc->Write();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Timing Plots

    dirTiming->cd();
    TH1D *hTDC_Silicon = new TH1D("hTDC_Silicon","TDC Silicon",4096,0,4096);
    TH1D *hTDC_GRETINA = new TH1D("hTDC_GRETINA","TDC GRETINA",4096,0,4096);
    TH1D *hTDC_RF = new TH1D("hTDC_RF","TDC RF",4096,0,4096);
    TH1D *hTDC_Silicon_Div = new TH1D("hTDC_Silicon_Div","TDC Silicon Divided",4096,0,4096);
    TH1D *hTDC_Silicon_GRETINA = new TH1D("hTDC_Silicon_GRETINA","TDC Silicon-GRETINA",4096,0,4096);
    TH1D *hTDC_Silicon_Delay = new TH1D("hTDC_Silicon_Delay","TDC Silicon Delayed",4096,0,4096);
    TH1D *hTDC_Silicon_Upstream = new TH1D("hTDC_Silicon_Upstream","TDC Silicon Upsteam",4096,0,4096);

    TCanvas *cTDC_Silicon = new TCanvas("cTDC_Silicon");
    data->Draw("tdcSilicon>>hTDC_Silicon");
    data->Draw("tdcSilicon_Div>>hTDC_Silicon_Div","","same");
    data->Draw("tdcSilicon_Delay>>hTDC_Silicon_Delay","","same");
    data->Draw("tdcSilicon_Upstream>>hTDC_Silicon_Upstream","","same");

    hTDC_Silicon->SetLineColor(1);
    hTDC_Silicon_Div->SetLineColor(2);
    hTDC_Silicon_Delay->SetLineColor(4);
    hTDC_Silicon_Upstream->SetLineColor(6);

    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   	legend->AddEntry(hTDC_Silicon,"TDC Silicon (CH#833)","l");
	legend->AddEntry(hTDC_Silicon_Div,"TDC Silicon Divided (CH#834)","l");
	legend->AddEntry(hTDC_Silicon_Delay,"TDC Silicon Delayed (CH#840)","l");
    legend->AddEntry(hTDC_Silicon_Upstream,"TDC Silicon Upsteam (CH#841))","l");
   	legend->Draw();

    cTDC_Silicon->Write();
    hTDC_Silicon->Write();
    hTDC_Silicon_Div->Write();
    hTDC_Silicon_Delay->Write();
    hTDC_Silicon_Upstream->Write();


    TCanvas *cTDC_GRETINA = new TCanvas("cTDC_GRETINA");
    hTDC_Silicon->Draw();
    data->Draw("tdcGRETINA>>hTDC_GRETINA","","same");
    data->Draw("tdcSilicon_GRETINA>>hTDC_Silicon_GRETINA","","same");

    hTDC_GRETINA->SetLineColor(2);
    hTDC_Silicon_GRETINA->SetLineColor(4);

    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
   	legend2->AddEntry(hTDC_Silicon,"TDC Silicon (CH#833)","l");
	legend2->AddEntry(hTDC_GRETINA,"TDC GRETINA (CH#839)","l");
	legend2->AddEntry(hTDC_Silicon_GRETINA,"TDC Silicon-GRETINA (CH#835)","l");
   	legend2->Draw();

    cTDC_GRETINA->Write();
    hTDC_GRETINA->Write();
    hTDC_Silicon_GRETINA->Write();

}
