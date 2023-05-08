Double_t CrystalBall(Double_t *x, Double_t *par) {

    return par[4]*ROOT::Math::crystalball_function(x[0], par[0], par[1], par[2], par[3]);
}

void mo95gamma_prob() {

    gROOT->SetBatch(true);

    TFile *infile = TFile::Open("./output.root");
    TFile *OutFile = new TFile("Prob_Gamma.root", "recreate");

    int bin;
    double sig = 1.8;
    double amp = 200.;
    double bg = 10.;
    double BackGround;

    double step = 0.05;

    char fname[32];
    char hname[32];
    char htitle[32];

    int singles_bin;
    const int n1 = 30;

    TH2D *hExc_Gam = (TH2D*)infile->Get("hExc_v_Egam");
    TH1D *hExc = (TH1D*)infile->Get("hExc_Subtracted_Sum");
    TH1D *hProj[n1];

    hExc->Rebin(5);    
    hExc->Write();

    hExc_Gam->RebinY(5);
    hExc_Gam->Write();

    step = 0.05;

    double Counts[n1]; 
    double ExMid[n1]; 
    double Counts_Error[n1];
    double Ex_Error[n1];
    double Prob[n1];
    double Prob_Error[n1];
    double Singles, Singles_Error;

    for(int i=0; i<n1; i++) {

        bin = i + 101 + 182;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts[i] = hProj[i]->Integral(201,208) - hProj[i]->Integral(197,201) - hProj[i]->Integral(208,211);
        Counts_Error[i] = sqrt(Counts[i]);

        ExMid[i] = 9.1 + i*step + step/2.;
        Ex_Error[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob[i] = Counts[i] / (Singles*0.166);
        Prob_Error[i] = Prob[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error[i]/Counts[i],2));

    }

    TGraphErrors *prob_204keV = new TGraphErrors(n1,ExMid,Prob,Ex_Error,Prob_Error);
    prob_204keV->SetMarkerStyle(20);
    prob_204keV->SetTitle("Coincidence Probability");
    prob_204keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_204keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_204keV->GetYaxis()->SetRangeUser(0.0001,1);
    prob_204keV->GetXaxis()->SetRangeUser(7,12);
    prob_204keV->SetLineColor(1);
    prob_204keV->Write();

    const int n2 = 24;
    double Counts2[n2];
    double Counts_Error2[n2];
    double ExMid2[n2];
    double Ex_Error2[n2];
    double Prob2[n2];
    double Prob_Error2[n2];   

    for(int i=0; i<n2; i++) {

        bin = i + 101 + 188;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts2[i] = hProj[i]->Integral(762,770) - hProj[i]->Integral(754,762); //(10./6.)*(hProj[i]->Integral(785,786) - hProj[i]->Integral(770,773));
        Counts_Error2[i] = sqrt(Counts2[i]);

        ExMid2[i] = 9.4 + i*step + step/2.;
        Ex_Error2[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob2[i] = Counts2[i] / (Singles*0.0838);
        Prob_Error2[i] = Prob2[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error2[i]/Counts2[i],2));

    }

    TGraphErrors *prob_766keV = new TGraphErrors(n2,ExMid2,Prob2,Ex_Error2,Prob_Error2);
    prob_766keV->SetTitle("Coincidence Probability");
    prob_766keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_766keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_766keV->GetYaxis()->SetRangeUser(0.001,10);
    prob_766keV->GetXaxis()->SetRangeUser(7,12);
    prob_766keV->SetLineColor(2);
    prob_766keV->SetMarkerColor(2);
    prob_766keV->SetMarkerStyle(21);
    prob_766keV->Write();


    const int n3 = 27;
    double Counts3[n3];
    double Counts_Error3[n3];
    double ExMid3[n3];
    double Ex_Error3[n3];
    double Prob3[n3];
    double Prob_Error3[n3];

    for(int i=0; i<n3; i++) {

        bin = i + 101 + 190;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts3[i] = hProj[i]->Integral(944,952) - hProj[i]->Integral(940,944) - hProj[i]->Integral(952,956);
        Counts_Error3[i] = sqrt(Counts3[i]);

        ExMid3[i] = 9.5 + i*step + step/2.;
        Ex_Error3[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob3[i] = Counts3[i] / (Singles*0.0746);
        Prob_Error3[i] = Prob3[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error3[i]/Counts3[i],2));

    }

    TGraphErrors *prob_948keV = new TGraphErrors(n3,ExMid3,Prob3,Ex_Error3,Prob_Error3);
    prob_948keV->SetTitle("Coincidence Probability");
    prob_948keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_948keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_948keV->GetYaxis()->SetRangeUser(0.001,1);
    prob_948keV->GetXaxis()->SetRangeUser(7,12);
    prob_948keV->SetLineColor(4);
    prob_948keV->SetMarkerColor(4);
    prob_948keV->SetMarkerStyle(22);
    prob_948keV->Write();
/*
    const int n4 = 110;
    double Counts4[n4];
    double Counts_Error4[n4];
    double ExMid4[n4];
    double Ex_Error4[n4];
    double Prob4[n4];
    double Prob_Error4[n4];


    for(int i=0; i<n4; i++) {

        bin = i + 51;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts4[i] = hProj[i]->Integral(1088,1094) - hProj[i]->Integral(1085,1088) - hProj[i]->Integral(1094,1097);
        Counts_Error4[i] = sqrt(Counts4[i]);

        ExMid4[i] = i*step + step/2.;
        Ex_Error4[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob4[i] = Counts4[i] / (Singles*0.0692);
        Prob_Error4[i] = Prob4[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error4[i]/Counts4[i],2));

    }

    TGraphErrors *prob_1091keV = new TGraphErrors(n4,ExMid4,Prob4,Ex_Error4,Prob_Error4);
    prob_1091keV->SetTitle("Coincidence Probability");
    prob_1091keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_1091keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_1091keV->GetYaxis()->SetRangeUser(0.001,1);
    prob_1091keV->GetXaxis()->SetRangeUser(7,12);
    prob_1091keV->SetLineColor(3);
    prob_1091keV->SetMarkerColor(3);
    prob_1091keV->SetMarkerStyle(23);
    prob_1091keV->Write();

    const int n5 = 105;
    double Counts5[n5];
    double Counts_Error5[n5];
    double ExMid5[n5];
    double Ex_Error5[n5];
    double Prob5[n5];
    double Prob_Error5[n5];

    for(int i=0; i<n5; i++) {

        bin = i + 51;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts5[i] = hProj[i]->Integral(808,816) - hProj[i]->Integral(800,808);
        Counts_Error5[i] = sqrt(Counts5[i]);

        ExMid5[i] = i*step + step/2.;
        Ex_Error5[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob5[i] = Counts5[i] / (Singles*0.081);
        Prob_Error5[i] = Prob5[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error5[i]/Counts5[i],2));

    }

    TGraphErrors *prob_812keV = new TGraphErrors(n5,ExMid5,Prob5,Ex_Error5,Prob_Error5);
    prob_812keV->SetTitle("Coincidence Probability");
    prob_812keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_812keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_812keV->GetYaxis()->SetRangeUser(0.001,1);
    prob_812keV->GetXaxis()->SetRangeUser(7,12);
    prob_812keV->SetLineColor(6);
    prob_812keV->SetMarkerColor(6);
    prob_812keV->SetMarkerStyle(24);
    prob_812keV->Write();

*/
    gROOT->SetBatch(false);
    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    prob_204keV->Draw("AP");
    prob_766keV->Draw("Psame");
    prob_948keV->Draw("Psame");
    //prob_1091keV->Draw("Psame");
    //prob_812keV->Draw("Psame");

    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetBorderSize(0);
   	legend->AddEntry(prob_204keV,"^{95}Mo 3/2_{1}^{+} #rightarrow 5/2_{g.s.}^{+} 204 keV","lep");
	legend->AddEntry(prob_766keV,"^{95}Mo 7/2_{1}^{+} #rightarrow 5/2_{g.s.}^{+} 765 keV","lep");
    legend->AddEntry(prob_948keV,"^{96}Mo 9/2_{1}^{+} #rightarrow 5/2_{g.s.}^{+} 948 keV","lep");
    //legend->AddEntry(prob_1091keV,"^{96}Mo 4_{2}^{+} #rightarrow 2_{1}^{+} 1091 keV","lep");
    //legend->AddEntry(prob_812keV,"^{96}Mo 6_{1}^{+} #rightarrow 4_{1}^{+} 812 keV","lep");
   	legend->Draw();
}