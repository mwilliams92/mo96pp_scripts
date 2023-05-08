Double_t CrystalBall(Double_t *x, Double_t *par) {

    return par[4]*ROOT::Math::crystalball_function(x[0], par[0], par[1], par[2], par[3]);
}

void gamma_prob() {

    gROOT->SetBatch(true);

    //TFile *infile = TFile::Open("./root_outputs/Sum_Coinc.root");
    TFile *infile = TFile::Open("./output.root");
    //TFile *infile2 = TFile::Open("./plotting_scripts/subtraction.root");
    //TTree *t1 = (TTree*)InFile->Get("t1");
    //TFile *InFile_s = TFile::Open("./root_outputs/SumRun_Singles.root");
    //TTree *t1_s = (TTree*)InFile_s->Get("t1");

    TFile *OutFile = new TFile("Prob_Gamma.root", "recreate");

    const int nsteps = 250;
    double Counts[nsteps];
   
    double ExMid[nsteps];
    
    double Counts_Error[nsteps];
    double Ex_Error[nsteps];
    double Prob[nsteps];
    double Prob_Error[nsteps];
    double Singles, Singles_Error;


    const int n2 = 100;
    double Counts2[n2];
    double Counts_Error2[n2];
    double ExMid2[n2];
    double Ex_Error2[n2];
    double Prob2[n2];
    double Prob_Error2[n2];

    int bin;

    TH1D *h[nsteps];
    TH1D *h2[nsteps];
    TF1 *fit[nsteps];
    TCanvas *c[nsteps];

    double sig = 1.8;
    double amp = 200.;
    double bg = 10.;
    double BackGround;

    double min = 6.0;
    double step = 0.05;
    double ExLow, ExHigh;

    char sExCut[32];
    char draw[32];
    char fname[32];
    char hname[32];
    char htitle[32];

    double mu[5] = {778., 849., 1091., 720., 766.};
    double range[5] = {6., 6., 6., 6., 6.};

    double EnergyStart[5] = {0.0, 9.0, 1.2, 1.2, 9.0};

    double Efficiency[5] = {0.083, 0.1, 0.1, 0.1, 0.1};

    int steps[5] = {240, 25, 90, 85, 30};

    int singles_bin;
    

    TH2D *hExc_Gam = (TH2D*)infile->Get("hExc_v_Egam");
    TH1D *hExc = (TH1D*)infile->Get("hExc_Subtracted_Sum");
    TH1D *hProj[nsteps];
    TGraphErrors* coinc[5];
    //TGraphErrors* prob[5];

    hExc->Rebin(5);    
    hExc->Write();

    hExc_Gam->RebinY(5);
    hExc_Gam->Write();

    step = 0.05;

    for(int i=0; i<250; i++) {

        bin = i + 101;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");
/*
        // Define Fit and Set Params
        sprintf(fname,"fit[%d]",i);
        fit[i] = new TF1(fname,"gaus(0)+pol0(3)",770,786);
        fit[i]->SetParameters(amp,778,sig,bg);
        fit[i]->SetParLimits(1,774,782);
        fit[i]->SetParLimits(2,1,2.5);
        hProj[i]->Fit(fname,"QR");
        fit[i]->Draw("same");
        hProj[i]->Write();

        BackGround = fit[i]->GetParameter(3); */
        //Counts[i] = (fit[i]->Integral(770,786) - BackGround*16);
        Counts[i] = hProj[i]->Integral(774,784) - hProj[i]->Integral(784,794); //(10./6.)*(hProj[i]->Integral(785,786) - hProj[i]->Integral(770,773));
        Counts_Error[i] = sqrt(Counts[i]);

        ExMid[i] = i*step + step/2.;
        Ex_Error[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob[i] = Counts[i] / (Singles*0.083*0.55);
        Prob_Error[i] = Prob[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error[i]/Counts[i],2));

    }

    TGraphErrors *prob_778keV = new TGraphErrors(250,ExMid,Prob,Ex_Error,Prob_Error);
    prob_778keV->SetMarkerStyle(20);
    prob_778keV->SetTitle("Coincidence Probability");
    prob_778keV->GetYaxis()->SetTitle("Couns / 50 keV");
    prob_778keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_778keV->GetYaxis()->SetRangeUser(0.0001,1);
    prob_778keV->GetXaxis()->SetRangeUser(7,12);
    prob_778keV->SetLineColor(1);
    prob_778keV->Write();

    step = 0.1;

    hExc_Gam->RebinY(2);
    hExc->Rebin(2);    

    for(int i=0; i<n2; i++) {

        bin = i + 51;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");
/*
        // Define Fit and Set Params
        sprintf(fname,"fit[%d]",i);
        fit[i] = new TF1(fname,"gaus(0)+pol0(3)",716,724);
        fit[i]->SetParameters(100,720,sig,bg);
        fit[i]->SetParLimits(1,716,724);
        fit[i]->SetParLimits(2,1,2.5);
        hProj[i]->Fit(fname,"QR");
        fit[i]->Draw("same");
        hProj[i]->Write();
*/
        //BackGround = fit[i]->GetParameter(3);
        //Counts[i] = (fit[i]->Integral(770,786) - BackGround*16);
        Counts2[i] = hProj[i]->Integral(716,724) - hProj[i]->Integral(708,716); //(10./6.)*(hProj[i]->Integral(785,786) - hProj[i]->Integral(770,773));
        Counts_Error2[i] = sqrt(Counts2[i]);

        ExMid2[i] = i*step + step/2.;
        Ex_Error2[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob2[i] = Counts2[i] / (Singles*0.0867*0.55);
        Prob_Error2[i] = Prob2[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error2[i]/Counts2[i],2));

    }

    TGraphErrors *prob_720keV = new TGraphErrors(n2,ExMid2,Prob2,Ex_Error2,Prob_Error2);
    prob_720keV->SetTitle("Coincidence Probability");
    prob_720keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_720keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_720keV->GetYaxis()->SetRangeUser(0.001,10);
    prob_720keV->GetXaxis()->SetRangeUser(7,12);
    prob_720keV->SetLineColor(2);
    prob_720keV->SetMarkerColor(2);
    prob_720keV->SetMarkerStyle(21);
    prob_720keV->Write();


    const int n3 = 110;
    double Counts3[n3];
    double Counts_Error3[n3];
    double ExMid3[n3];
    double Ex_Error3[n3];
    double Prob3[n3];
    double Prob_Error3[n3];

    for(int i=0; i<n3; i++) {

        bin = i + 51;
        sprintf(hname,"h[%d]",i);
        hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin,"goff");

        Counts3[i] = hProj[i]->Integral(844,854) - hProj[i]->Integral(839,844) - hProj[i]->Integral(854,859);
        Counts_Error3[i] = sqrt(Counts3[i]);

        ExMid3[i] = i*step + step/2.;
        Ex_Error3[i] = step/2;

        Singles = hExc->Integral(bin,bin);

        Prob3[i] = Counts3[i] / (Singles*0.0791);
        Prob_Error3[i] = Prob3[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error3[i]/Counts3[i],2));

    }

    TGraphErrors *prob_849keV = new TGraphErrors(n3,ExMid3,Prob3,Ex_Error3,Prob_Error3);
    prob_849keV->SetTitle("Coincidence Probability");
    prob_849keV->GetYaxis()->SetTitle("Couns / 100 keV");
    prob_849keV->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob_849keV->GetYaxis()->SetRangeUser(0.001,1);
    prob_849keV->GetXaxis()->SetRangeUser(7,12);
    prob_849keV->SetLineColor(4);
    prob_849keV->SetMarkerColor(4);
    prob_849keV->SetMarkerStyle(22);
    prob_849keV->Write();

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


    gROOT->SetBatch(false);
    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    prob_778keV->Draw("AP");
    prob_720keV->Draw("Psame");
    prob_849keV->Draw("Psame");
    prob_1091keV->Draw("Psame");
    prob_812keV->Draw("Psame");

    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetBorderSize(0);
   	legend->AddEntry(prob_778keV,"^{96}Mo 2_{1}^{+} #rightarrow 0_{g.s.}^{+} 778 keV","lep");
	legend->AddEntry(prob_720keV,"^{96}Mo 2_{2}^{+} #rightarrow 2_{1}^{+} 720 keV","lep");
    legend->AddEntry(prob_849keV,"^{96}Mo 4_{1}^{+} #rightarrow 4_{1}^{+} 849 keV","lep");
    legend->AddEntry(prob_1091keV,"^{96}Mo 4_{2}^{+} #rightarrow 2_{1}^{+} 1091 keV","lep");
    legend->AddEntry(prob_812keV,"^{96}Mo 6_{1}^{+} #rightarrow 4_{1}^{+} 812 keV","lep");
   	legend->Draw();




/*
    for (int j=0; j<1; j++) {

        for(int i=0; i<steps[j]; i++) {

            sprintf(hname,"h[%d]",i);
            //bin = 1000 + i*5;
            bin = 500 + EnergyStart[j]*100 + i*5;
            // 1 bin = 0.01 MeV. // 0 MeV = bin 500
            hProj[i] = hExc_Gam->ProjectionX(hname,bin,bin+5,"goff");
           

            // Define Fit and Set Params
            sprintf(fname,"fit[%d]",i);
            fit[i] = new TF1(fname,"gaus(0)+pol0(3)",mu[j]-range[j],mu[j]+range[j]);
            fit[i]->SetParameters(amp,mu[j],sig,bg);
            fit[i]->SetParLimits(1,0.98*mu[j],1.02*mu[j]);
            fit[i]->SetParLimits(2,0.5*sig,1.5*sig);
            hProj[i]->Fit(fname,"QR");
            fit[i]->Draw("same");

            hProj[i]->Write();
            //fit[i]->Write();

            BackGround = fit[i]->GetParameter(3);
            Counts[i] = (fit[i]->Integral(mu[j]-range[j],mu[j]+range[j]) - BackGround*range[j]); 
            Counts_Error[i] = sqrt(Counts[i]);
            ExMid[i] = bin*0.01 + 0.025 - 5.;
            Ex_Error[i] = step/2;

            //cout << i << ", " << Counts[i] << ", " << ExMid[i] << endl;
        
            singles_bin = 101 + i;

            Singles = hExc->GetBinContent(singles_bin);

            cout << ExMid[i] << ", " << singles_bin << ", " << Singles << endl;

            //Singles = hExc->Integral(bin,bin+10);
            //Singles = 4844*ExMid[i] - 23564; //2714.3
            Singles_Error = sqrt(Singles);

            Prob[i] = Counts[i] / (Singles*Efficiency[j]);
            Prob_Error[i] = Prob[i] * sqrt(pow(Singles_Error/Singles,2) + pow(Counts_Error[i]/Counts[i],2));

            //hProj[i]->Write();
            //fit[i]->Write();

        }

        coinc[j] = new TGraphErrors(steps[j],ExMid,Counts,Ex_Error,Counts_Error);
        coinc[j]->Draw("AP");
        coinc[j]->SetMarkerStyle(20);
        coinc[j]->Write();

        prob[j] = new TGraphErrors(steps[j],ExMid,Prob,Ex_Error,Prob_Error);
        prob[j]->Draw("AP");
        prob[j]->SetMarkerStyle(20);
        prob[j]->Write();

    }

    hExc_Gam->Write();

    gROOT->SetBatch(false);

    TCanvas *c1 = new TCanvas("c1");
    c1->SetLogy();
    coinc[0]->Draw("ALP");
    coinc[0]->SetMarkerStyle(20);
    coinc[0]->SetTitle("Coincidence Counts (no efficiency correction)");
    coinc[0]->GetYaxis()->SetTitle("Couns / 50 keV");
    coinc[0]->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    coinc[0]->GetYaxis()->SetRangeUser(10,10000);
    coinc[0]->SetLineColor(1);

    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    c2->SetGridy();
    c2->SetGridx();
    prob[0]->Draw("AP");
    prob[0]->SetMarkerStyle(20);
    prob[0]->SetTitle("Coincidence Probability");
    prob[0]->GetYaxis()->SetTitle("Couns / 50 keV");
    prob[0]->GetXaxis()->SetTitle("Excitation Energy (MeV)");
    prob[0]->GetYaxis()->SetRangeUser(0.001,10);
    prob[0]->GetXaxis()->SetRangeUser(7,12);
    prob[0]->SetLineColor(1);
*/
/*
    coinc[2]->Draw("LPsame");
    coinc[2]->SetMarkerStyle(21);
    coinc[2]->SetMarkerColor(2);
    coinc[2]->SetLineColor(2);

    coinc[3]->Draw("LPsame");
    coinc[3]->SetMarkerStyle(21);
    coinc[3]->SetMarkerColor(4);
    coinc[3]->SetLineColor(4);

    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
   	legend->AddEntry(coinc[0],"^{96}Mo 2_{1}^{+} #rightarrow 0_{g.s.}^{+} 778 keV","lep");
	legend->AddEntry(coinc[2],"^{96}Mo 4_{2}^{+} #rightarrow 2_{1}^{+} 1091 keV","lep");
    legend->AddEntry(coinc[3],"^{96}Mo 2_{2}^{+} #rightarrow 2_{1}^{+} 720 keV","lep");
   	legend->Draw();

    TCanvas *c2 = new TCanvas("c2");
    c2->SetLogy();
    coinc[1]->Draw("ALP");
    coinc[1]->SetMarkerStyle(20);
    coinc[1]->GetYaxis()->SetRangeUser(10,10000);

    coinc[1]->SetTitle("Coincidence Counts (no efficiency correction)");
    coinc[1]->GetYaxis()->SetTitle("Couns / 100 keV");
    coinc[1]->GetXaxis()->SetTitle("Excitation Energy (MeV)");

    coinc[4]->Draw("LPsame");
    coinc[4]->SetMarkerStyle(21);
    coinc[4]->SetMarkerColor(2);

    //coinc[3]->Draw("Psame");
    //coinc[3]->SetMarkerStyle(21);
    //coinc[3]->SetMarkerColor(4);

    TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);
   	legend2->AddEntry(coinc[0],"^{95}Mo 3/2_{1}^{+} #rightarrow 5/2_{g.s.}^{+} 204 keV","lep");
	legend2->AddEntry(coinc[2],"^{95}Mo 7/2_{1}^{+} #rightarrow 5/2_{g.s.}^{+} 766 keV","lep");
    //legend->AddEntry(coinc[3],"2_{2}^{+} #rightarrow 2_{1}^{+} 720 keV","lep");
   	legend2->Draw();

    c1->Write();
    c2->Write();
*/

}