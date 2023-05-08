{
    TFile *mo96data = TFile::Open("../output.root");
    TFile *subtraction = TFile::Open("subtraction.root","RECREATE");

    TH2D *hExc_Ang = (TH2D*)mo96data->Get("hExcAng");

    int start_ang = 66;
    int end_ang = 88;
    int angle;

    TH1D *hExc[23];
    TH1D *hExc_Subtracted[23];
    TH1D *hExc_Subtracted_Sum = new TH1D("hExc_Subtracted_Sum","Subtracted Excitation Energy",2000,-5,15);

    TH2D *hExc_v_Ang_Sub = new TH2D("hExc_v_Ang_Sub","Subtracted Excitation Energy v Angle",180,0,180,2000,-5,15);

    char hname[8], hname2[8];
    double f1_low, f1_mid, f1_high, f2_low, f2_mid, f2_mid2, f2_high, f3_low, f3_mid, f3_high, f4_low, f4_mid, f4_high;
    double fitval, subtract;
    double step = 0.02;
    double step2 = 0.015;
    double sigma = 0.07;

    //double start = 6.2

    hExc_Ang->Write();
    TH1D *hExc_Raw = (TH1D*) hExc_Ang->ProjectionY(hname,66,88);
    hExc_Raw->Write();
   

    for (int i=0; i<23; i++) {

        angle = start_ang+i;
        
        sprintf(hname,"hExc[%d]",i);

        hExc[i] = (TH1D*) hExc_Ang->ProjectionY(hname,angle,angle);
        //hExc[i]->Rebin(2);

        f1_mid = 6.84 + i*step;
        f1_low = f1_mid - 0.5;
        f1_high = f1_mid + 0.5;        

        TF1 *f1 = new TF1("f1","gaus(0)+pol1(3)",f1_low,f1_high);
        f1->SetParameter(0,800.);
        f1->SetParLimits(0,0,5000.);
        f1->SetParameter(1,f1_mid);
        f1->SetParLimits(1,f1_mid-1.0,f1_mid+1.0);
        f1->SetParameter(2,sigma);
        f1->SetParLimits(2,sigma*0.5,sigma*2.0);
        f1->SetParameter(3,50);
        f1->SetParameter(4,1);
        hExc[i]->Fit("f1","R+");

        f2_mid = 7.82 + i*step;
        f2_low = f2_mid - 0.5;
        f2_mid2 = f2_mid - 0.18;
        f2_high = f2_mid + 0.5;

        TF1 *f2 = new TF1("f2","gaus(0)+gaus(3)+pol1(6)",f2_low,f2_high);
        f2->SetParameter(0,2000.);
        f2->SetParLimits(0,0,5000.);
        f2->SetParameter(1,f2_mid);
        f2->SetParameter(2,sigma);
        f2->SetParLimits(2,sigma*0.8,sigma*1.2);
        f2->SetParameter(3,1000);
        f2->SetParLimits(3,0,2000);
        f2->SetParameter(4,f2_mid2);
        f2->SetParLimits(4,f2_mid2*0.9,f2_mid2*1.1);
        f2->SetParameter(5,sigma);
        f2->SetParLimits(5,sigma*0.5,sigma*2.0);
        f2->SetParameter(6,-50);
        f2->SetParameter(7,50);
        hExc[i]->Fit("f2","R+");

        f3_mid = 9.54 + i*step2;
        f3_low = f3_mid - 0.5;
        f3_high = f3_mid + 0.5;

        TF1 *f3 = new TF1("f3","gaus(0)+pol1(3)",f3_low,f3_high);
        f3->SetParameter(0,1000.);
        f3->SetParameter(1,f3_mid);
        f3->SetParLimits(1,f3_mid-0.1,f3_mid+0.1);
        f3->SetParameter(2,sigma);
        f3->SetParLimits(2,sigma*0.8,sigma*1.3);
        f3->SetParameter(3,-50);
        f3->SetParameter(4,50);
        hExc[i]->Fit("f3","R+");

        f4_mid = 10.6 + i*step;
        f4_low = f4_mid - 0.4;
        f4_high = f4_mid + 0.4;

        TF1 *f4 = new TF1("f4","gaus(0)+pol1(3)",f4_low,f4_high);
        f4->SetParameter(0,3000.);
        f4->SetParameter(1,f4_mid);
        f4->SetParameter(2,sigma);
        f4->SetParLimits(2,0.1*0.8,0.1*1.5);
        f4->SetParameter(3,-50);
        f4->SetParameter(4,50);
        hExc[i]->Fit("f4","R+");

        TF1 *f5 = new TF1("f5","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)");
        double par[15];
        par[0] = f1->GetParameter(0);
        par[1] = f1->GetParameter(1);
        par[2] = f1->GetParameter(2);
        par[3] = f2->GetParameter(0);
        par[4] = f2->GetParameter(1);
        par[5] = f2->GetParameter(2);
        par[6] = f2->GetParameter(3);
        par[7] = f2->GetParameter(4);
        par[8] = f2->GetParameter(5);
        par[9] = f3->GetParameter(0);
        par[10] = f3->GetParameter(1);
        par[11] = f3->GetParameter(2);
        par[12] = f4->GetParameter(0);
        par[13] = f4->GetParameter(1);
        par[14] = f4->GetParameter(2);

        for (int j=0; j<15; j++) {
            f5->FixParameter(j,par[j]);
        }

        f5->SetLineColor(1);
        f5->SetLineStyle(2);
        f5->SetNpx(1000);
        hExc[i]->Fit("f5");

        //cout << "Debug1" << endl;

        //TCanvas *c2 = new TCanvas("c2");
        //hExc_Subtracted[i] = (TH1D*) hExc[i]->Clone();
        sprintf(hname2,"hSub[%d]",i);
        hExc_Subtracted[i] = new TH1D(hname2,"Subtracted",2000,-5,15);
    
        for (int j=0; j<2000; j++) {

            fitval = f5->Eval(hExc[i]->GetBinCenter(j));
            subtract = hExc[i]->GetBinContent(j) - fitval;
            hExc_Subtracted[i]->SetBinContent(j,subtract);

            hExc_v_Ang_Sub->SetBinContent(i+start_ang,j,subtract);

        }
        hExc_Subtracted[i]->SetLineColor(2);

        hExc[i]->Write();
        hExc_Subtracted[i]->Write();

        hExc_Subtracted_Sum->Add(hExc_Subtracted[i]);

        //f1->Delete();
        f2->Delete();
        f3->Delete();
        f4->Delete();
        f5->Delete();
    
    }

    hExc_Subtracted_Sum->Write();
    hExc_v_Ang_Sub->Write();

}