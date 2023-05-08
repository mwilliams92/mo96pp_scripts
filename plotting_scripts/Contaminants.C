{


    TFile *mo96data = TFile::Open("../root_outputs/SumRun_Combined.root");
    TFile *contaminants = TFile::Open("./contaminants.root");

    TCanvas *c1 = new TCanvas();

    TH1D *hExc_O16 = (TH1D*)mo96data->Get("hExc_O16");
    //hExc_O16->Smooth(10);
    hExc_O16->Draw();
    hExc_O16->Smooth(1);
    TH1D *hExc_O16_6MeV = (TH1D*)contaminants->Get("hExc_O16_6MeV");
    hExc_O16_6MeV->Scale(2.2685);
    hExc_O16_6MeV->Draw("HIST+same");

    TF1 *f1 = new TF1("f1","gaus(0)+pol1(3)",5.5,6.7);
    f1->SetParameter(0,270);
    f1->SetParameter(1,6.15);
    f1->SetParameter(2,0.061);
    hExc_O16->Fit("f1","R");

    double A = f1->GetParameter(0);
    double mu = f1->GetParameter(1);
    double integral = A*mu*2*TMath::Pi();
    double A_er = f1->GetParError(0);
    double mu_er = f1->GetParError(1);
    double integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));
    
    cout << integral <<  " +/- " << integral_er << endl;

    TF1 *f2 = new TF1("f2","gaus",5.5,6.7);
    f2->SetParameter(0,270);
    f2->SetParameter(1,6.15);
    f2->SetParameter(2,0.061);
    hExc_O16_6MeV->Fit("f2","R");
    f2->Draw("same");

    A = f2->GetParameter(0);
    mu = f2->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f2->GetParError(0);
    mu_er = f2->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));
    
    cout << integral <<  " +/- " << integral_er << endl;

//////////////////////////////////////////////////////////////////

    TCanvas *c2 = new TCanvas();
    TH1D *hExc_O16_copy = (TH1D*) hExc_O16->Clone();
    hExc_O16_copy->Add(hExc_O16_6MeV,-1);
    hExc_O16_copy->Draw("HIST");
    hExc_O16_copy->Rebin(2);

    TF1 *f3 = new TF1("f3","gaus(0)+pol1(3)",7.05,7.7);
    f3->SetParameter(0,2500);
    f3->SetParameter(1,7.14);
    f3->FixParameter(2,0.077);
    f3->SetParameter(3,-2000);
    f3->SetParameter(4,500);
    hExc_O16_copy->Fit("f3","R");
    f3->Draw("same");

    A = f3->GetParameter(0);
    mu = f3->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f3->GetParError(0);
    mu_er = f3->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));

    cout << integral <<  " +/- " << integral_er << endl;

    TH1D *hExc_O16_7MeV = (TH1D*)contaminants->Get("hExc_O16_7MeV");
    hExc_O16_7MeV->Scale(1.31);
    hExc_O16_7MeV->Draw("HIST+same");
    hExc_O16_7MeV->Rebin(2);

    TF1 *f4 = new TF1("f4","gaus(0)",6.8,7.5);
    f4->SetParameter(0,2500);
    f4->SetParameter(1,7.14);
    f4->SetParameter(2,0.07);
    hExc_O16_7MeV->Fit("f4","R+");
    f4->Draw("same");

    A = f4->GetParameter(0);
    mu = f4->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f4->GetParError(0);
    mu_er = f4->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));
    
    cout << integral <<  " +/- " << integral_er << endl;

//////////////////////////////////////////////////////////////////////

    TCanvas *c3 = new TCanvas();
    TH1D *hExc_O16_copy2 = (TH1D*) hExc_O16_copy->Clone();
    hExc_O16_copy2->Add(hExc_O16_7MeV,-1);
    hExc_O16_copy2->Draw("HIST");
    //hExc_O16_copy2->Rebin(2);

    TF1 *f5 = new TF1("f5","gaus(0)+pol1(3)",6.5,7.1);
    f5->SetParameter(0,500);
    f5->SetParameter(1,6.95);
    f5->SetParameter(2,0.07);
    f5->SetParameter(3,-2000);
    f5->SetParameter(4,500);
    hExc_O16_copy2->Fit("f5","R+");
    f5->Draw("same");

    A = f5->GetParameter(0);
    mu = f5->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f5->GetParError(0);
    mu_er = f5->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));

    cout << integral <<  " +/- " << integral_er << endl;

    TH1D *hExc_O16_7MeV_2 = (TH1D*)contaminants->Get("hExc_O16_7MeV_2");
    hExc_O16_7MeV_2->Scale(0.312);
    hExc_O16_7MeV_2->Draw("HIST+same");
    hExc_O16_7MeV_2->Rebin(2);

    TF1 *f6 = new TF1("f6","gaus(0)",6.5,7.5);
    f6->SetParameter(0,2500);
    f6->SetParameter(1,6.96);
    f6->SetParameter(2,0.07);
    hExc_O16_7MeV_2->Fit("f6","R+");
    f6->Draw("same");

    A = f6->GetParameter(0);
    mu = f6->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f6->GetParError(0);
    mu_er = f6->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));
    
    cout << integral <<  " +/- " << integral_er << endl;

//////////////////////////////////////////////////////////////////////

    TCanvas *c4 = new TCanvas();
    TH1D *hExc_O16_copy3 = (TH1D*) hExc_O16_copy2->Clone();
    hExc_O16_copy3->Add(hExc_O16_7MeV_2,-1);
    hExc_O16_copy3->GetXaxis()->SetRangeUser(5.5,8.5);
    hExc_O16_copy3->Draw("HIST");
    
    TF1 *f7 = new TF1("f7","gaus(0)+pol1(3)",8.5,9.3);
    f7->SetParameter(0,500);
    f7->SetParameter(1,8.9);
    f7->SetParameter(2,0.07);
    f7->SetParameter(3,-2000);
    f7->SetParameter(4,500);
    hExc_O16_copy2->Fit("f7","R+");
    f7->Draw("same");

    A = f7->GetParameter(0);
    mu = f7->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f7->GetParError(0);
    mu_er = f7->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));

    cout << integral <<  " +/- " << integral_er << endl;

    TH1D *hExc_O16_9MeV = (TH1D*)contaminants->Get("hExc_O16_9MeV");
    hExc_O16_9MeV->Scale(0.439);
    hExc_O16_9MeV->Draw("HIST+same");
    hExc_O16_9MeV->Rebin(2);

    TF1 *f8 = new TF1("f8","gaus(0)",8.5,9.5);
    f8->SetParameter(0,2500);
    f8->SetParameter(1,6.96);
    f8->SetParameter(2,0.07);
    hExc_O16_9MeV->Fit("f8","R+");
    f8->Draw("same");

    A = f8->GetParameter(0);
    mu = f8->GetParameter(1);
    integral = A*mu*2*TMath::Pi();
    A_er = f8->GetParError(0);
    mu_er = f8->GetParError(1);
    integral_er = integral*sqrt(pow(A_er/A,2)+pow(mu_er/mu,2));
    
    cout << integral <<  " +/- " << integral_er << endl;

//////////////////////////////////////////////////////////////////////

    TCanvas *c5 = new TCanvas();
    hExc_O16->Draw("E0");
    hExc_O16->Rebin(2);
    TH1D *hExc_O16_copy4 = (TH1D*) hExc_O16_copy3->Clone();
    hExc_O16_copy4->Add(hExc_O16_9MeV,-1);
    //hExc_O16_copy4->GetXaxis()->SetRangeUser(5.5,8.5);
    hExc_O16_copy4->Draw("E0+same");
    hExc_O16_copy4->SetLineColor(2);

}
