{

    TFile *fdata = TFile::Open("../root_outputs/SumRun_59-93.root", "READ");
    TTree *tsingles = (TTree*)fdata->Get("t1");
    TTree *tcoinc = (TTree*)fdata->Get("t2");

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

    TCut barrel_protons = "cutg";
    TCut angle_cut = "ThetaLab>65 && ThetaLab<88";
    TCut barrel = "EndCap==0 && Upstream==0";
    TCut O16_State = "Exc_O16>8.75 && Exc_O16<9.05";

    TCanvas *c1 = new TCanvas("c1");
    tsingles->Draw("Exc_O16>>h1(2000,-5,15)",barrel + angle_cut);

    TF1 *f1 = new TF1("f1","gaus(0)+pol1(3)",8.6,9.2);
    f1->SetParameter(0,2000.);
    f1->SetParameter(1,8.835);
    f1->SetParameter(2,0.07);
    f1->SetParameter(3,500);
    f1->SetParameter(4,500);
    h1->Fit("f1","R+");

    double amp, sigma, bin;
    double constant = 0.3989;
    bin = 100.;            // bins per MeV
    amp = f1->GetParameter(0);
    sigma = f1->GetParameter(2);
    double singles = (bin * amp * sigma) / constant;  

    //cout << "Singles Counts = " << singles << endl;

    TCanvas *c2 = new TCanvas("c2");
    tcoinc->Draw("Exc_O16:Egam>>h2(10000,0,10000,2000,-5,15)",barrel + angle_cut,"colz");

    TCanvas *c3 = new TCanvas("c3");
    tcoinc->Draw("Egam>>h3(5000,0,5000)",barrel + angle_cut + O16_State);
    h3->Rebin(10);

    TF1 *f2 = new TF1("f2","gaus(0)+pol1(3)",2650,2850);
    f2->SetParameter(0,60.);
    f2->SetParLimits(0,20.,100.);
    f2->SetParameter(1,2734);
    f2->SetParLimits(1,2700.,2750.);
    f2->SetParameter(2,14);
    f2->SetParLimits(2,5.,20.);
    f2->SetParameter(3,100);
    f2->SetParameter(4,-0.03);
    h3->Fit("f2","R+");

    bin = 0.1;            // bins per keV
    amp = f2->GetParameter(0);
    sigma = f2->GetParameter(2);
    double coincidences = (bin * amp * sigma) / constant;

    double branching = 0.777;
    double efficiency = 0.04;
    double corrected_coincidences = coincidences / (branching * efficiency);


    double coinc_eff = corrected_coincidences / singles;

    cout << "O16 state singles counts = " <<  singles << endl;
    cout << "O16 state coincidence counts = " <<  coincidences << endl;
    cout << "O16 state corrected coincidence counts = " <<  corrected_coincidences << endl;
    cout << "O16 state coincidence Efficiency = " << coinc_eff << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*

    TCanvas *c4 = new TCanvas("c4");
    tsingles->Draw("ExcitationEnergy_C12>>h4(2000,-5,15)","ThetaLab>50 && ThetaLab<55");

    TF1 *f3 = new TF1("f3","gaus(0)+pol2(3)",4.1,4.8);
    f3->SetParameter(0,8000.);
    f3->SetParameter(1,4.44);
    f3->SetParameter(2,0.07);
    f3->SetParameter(3,500);
    f3->SetParameter(4,500);
    f3->SetParameter(5,0.1);
    h4->Fit("f3","R+");

    bin = 2000./20.;            // bins per MeV
    amp = f3->GetParameter(0);
    sigma = f3->GetParameter(2);
    singles = (bin * amp * sigma) / constant;  

    TCanvas *c5 = new TCanvas("c5");
    tcoinc->Draw("ExcitationEnergy_C12:GammaEnergy>>h5(10000,0,10000,2000,-5,15)","ThetaLab>50 && ThetaLab<55","colz");

    TCanvas *c6 = new TCanvas("c6");
    tcoinc->Draw("GammaEnergy>>h6(5000,0,5000","ThetaLab>50 && ThetaLab<55 && ExcitationEnergy_C12>4.25 && ExcitationEnergy_C12<4.61");
    h5->Rebin(4);

    TF1 *f4 = new TF1("f4","gaus(0)+pol1(3)",4300,4500);
    f4->SetParameter(0,150.);
    f4->SetParLimits(0,50.,400.);
    f4->SetParameter(1,2735);
    f4->SetParLimits(1,2700.,2750.);
    f4->SetParameter(2,15);
    f4->SetParLimits(2,10.,30.);
    f4->SetParameter(3,150);
    f4->SetParameter(4,-0.1);
    h6->Fit("f4","R+");

    bin = 1./4.;            // bins per keV
    amp = f4->GetParameter(0);
    sigma = f4->GetParameter(2);
    coincidences = h6->Integral(4360,4505);

    branching = 1.0;
    efficiency = 0.026;
    corrected_coincidences = coincidences / (branching * efficiency);

    coinc_eff = corrected_coincidences / singles;

    cout << "C12 state singles counts = " <<  singles << endl;
    cout << "C12 state coincidence counts = " <<  coincidences << endl;
    cout << "C12 state corrected coincidence counts = " <<  corrected_coincidences << endl;
    cout << "C12 state coincidence Efficiency = " << coinc_eff << endl;
    */
}