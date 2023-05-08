{
    // Open data files
    TFile *fdata = TFile::Open("./root_outputs/SumRuns_59-93.root", "READ");            // data file
    TFile *fsubtracted = TFile::Open("./plotting_scripts/subtraction.root", "READ");    // subtracted singles

    // Open output file
    TFile *outfile = new TFile("output.root","RECREATE");                       // output file

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

    // Creat raw excitation energy histogram to be used by subtraction program
    TH2D *hExcAng = new TH2D("hExcAng","Excitation Energy vs Angle",180,0,180,2000,-5,15);
    TTree *tsingles = (TTree*)fdata->Get("t1");
    tsingles->Draw("Exc:ThetaLab>>hExcAng",barrel + barrel_protons + angle_cut,"goff");
    hExcAng->Write();

    // Create 2D histogram of Exc vs Egam for events within the subtractable range in the barrel
    TH2D *hExc_v_Egam = new TH2D("hExc_v_Egam","Excitation vs Gamma Energy",10000,0,10000,2000,-5,15);
    TTree *tcoinc = (TTree*)fdata->Get("t2");
    tcoinc->Draw("Exc:Egam>>hExc_v_Egam",barrel + barrel_protons + angle_cut,"goff");
    hExc_v_Egam->Write();

    // Grab and write subtracted singles histogram for barrel events
    TH1D *hExc = (TH1D*)fsubtracted->Get("hExc_Subtracted_Sum");
    hExc->Write();

    //Check
    TCanvas *c1 = new TCanvas("c1");
    tsingles->Draw("Exc>>h1(2000,-5,15)",barrel + barrel_protons + angle_cut);
    hExc->SetLineColor(2);
    hExc->Draw("same");

}