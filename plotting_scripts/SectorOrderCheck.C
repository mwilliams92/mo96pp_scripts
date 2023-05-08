{

	gROOT->SetBatch(kTRUE);

	TChain* Chain = new TChain ("data");
/*  Chain->Add("../OutputFolder/Run0030.root");
	Chain->Add("../OutputFolder/Run0031.root");
	Chain->Add("../OutputFolder/Run0032.root");
	Chain->Add("../OutputFolder/Run0033.root");
	Chain->Add("../OutputFolder/Run0034.root");
	Chain->Add("../OutputFolder/Run0035.root");
	Chain->Add("../OutputFolder/Run0037.root");
	Chain->Add("../OutputFolder/Run0038.root");
	Chain->Add("../OutputFolder/Run0039.root");
	Chain->Add("../OutputFolder/Run0040.root");
	Chain->Add("../OutputFolder/Run0041.root");
	Chain->Add("../OutputFolder/Run0042.root");
	Chain->Add("../OutputFolder/Run0043.root");
	Chain->Add("../OutputFolder/Run0044.root");
	Chain->Add("../OutputFolder/Run0045.root");
	Chain->Add("../OutputFolder/Run0046.root");
	Chain->Add("../OutputFolder/Run0047.root");
	Chain->Add("../OutputFolder/Run0048.root");
	Chain->Add("../OutputFolder/Run0049.root");
*/
	//Chain->Add("../OutputFolder/cal228th_315deg_US_BR.root");
	//Chain->Add("../OutputFolder/cal228th_45deg_US_BL.root");

	Chain->Add("../OutputFolder/Run0080.root"); // CD2


	TFile *root_outfile = new TFile("OrderCheck.root", "RECREATE");

	TH2D *hOrder[12];

	char cut[4096];
	char draw[4096];
	char hname[4096];


	// Plot Position vs Sector for Each detector. Counts should increase with position. Position should increase with sector.
	for (int i=0; i<12; i++) {

			cout << "Detector " << i << endl; 

			sprintf(hname,"hOrder[%d]", i);
			sprintf(draw,"SX3StripPositionCal:SX3Sector>>hOrder[%d]",i);
			sprintf(cut,"SX3Upstream==0 && SX3Det==%d", i);

			hOrder[i] = new TH2D(hname,hname,8,-4,4,100,0,1);
			Chain->Draw(draw,cut,"goff");		
			
	}
	
	root_outfile->Write();
	root_outfile->Close();	


}
