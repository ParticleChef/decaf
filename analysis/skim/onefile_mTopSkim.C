void mTopSkim() {

    std::ifstream fileList("EGM18Dlist.txt");
    std::string fileName;

	TChain* inChain = new TChain("Events");

    while(std::getline(fileList, fileName)){
        inChain->Add(fileName.c_str());
    }
	//inChain->Add("/data/data/2018D/EGammaD/DATA_2018_NanoAODv9_1000.root");
	fileList.close();

    cout << "TotalNEvents " << inChain->GetEntries() << endl;
	unsigned int NumEvents = inChain->GetEntries();
	int per99 = NumEvents / 100;
	int per100 = 0;


    Bool_t goodVertices, globalSuperTightHalo2016Filter, HBHENoiseFilter, HBHENoiseIsoFilter, EcalDeadCellTriggerPrimitiveFilter, BadPFMuonFilter, BadPFMuonDzFilter, ecalBadCalibFilter;
    inChain->SetBranchAddress("Flag_goodVertices", &goodVertices);
    inChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &globalSuperTightHalo2016Filter);
    inChain->SetBranchAddress("Flag_HBHENoiseFilter", &HBHENoiseFilter);
    inChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &HBHENoiseIsoFilter);
    inChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter);
    inChain->SetBranchAddress("Flag_BadPFMuonFilter", &BadPFMuonFilter);
    inChain->SetBranchAddress("Flag_BadPFMuonDzFilter", &BadPFMuonDzFilter);
    inChain->SetBranchAddress("Flag_ecalBadCalibFilter", &ecalBadCalibFilter);

    Bool_t Ele32_WPTight_Gsf, Photon200, PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, PFMETNoMu120_PFMHTNoMu120_IDTight;
    inChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &Ele32_WPTight_Gsf);
    inChain->SetBranchAddress("HLT_Photon200", &Photon200);
    inChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60);
    inChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &PFMETNoMu120_PFMHTNoMu120_IDTight);

    const unsigned int maxEvents = 500000;
    unsigned int eventCount = 0;
    unsigned int fileCount = 0;

	//TFile* outFile = new TFile("out.root", "recreate");
	//TTree* outTree = inChain->CloneTree(0);
	TFile* outFile = nullptr;
	TTree* outTree = nullptr;
	int nPass = 0;
	for(unsigned int eventLoop=0; eventLoop < NumEvents; eventLoop++) {
		inChain->GetEntry(eventLoop);
		//if(eventLoop % per99 == 10) cout << "processing " << per100++ << " %" << endl;
		if((goodVertices & globalSuperTightHalo2016Filter & HBHENoiseFilter & HBHENoiseIsoFilter & EcalDeadCellTriggerPrimitiveFilter & BadPFMuonFilter & BadPFMuonDzFilter * ecalBadCalibFilter) &
            (Ele32_WPTight_Gsf | Photon200 | PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 | PFMETNoMu120_PFMHTNoMu120_IDTight)) {
            if(eventCount % maxEvents == 0){
                if(outFile){
                    outFile->Write();
                    outFile->Close();
                } // if outFile is not none
                std::ostringstream oss;
                oss << "/data/mTopnTuples/monotop2018/EGammaD/DATA_2018_NanoAODv9_skim_" << fileCount++ << ".root";
                outFile = new TFile(oss.str().c_str(), "recreate");
                outTree = inChain->CloneTree(0);
                std::cout<<"fileCount: "<<fileCount<<std::endl;
            } // if eventCount is n*maxEventsPerFile

			outTree->Fill();
            eventCount++;
			nPass++;
		} // skim variables condition
	} // eventLoop	
    if(outFile){
	    cout << "skim0 nPass " << nPass << " outTree " << outTree->GetEntries() << endl;
	    outFile->Write();
        outFile->Close();
    }
}

/****
void skimNano1() {

	TChain* inChain = new TChain("Events");
	inChain->Add("orig/dataNanoAOD.root");
	cout << "TotalNEvents " << inChain->GetEntries() << endl;

	//inChain->Draw(">>myList","MET_pt > 150.0","entryList");
	inChain->Draw(">>myList","nElectron > 1","entryList");
	TEntryList *myList = (TEntryList*)gDirectory->Get("myList");
	myList->Print("all");

	TFile* outFile = new TFile("out1.root", "recreate");
	inChain->SetEntryList(myList);
	TTree* outTree = inChain->CopyTree("");
	cout << " skim1 outTree " << outTree->GetEntries() << endl;
	outFile->Write();	

}
****/

