void evtrunlumi() {

	TChain* inChain = new TChain("Events");
	inChain->Add("Final_EGM_id/event_ids__2018__egamma__CR_Gamma_fail.root");
	inChain->Add("Final_EGM_id/event_ids__2018__egamma__CR_Gamma_pass.root");

    cout << "TotalNEvents " << inChain->GetEntries() << endl;
	unsigned int NumEvents = inChain->GetEntries();


    Long64_t evt, run, lumi;
    inChain->SetBranchAddress("Evt_ID", &evt);
    inChain->SetBranchAddress("Evt_Run", &run);
    inChain->SetBranchAddress("Evt_Lumi", &lumi);

    std::ofstream outFile("egm_GCR_Moritz.txt");

    outFile << "event number\t" << "run number\t" << "lumi number" << std::endl;

    Long64_t nEntries = inChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        inChain->GetEntry(i);
        outFile << evt << "\t" << run << "\t" << lumi << std::endl;
    }

    outFile.close();

}

