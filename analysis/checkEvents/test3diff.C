#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm> // For std::find
#include <chrono> // For measure timing

// Function to load event numbers from a text file into a vector
//std::vector<std::tuple<ULong64_t,UInt_t, UInt_t>> loadEventList(const char* filename) {
//    std::vector<std::tuple<ULong64_t,UInt_t, UInt_t>> eventList;
std::vector<std::tuple<ULong64_t,UInt_t>> loadEventList(const char* filename) {
    std::vector<std::tuple<ULong64_t, UInt_t>> eventList;
    std::ifstream infile(filename);
    UInt_t run;
    ULong64_t event;
    UInt_t cnt = 0;
    //UInt_t lumi;

    //while (infile >>event>>run>>lumi) {
    //    eventList.emplace_back(event, run, lumi);
    //}
    while (infile >>event>>run) {
        if(cnt == 0) {
            std::cout<<"event: "<<event<<"   run: "<<run<<endl;
            cnt++;
        }
        eventList.emplace_back(event, run);
    }
    return eventList;
}

// Function to save selected events based on event numbers list
void saveSelectedEventsFromList(const char* inputFileList, const char* outputFileName, const char* treeName, const char* eventListFile) {
    // Load the event numbers from the text file
    std::vector<std::tuple<ULong64_t, UInt_t>> eventList = loadEventList(eventListFile);

    // Create a TChain to combine multiple input files
    TChain chain(treeName);
    std::ifstream infile(inputFileList);
    TString fileName;

    while (infile >> fileName) {
        chain.Add(fileName);
    }

    if (chain.GetEntries() == 0) {
        std::cerr << "Error: No entries found in the input files." << std::endl;
        return;
    }

    // Create a new file to save the selected events
    TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || !outputFile->IsOpen()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        return;
    }

    // Create a new TTree to store selected events
    TTree *outputTree = chain.CloneTree(0); // Clone the tree structure, but not the data

    // Variables to hold event data, adjust according to your TTree structure
    ULong64_t event;//, run, lumi;
    UInt_t run;
    UInt_t lumi;

    chain.SetBranchAddress("run", &run);
    chain.SetBranchAddress("event", &event);
    //chain.SetBranchAddress("luminosityBlock", &lumi);

    // Loop over all events and apply the selection criteria
    Long64_t nEntries = chain.GetEntries();
    Long64_t matchedEntries = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);

        // Check if the event number is in the list
        //if (std::find(eventList.begin(), eventList.end(), std::make_tuple(event,run,lumi)) != eventList.end()) {
        if (std::find(eventList.begin(), eventList.end(), std::make_tuple(event,run)) != eventList.end()) {
            outputTree->Fill(); // Save the event if it matches the list
            matchedEntries++;
        }
    }

    std::cout << "Number of Matched entries " << matchedEntries << std::endl;

    // Write the output tree to the file
    outputTree->Write();
    outputFile->Close();

    std::cout << "Selected events have been saved to " << outputFileName << std::endl;
}

// Example usage function
void test3diff(const char* input, const char* output, const char* tree, const char* textfile) {
    auto start_time = std::chrono::high_resolution_clock::now();

   // const char* inputFileList = "inputrootlist.txt"; // List of input ROOT files
   // const char* outputFileName = "selectedEvents.root"; // Output ROOT file name
   // const char* treeName = "Events"; // Name of the TTree in input ROOT files
   // const char* eventListFile = "EventNo_mTop.txt"; // Event number list file

    saveSelectedEventsFromList(input, output, tree, textfile);

    auto end_time = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
}

//// Example usage function
//void test2() {
//    auto start_time = std::chrono::high_resolution_clock::now();
//
//    const char* inputFileList = "inputrootlist.txt"; // List of input ROOT files
//    const char* outputFileName = "selectedEvents.root"; // Output ROOT file name
//    const char* treeName = "Events"; // Name of the TTree in input ROOT files
//    const char* eventListFile = "EventNo_mTop.txt"; // Event number list file
//
//    saveSelectedEventsFromList(inputFileList, outputFileName, treeName, eventListFile);
//
//    auto end_time = std::chrono::high_resolution_clock::now();
//    
//    std::chrono::duration<double> elapsed = end_time - start_time;
//    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;
//}
