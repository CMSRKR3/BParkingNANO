#include <iostream>
#include <string>
#include <vector>
#include <TFile.h>
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TMath.h>
#include <TTree.h>
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector2.h>               // 2D vector class
#include <TMath.h>                  // ROOT math library
#include <TLeaf.h>
#include <list>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <typeinfo>
#include <TGraph.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TNtuple.h> 
#include <TH1.h>
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RVec.hxx"
#include <TFile.h>
#include <TH1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <iostream>

using namespace std;

int bparkingplotter(){
  TTree* eventTree=0;
  TFile *inFile = TFile::Open("/eos/user/j/jodedra/BPARKINGNANOSTUFF/CMSSW_12_4_0_pre3/src/PhysicsTools/BParkingNano/test/BParkNANOCharmoniumDataDataset10000evts_data_124X.root", "READ");
  TFile *f = new TFile("outputDATA200000evts.root", "RECREATE");
  TTree *T = new TTree("myTree","a tree");
  Float_t bmass = 0;
  T->Branch("bmass", &bmass, "bmass/Float_t");
  
  assert(inFile);
  TTreeReader fReader;
  TTreeReaderArray<Float_t> mass = {fReader, "BToKMuMu_fit_mass"};
  TTreeReaderValue<UInt_t> run = {fReader, "run"};

  eventTree = (TTree*)inFile->Get("Events");
  assert(eventTree);
  fReader.SetTree(eventTree);
  std::cout<<eventTree->GetEntries()<<endl;
  int counter = 0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++){
    fReader.SetLocalEntry(ientry);
    for(UInt_t i=0; i<mass.GetSize();i++){
	      std::cout<<mass[i]<<endl;
	      bmass = (Float_t)mass[i];
        counter++;
	      T->Fill(); 	
    }	
  }
  std::cout<<counter<<endl;
  f->Write();
  f->Close();
  return 0;
  
}
