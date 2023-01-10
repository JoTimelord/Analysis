//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan  8 20:39:44 2023 by ROOT version 6.22/09
// from TTree Events/
// found on file: ../tempdata/VBSOSWWH_Inclusive_4f_LO_10k.root
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Long64_t> NUP = {fReader, "NUP"};
   TTreeReaderValue<Long64_t> IDPRUP = {fReader, "IDPRUP"};
   TTreeReaderValue<Double_t> XWGUP = {fReader, "XWGUP"};
   TTreeReaderValue<Double_t> SCALUP = {fReader, "SCALUP"};
   TTreeReaderValue<Double_t> AQEDUP = {fReader, "AQEDUP"};
   TTreeReaderValue<Double_t> AQCDUP = {fReader, "AQCDUP"};
   TTreeReaderValue<Int_t> nIDUP = {fReader, "nIDUP"};
   TTreeReaderArray<Long64_t> IDUP = {fReader, "IDUP"};
   TTreeReaderValue<Int_t> nISTUP = {fReader, "nISTUP"};
   TTreeReaderArray<Long64_t> ISTUP = {fReader, "ISTUP"};
   TTreeReaderValue<Int_t> nMOTHUP1 = {fReader, "nMOTHUP1"};
   TTreeReaderArray<Long64_t> MOTHUP1 = {fReader, "MOTHUP1"};
   TTreeReaderValue<Int_t> nMOTHUP2 = {fReader, "nMOTHUP2"};
   TTreeReaderArray<Long64_t> MOTHUP2 = {fReader, "MOTHUP2"};
   TTreeReaderValue<Int_t> nICOLUP1 = {fReader, "nICOLUP1"};
   TTreeReaderArray<Long64_t> ICOLUP1 = {fReader, "ICOLUP1"};
   TTreeReaderValue<Int_t> nICOLUP2 = {fReader, "nICOLUP2"};
   TTreeReaderArray<Long64_t> ICOLUP2 = {fReader, "ICOLUP2"};
   TTreeReaderValue<Int_t> nP_X = {fReader, "nP_X"};
   TTreeReaderArray<Double_t> P_X = {fReader, "P_X"};
   TTreeReaderValue<Int_t> nP_Y = {fReader, "nP_Y"};
   TTreeReaderArray<Double_t> P_Y = {fReader, "P_Y"};
   TTreeReaderValue<Int_t> nP_Z = {fReader, "nP_Z"};
   TTreeReaderArray<Double_t> P_Z = {fReader, "P_Z"};
   TTreeReaderValue<Int_t> nE = {fReader, "nE"};
   TTreeReaderArray<Double_t> E = {fReader, "E"};
   TTreeReaderValue<Int_t> nM = {fReader, "nM"};
   TTreeReaderArray<Double_t> M = {fReader, "M"};
   TTreeReaderValue<Int_t> nVTIMUP = {fReader, "nVTIMUP"};
   TTreeReaderArray<Double_t> VTIMUP = {fReader, "VTIMUP"};
   TTreeReaderValue<Int_t> nSPINUP = {fReader, "nSPINUP"};
   TTreeReaderArray<Double_t> SPINUP = {fReader, "SPINUP"};
   TTreeReaderArray<Double_t> rwgt = {fReader, "rwgt"};


   Selector(TTree * /*tree*/ =0) { }
   virtual ~Selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Selector,0);

};

#endif

#ifdef Selector_cxx
void Selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Selector_cxx
