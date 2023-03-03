// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

// Select two VBF jets by the setting maximum delta eta
void selectVBFs(int& i, int& j, std::vector<LorentzVector> jets) {
    double deltaEta=-999.999;
    for (unsigned int jetsi = 0; jetsi < jets.size(); jetsi++) {
        for (unsigned int jetsj = jetsi+1; jetsj < jets.size(); jetsj++) {
            double this_deltaE=fabs(jets[jetsi].Eta()-jets[jetsj].Eta());
            if (this_deltaE > deltaEta) {
                deltaEta=this_deltaE;
                i=jetsi;
                j=jetsj;
            }
        }
    }
}

// Function to select opposite-sign leptons (electrons & muons); the two leptons have to be the same type
// Select Z->l+l-
bool eq2OSLeptonsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow, double pt=30, bool sf=true) {
    int n_leptons = 0;
    std::vector<LorentzVector> elec{};
    std::vector<LorentzVector> muon{};
    std::vector<int> lepindx{};
    std::vector<int> pdgids{}; // Save to cutflow.globals
    int os=1;
    bool logic=false;
    // Select electrons
    for (unsigned int elec_i = 0; elec_i < nt.Electron_pt().size(); ++elec_i)
    {
        double this_pt=nt.Electron_pt().at(elec_i);
        LorentzVector this_elec=nt.Electron_p4().at(elec_i);
        int this_elec_id=nt.Electron_jetIdx()[elec_i];
        int this_pdgid=nt.Electron_pdgId()[elec_i];
        // Use IDtight so that the requirement is for analysis
        if (this_pt > 30 && ttH_UL::electronID(elec_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            elec.push_back(this_elec);
            lepindx.push_back(this_elec_id);
            os=os*this_pdgid;
            pdgids.push_back(this_pdgid);
        }
    }
    // Select Muons
    for (unsigned int muon_i = 0; muon_i < nt.Muon_pt().size(); ++muon_i) {
        double this_pt=nt.Muon_pt().at(muon_i);
        LorentzVector this_muon=nt.Muon_p4().at(muon_i);
        int this_muon_id=nt.Muon_jetIdx()[muon_i];
        int this_pdgid=nt.Muon_pdgId()[muon_i];
        if (this_pt > 30 && ttH_UL::muonID(muon_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            muon.push_back(this_muon);
            lepindx.push_back(this_muon_id);
            os=os*this_pdgid;
            pdgids.push_back(this_pdgid);
        }
    }
    // Require events with exactly two opposite sign leptons of same flavor (muon/antimuon, electron/positron pairs)
    /* count the number of electrons (n_electrons veto == 2) */ 
    if (sf) {logic=(n_leptons==2 && (elec.size()==2 || muon.size()==2) && os<0);}
    else {logic=(n_leptons==2 && os<0);}
    if (logic) {
        LorentzVector ld, sd;
        int ldid=-999;
        int sdid=-999;
        // Select leading and subleading leptons with their pdgids
        if (elec.size()==2) {
            ld=(elec.at(0)).Pt()>(elec.at(1)).Pt() ? elec.at(0):elec.at(1);
            sd=(elec.at(1)).Pt()>(elec.at(0)).Pt() ? elec.at(0):elec.at(1);
            ldid=(elec.at(0)).Pt()>(elec.at(1)).Pt() ? pdgids.at(0):pdgids.at(1);
            sdid=(elec.at(1)).Pt()>(elec.at(0)).Pt() ? pdgids.at(0):pdgids.at(1);
        }
        else if (muon.size()==2) {
            ld=(muon.at(0)).Pt()>(muon.at(1)).Pt() ? muon.at(0):muon.at(1);
            sd=(muon.at(1)).Pt()>(muon.at(0)).Pt() ? muon.at(0):muon.at(1);
            ldid=(muon.at(0)).Pt()>(muon.at(1)).Pt() ? pdgids.at(0):pdgids.at(1);
            sdid=(muon.at(1)).Pt()>(muon.at(0)).Pt() ? pdgids.at(0):pdgids.at(1);
        }
        cutflow.globals.setVal<LorentzVector>("sd_lep_p4", sd);
        cutflow.globals.setVal<LorentzVector>("ld_lep_p4", ld);
        cutflow.globals.setVal<std::vector<int>>("jetidx", lepindx); // JetIdx
        cutflow.globals.setVal<int>("ld_lep_pdgid",ldid);
        cutflow.globals.setVal<int>("sd_lep_pdgid",sdid);
    }
    // Require only two lepton events to pass
    return logic;
}


// Preselection functions
// Function to select VBF system mass >= 500
bool mjjGt500(Nano& nt, Arbol& arbol, Cutflow& cutflow, double mjj=500) {
    bool logic = (cutflow.globals.getVal<double>("mjj")>=mjj);
    return logic;
}

// Function to select delta eta between VBF jets >= 3, then fill the selections
bool deltaEtaGt3(Nano& nt, Arbol& arbol, Cutflow& cutflow, double mindeta=3) {
    bool logic = (cutflow.globals.getVal<double>("detajj")>=mindeta);
    return logic;
}