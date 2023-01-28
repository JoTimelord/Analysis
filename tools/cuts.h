// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

// Function to select VBFjets
void selectVBFs(int& i, int& j, std::vector<LorentzVector> jets) {
    double largestEta1 = -999;
    double largestEta2 = 999;
    for (unsigned int jetsi = 0; jetsi < jets.size(); jetsi++) {
        double jeteta=jets[jetsi].Eta();
        if (jeteta > 0 && jeteta >= largestEta1) {
            i = jetsi;
            largestEta1 = jeteta;
        }
        else if (jeteta < largestEta2 && jeteta < 0) {
            j = jetsi;
            largestEta2 = jeteta;
        }
    }
}

// Function to select jets with largest pt's
// void selectLDs(int& i, int &)
// Function to select leptons (electrons & muons); the two leptons can be mix and match
bool eq2ElectronsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_leptons = 0;
    float maxpt1=-999.999;
    float maxpt2=-9999.99;
    int ldid=-999;
    int sdid=-999;
    std::vector<LorentzVector> leps{};
    std::vector<int> lepindx{};
    // Select electrons
    for (unsigned int elec_i = 0; elec_i < nt.Electron_pt().size(); ++elec_i)
    {
        float pt=nt.Electron_pt().at(elec_i);
        // Use IDtight so that the requirement is for analysis
        if (pt > 30 && ttH_UL::electronID(elec_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            leps.push_back(nt.Electron_p4().at(elec_i));
            lepindx.push_back(nt.Electron_jetIdx()[elec_i]);
        }
    }
    // Select Muons
    for (unsigned int muon_i = 0; muon_i < nt.Muon_pt().size(); ++muon_i) {
        float pt=nt.Muon_pt().at(muon_i);
        if (pt > 30 && ttH_UL::muonID(muon_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            leps.push_back(nt.Muon_p4().at(muon_i));
            lepindx.push_back(nt.Muon_jetIdx()[muon_i]);
        }
    }
    // Select leading and subleading leptons
    for (long unsigned int i=0; i<leps.size(); i++) {
        float pt=(leps.at(i)).Pt();
        if (pt > maxpt1) {
            maxpt2=maxpt1;
            maxpt1=pt;
            sdid=ldid;
            ldid=i;
        }
        else if (pt > maxpt2 && pt <= maxpt1) {
            maxpt2=pt;
            sdid=i;
        }
    }
    // Require only two lepton events to pass
    if (n_leptons == 2) {
        cutflow.globals.setVal<LorentzVector>("sd_lep_p4", leps.at(sdid));
        cutflow.globals.setVal<LorentzVector>("ld_lep_p4", leps.at(ldid));
        cutflow.globals.setVal<std::vector<int>>("jetidx", lepindx);
    }
    return (n_leptons = 2);
}

// Function to select events with at least one fatjet (boosted Higgs); check overlap against leptons; check hbb score >= 0.9 or wjj >= 0.95
bool geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_fatjet=0;
    std::vector<LorentzVector> fatjets{};
    float highesthbb=-999.999;
    float highestwjj=-999.999;
    int hbbindx=0;
    int wjjindx=0;
    for (unsigned int ifatjet = 0; ifatjet < nt.FatJet_pt().size(); ++ifatjet) {
        // Check against overlap with leptons by checking difference in pseudorapidity (Eta) and Phi
        bool isOverlap = false;
        if (ROOT::Math::VectorUtil::DeltaR(cutflow.globals.getVal<LorentzVector>("sd_lep_p4"), nt.FatJet_p4()[ifatjet]) < 0.8) {
            isOverlap = true;
            break;
        }
        if (ROOT::Math::VectorUtil::DeltaR(cutflow.globals.getVal<LorentzVector>("ld_lep_p4"), nt.FatJet_p4()[ifatjet]) < 0.8) {
            isOverlap = true;
            break;
        }
        if (isOverlap) {continue; }
        // Check ak8 jets with softdropmass larger than 40
        if (nt.FatJet_msoftdrop()[ifatjet] >=40 ) {
            n_fatjet++;
            // Select boosted higgs
            if (nt.FatJet_particleNet_HbbvsQCD()[ifatjet] > highesthbb) { 
                hbbindx = ifatjet; 
                highesthbb = nt.FatJet_particleNet_HbbvsQCD()[ifatjet];
            }
            if (nt.FatJet_particleNet_WvsQCD()[ifatjet] > highestwjj) {
                wjjindx = ifatjet;
                highestwjj = nt.FatJet_particleNet_WvsQCD()[ifatjet];
            }
        }
    }
    
    // Select Logic: Pass events with at least one fatjets
    bool logic = (n_fatjet >= 1 && (highestwjj >= 0.95 || highesthbb >= 0.9));
    if (logic) {
        LorentzVector p4;
        p4.SetCoordinates(nt.FatJet_pt()[hbbindx], nt.FatJet_eta()[hbbindx], nt.FatJet_phi()[hbbindx], nt.FatJet_msoftdrop()[hbbindx]);
        cutflow.globals.setVal<LorentzVector>("fatjet_p4", p4);
    }
    return (logic);
}

// Function to select events with at least two ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against fatjets
// Select VBF jets
/* Is there a tight jet requirement? */
bool geq2JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    std::vector<LorentzVector> jets{};
    int n_jets=0;
    int VBF1=0;
    int VBF2=1;
    int ld=0;
    int sd=1;
    for (unsigned int ijet = 0; ijet < nt.Jet_pt().size(); ++ijet) {
        // Overlap check against leptons
        bool isOverlap=false;
        if (cutflow.globals.getVal<std::vector<int>>("jetidx").at(0)== -1 || cutflow.globals.getVal<std::vector<int>>("jetidx").at(0)== -1) {
            if (ROOT::Math::VectorUtil::DeltaR(cutflow.globals.getVal<LorentzVector>("ld_lep_p4"), nt.Jet_p4().at(ijet)) < 0.4 || ROOT::Math::VectorUtil::DeltaR(cutflow.globals.getVal<LorentzVector>("sd_lep_p4"), nt.Jet_p4().at(ijet)) < 0.4)
            {
                isOverlap = true;
                break;
            }
        }
        else if (cutflow.globals.getVal<std::vector<int>>("jetidx").at(0)== (int) ijet || cutflow.globals.getVal<std::vector<int>>("jetidx").at(1)== (int) ijet)
        {
            isOverlap = true;
            break;
        }
        
        // Overlap check against fatjet
        if (ROOT::Math::VectorUtil::DeltaR(cutflow.globals.getVal<LorentzVector>("fatjet_p4"), nt.Jet_p4().at(ijet)) <= 0.8) {
            isOverlap = true;
            break;
        }
        if (isOverlap) {continue; }
        if (nt.Jet_pt().at(ijet) < 30.0) {continue; }
        n_jets++;
        jets.push_back(nt.Jet_p4().at(ijet));
    }
    if (n_jets >= 2) { 
        selectVBFs(VBF1, VBF2, jets);
        cutflow.globals.setVal<LorentzVector>("ld_vbf_p4", nt.Jet_p4().at(VBF1));
        cutflow.globals.setVal<LorentzVector>("sd_vbf_p4", nt.Jet_p4().at(VBF2));
    }
    return (n_jets >= 2);
}


class AnalysisCut : public Cut
{
public:
    Arbol& arbol;
    Nano& nt;
    Utilities::Variables& globals;

    AnalysisCut(std::string name, Nano& _nt, Arbol& _arbol, Cutflow& _cutflow)
    : Cut(name), arbol(_arbol), nt(_nt), globals(_cutflow.globals)
    {
        // Do nothing
    }

    virtual ~AnalysisCut() {}
};

class SelectLHEVariables : public AnalysisCut
{
public:
    SelectLHEVariables(std::string name, Nano& nt, Arbol& arbol, Cutflow& cutflow)
    : AnalysisCut(name, nt, arbol, cutflow)
    {
        // Do nothing
    }
    
    // LHE levels are truth-level and doesn't require selections so this function always returns true
    // It simply sets LHELevel variables in arbol tree
    bool evaluate() {
        unsigned int particleNo=nt.LHEPart_status().size();
        std::vector<LorentzVector> jets;
        LorentzVector W,H,Z;
        int VBFindx1=0;
        int VBFindx2=1;
        
        // Selection based on pdgid
        for (unsigned int i=0; i<particleNo; i++) {
            // Check if it's outgoing particles
            if (nt.LHEPart_status()[i]==1) {
                Int_t pdgid=TMath::Abs(nt.LHEPart_pdgId()[i]);
                LorentzVector p4;
                p4.SetCoordinates(nt.LHEPart_pt()[i],nt.LHEPart_eta()[i],nt.LHEPart_phi()[i],nt.LHEPart_mass()[i]);
                switch (pdgid) {
                    case 24:
                        W=p4;
                        break;
                    case 25:
                        H=p4;
                        break;
                    case 23:
                        Z=p4;
                        break;
                    default:
                        jets.push_back(p4);
                        break;
                } 
            }
        }
        
        // Sanity check to make sure the LHE level has no corrupted info
        if (jets.size()<2) {return false;}

        // Select the VBF jets with the largest etas in both directions
        selectVBFs(VBFindx1, VBFindx2, jets);
        
        // Set variables to arbol tree, 'mvvh', 'lt', 'st', 'mjj', 'detajj'
        arbol.setLeaf<double>("mvvh", (W+H+Z).M());
        arbol.setLeaf<double>("lt", Z.Pt());
        arbol.setLeaf<double>("st", (W+H+Z).Pt());
        arbol.setLeaf<double>("mjj", (jets[VBFindx1]+jets[VBFindx2]).M());
        arbol.setLeaf<double>("detajj", (TMath::Abs(jets[VBFindx1].Eta()-jets[VBFindx2].Eta())));

        // Cut logic 
        return true;
    }

    double weight()
    {
        /* Event weight (applied only if Cut::evaluate() returns true) */
        return 1.;
    };

    virtual ~SelectLHEVariables() {}
};