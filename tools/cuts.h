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

// Function to select two jets with minimum dR for W reconstruction
void selectWjets(int& i, int& j, int VBF1, int VBF2, std::vector<LorentzVector> jets) {
    double mindR=999;
    for (int ijet=0; ijet<jets.size(); ijet++) {
        if (ijet!=VBF1 && ijet!=VBF2) {
            for (int jjet=ijet+1; jjet<jets.size(); jjet++) {
                if (jjet!=VBF1 && jjet!=VBF2) {
                    if (ROOT::Math::VectorUtil::DeltaR(jets.at(ijet),jets.at(jjet)) < mindR) {
                        mindR=ROOT::Math::VectorUtil::DeltaR(jets.at(ijet),jets.at(jjet));
                        i=ijet;
                        j=jjet;
                    }
                }
            }
        }
    }
}

// Function to fill reconstruction level arbol tree
void fillTree(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    LorentzVector Z=cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    LorentzVector H=cutflow.globals.getVal<LorentzVector>("fatjet_p4");
    LorentzVector W=cutflow.globals.getVal<LorentzVector>("ld_jet_p4")+cutflow.globals.getVal<LorentzVector>("sd_jet_p4");
    arbol.setLeaf<LorentzVector>("vbfjet1_LV",cutflow.globals.getVal<LorentzVector>("ld_vbf_p4"));
    arbol.setLeaf<LorentzVector>("vbfjet2_LV",cutflow.globals.getVal<LorentzVector>("sd_vbf_p4"));
    arbol.setLeaf<LorentzVector>("wjet1_LV",cutflow.globals.getVal<LorentzVector>("ld_jet_p4"));
    arbol.setLeaf<LorentzVector>("wjet2_LV",cutflow.globals.getVal<LorentzVector>("sd_jet_p4"));
    arbol.setLeaf<LorentzVector>("lep1_LV",cutflow.globals.getVal<LorentzVector>("ld_lep_p4"));
    arbol.setLeaf<LorentzVector>("lep2_LV",cutflow.globals.getVal<LorentzVector>("sd_lep_p4"));
    arbol.setLeaf<int>("lep1_ID", );
    arbol_.newBranch<int>("lep2_ID", -999);
    arbol_.newBranch<double>("mll", -999);
    arbol_.newBranch<double>("ptll", -999.0);
    arbol_.newBranch<double>("met", -999.0);
    arbol_.newBranch<int>("n_ak4", -999);
    arbol_.newBranch<int>("n_ak8", -999);
    arbol_.newBranch<double>("hbb_score", -999);
    arbol_.newBranch<LorentzVector>("hbb_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("xbb_score", -999);
    arbol_.newBranch<double>("hbb_pnetmass", -999);
    arbol_.newBranch<double>("hbb_pnetmass", -999);
    arbol_.newBranch<bool>("b_veto", false);
    arbol.setLeaf<double>("mvvh", (W+H+Z).M());
    arbol.setLeaf<double>("lt", Z.Pt()+cutflow.globals.getVal<float>("met_pt"));
    arbol.setLeaf<double>("st", Z.Pt()+H.Pt()+W.Pt());
    arbol.setLeaf<double>("mjj", (cutflow.globals.getVal<LorentzVector>("ld_vbf_p4")+cutflow.globals.getVal<LorentzVector>("sd_vbf_p4")).M());
    arbol.setLeaf<double>("detajj", (TMath::Abs((cutflow.globals.getVal<LorentzVector>("ld_vbf_p4")).Eta()-(cutflow.globals.getVal<LorentzVector>("sd_vbf_p4")).Eta())));
}

// Function to select opposite-sign leptons (electrons & muons); the two leptons have to be the same type
// Select Z->l+l-
bool eq2ElectronsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_leptons = 0;
    float maxpt1=-999.999;
    float maxpt2=-9999.99;
    std::vector<LorentzVector> elec{};
    std::vector<LorentzVector> muon{};
    std::vector<int> lepindx{};
    int os=1;
    // Select electrons
    for (unsigned int elec_i = 0; elec_i < nt.Electron_pt().size(); ++elec_i)
    {
        float pt=nt.Electron_pt().at(elec_i);
        // Use IDtight so that the requirement is for analysis
        if (pt > 30 && ttH_UL::electronID(elec_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            elec.push_back(nt.Electron_p4().at(elec_i));
            lepindx.push_back(nt.Electron_jetIdx()[elec_i]);
            os=os*nt.Electron_pdgId()[elec_i];
        }
    }
    // Select Muons
    for (unsigned int muon_i = 0; muon_i < nt.Muon_pt().size(); ++muon_i) {
        float pt=nt.Muon_pt().at(muon_i);
        if (pt > 30 && ttH_UL::muonID(muon_i, ttH::IDtight, nt.year())) { 
            n_leptons++; 
            muon.push_back(nt.Muon_p4().at(muon_i));
            lepindx.push_back(nt.Muon_jetIdx()[muon_i]);
            os=os*nt.Muon_pdgId()[muon_i];
        }
    }
    // Require events with exactly two opposite sign leptons of same flavor (muon/antimuon, electron/positron pairs)
    /* count the number of electrons (n_electrons veto == 2) */ 
    bool logic=(n_leptons==2 && (elec.size()==2 || muon.size()==2) && os<0);
    if (not(logic)) {return false;}
    if (logic) {
        LorentzVector ld, sd;
        // Select leading and subleading leptons
        if (elec.size()==2) {
            ld=(elec.at(0)).Pt()>(elec.at(1)).Pt() ? elec.at(0):elec.at(1);
            sd=(elec.at(1)).Pt()>(elec.at(0)).Pt() ? elec.at(0):elec.at(1);
        }
        else if (muon.size()==2) {
            ld=(muon.at(0)).Pt()>(muon.at(1)).Pt() ? muon.at(0):muon.at(1);
            sd=(muon.at(1)).Pt()>(muon.at(0)).Pt() ? muon.at(0):muon.at(1);
        }
        cutflow.globals.setVal<LorentzVector>("sd_lep_p4", sd);
        cutflow.globals.setVal<LorentzVector>("ld_lep_p4", ld);
        cutflow.globals.setVal<std::vector<int>>("jetidx", lepindx);
    }
    
    // Require only two lepton events to pass
    return logic;
}

// Function to select events with at least one fatjet (boosted Higgs); check overlap against leptons
// Requirement on fatjet: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
bool geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_fatjet=0;
    std::vector<LorentzVector> fatjets{};
    float highesthbb=-999.999; // mass-correlated score
    float highestxbb=-999.999; // mass-decorrelated score
    int hbbindx=-999;
    int xbbindx=-999;
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
        // Check ak8 jets with softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
        if (nt.FatJet_msoftdrop()[ifatjet] >=40 && nt.FatJet_pt()[ifatjet]>250 && TMath::Abs(nt.FatJet_eta())<2.5 && nt.FatJet_jetId>0) {
            n_fatjet++;
            // Select boosted higgs by hbb
            if (nt.FatJet_particleNet_HbbvsQCD()[ifatjet] > highesthbb) { 
                hbbindx = ifatjet; 
                highesthbb = nt.FatJet_particleNet_HbbvsQCD()[ifatjet];
            }
            // Select boosted Higgs by xbb
            if (nt.FatJet_particleNetMD_Xbb()[ifatjet] > highestxbb) {
                xbbindx = ifatjet;
                highestxbb = nt.FatJet_particleNetMD_Xbb()[ifatjet];
            }
        }
    }
    // Select Logic: Pass events with at least one fatjets
    bool logic = (n_fatjet >= 1 && highesthbb >= 0.9);
    if (logic) {
        LorentzVector p4;
        p4.SetCoordinates(nt.FatJet_pt()[hbbindx], nt.FatJet_eta()[hbbindx], nt.FatJet_phi()[hbbindx], nt.FatJet_msoftdrop()[hbbindx]);
        cutflow.globals.setVal<LorentzVector>("fatjet_p4", p4);
        cutflow.globals.setVal<double>("hbb",highesthbb);
        cutflow.globals.setVal<double>("xbb",highestxbb);
        cutflow.globals.setVal<bool>("samejet", (hbbindx==xbbindx));
    }
    return logic;
}


// Function to make checks on fatjet score
bool fatJetHbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore=0.5) {
    bool logic = (cutflow.globals.getVal<double>("hbb") >= setscore);
    return logic;
}

// Function to select events with at least four ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against fatjets
// Select VBF jets
bool geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    std::vector<LorentzVector> jets{};
    int n_jets=0;
    int VBF1=0;
    int VBF2=1;
    int ld=0;
    int sd=1;
    int w1=2;
    int w2=3;
    bool pass=false;
    for (unsigned int ijet = 0; ijet < nt.Jet_pt().size(); ++ijet) {
        // Tight jet id requirement
        if (nt.year() == 2016)
        {
            if (nt.Jet_jetId()[ijet] < 1) // For 2016 apparently it's >= 1
                continue;
        }
        else
        {
            if (nt.Jet_jetId()[ijet] < 2) // "Tight" ID requirement while for others >= 2
                continue;
        }
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
    if (n_jets >= 4) { 
        selectVBFs(VBF1, VBF2, jets);
        selectWjets(w1, w2, VBF1, VBF2, jets);
        /* put into another lambda cut with mjj */
        pass=(TMath::Abs((nt.Jet_p4().at(VBF1)).Eta()-(nt.Jet_p4().at(VBF2)).Eta()) >=3);
    }
    if (pass) {
        cutflow.globals.setVal<LorentzVector>("ld_vbf_p4", nt.Jet_p4().at(VBF1));
        cutflow.globals.setVal<LorentzVector>("sd_vbf_p4", nt.Jet_p4().at(VBF2));
        cutflow.globals.setVal<LorentzVector>("ld_jet_p4", nt.Jet_p4().at(w1));
        cutflow.globals.setVal<LorentzVector>("sd_jet_p4", nt.Jet_p4().at(w2));
        //cutflow.globals.setVal<bool>("medium_b_veto", (nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_medium))
        cutflow.globals.setVal<float>("met_pt", nt.MET_pt());
    }
    return pass;
}


// Preselection functions
// Function to select VBF system mass >= 500
bool mjjGt500(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    cutflow.globals.setVal<double>("mjj", (cutflow.globals.getVal<LorentzVector>("ld_vbf_p4") + cutflow.globals.setVal<LorentzVector>("sd_vbf_p4")).M());
    bool logic = (cutflow.globals.getVal<double>("mjj")>=500);
    return logic;
}

// Function to select delta eta between VBF jets >= 3, then fill the selections
bool deltaEtaGt3(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    cutflow.globals.setVal<double>("detajj", TMath::Abs(cutflow.globals.getVal<LorentzVector>("ld_vbf_p4").eta()-cutflow.globals.getVal<LorentzVector>("sd_vbf_p4").eta()));
    bool logic = (cutflow.globals.getVal<double>("detajj")>=3);
}


// Function to calculate and select ST
bool stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    return false;    
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
        arbol.setLeaf<double>("lhe_mvvh", (W+H+Z).M());
        arbol.setLeaf<double>("lhe_lt", Z.Pt());
        arbol.setLeaf<double>("lhe_st", W.Pt()+H.Pt()+Z.Pt());
        arbol.setLeaf<double>("lhe_mjj", (jets[VBFindx1]+jets[VBFindx2]).M());
        arbol.setLeaf<double>("lhe_detajj", (TMath::Abs(jets[VBFindx1].Eta()-jets[VBFindx2].Eta())));

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