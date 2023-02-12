// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

// Function to select VBFjets with the largest delta Eta
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

// Function to select two jets with minimum dR for W reconstruction
void selectWjets(int& i, int& j, int VBF1, int VBF2, std::vector<LorentzVector> jets) {
    double mindR=999;
    for (int ijet=0; ijet<jets.size(); ijet++) {
        if (ijet!=VBF1 && ijet!=VBF2) {
            for (int jjet=ijet+1; jjet<jets.size(); jjet++) {
                if (jjet!=VBF1 && jjet!=VBF2) {
                    double this_deltaR=fabs(ROOT::Math::VectorUtil::DeltaR(jets.at(ijet),jets.at(jjet)));
                    if (this_deltaR < mindR) {
                        mindR=this_deltaR;
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
    LorentzVector H=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    LorentzVector W=cutflow.globals.getVal<LorentzVector>("ld_jet_p4")+cutflow.globals.getVal<LorentzVector>("sd_jet_p4");
    arbol.setLeaf<LorentzVector>("vbfjet1_LV",cutflow.globals.getVal<LorentzVector>("ld_vbf_p4"));
    arbol.setLeaf<LorentzVector>("vbfjet2_LV",cutflow.globals.getVal<LorentzVector>("sd_vbf_p4"));
    arbol.setLeaf<LorentzVector>("wjet1_LV",cutflow.globals.getVal<LorentzVector>("ld_jet_p4"));
    arbol.setLeaf<LorentzVector>("wjet2_LV",cutflow.globals.getVal<LorentzVector>("sd_jet_p4"));
    arbol.setLeaf<LorentzVector>("lep1_LV",cutflow.globals.getVal<LorentzVector>("ld_lep_p4"));
    arbol.setLeaf<LorentzVector>("lep2_LV",cutflow.globals.getVal<LorentzVector>("sd_lep_p4"));
    arbol.setLeaf<int>("lep1_ID",cutflow.globals.getVal<int>("ldid"));
    arbol.setLeaf<int>("lep2_ID",cutflow.globals.getVal<int>("sdid"));
    arbol.setLeaf<double>("mll",(cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4")).M());
    arbol.setLeaf<double>("ptll",(cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4")).Pt());
    arbol.setLeaf<double>("met", cutflow.globals.getVal<float>("met_pt"));
    arbol.setLeaf<int>("n_ak4", cutflow.globals.getVal<int>("n_ak4"));
    arbol.setLeaf<int>("n_ak8", cutflow.globals.getVal<int>("n_ak8"));
    arbol.setLeaf<double>("hbb_score", cutflow.globals.getVal<double>("hbb"));
    arbol.setLeaf<LorentzVector>("hbb_LV", cutflow.globals.getVal<LorentzVector>("hbbjet_p4"));
    arbol.setLeaf<double>("xbb_score", cutflow.globals.getVal<double>("xbb"));
    arbol.setLeaf<double>("hbb_pnetmass", cutflow.globals.getVal<double>("pn_mass"));
    arbol.setLeaf<bool>("b_veto", cutflow.globals.getVal<bool>("medium_b_veto"));
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
    std::vector<LorentzVector> elec{};
    std::vector<LorentzVector> muon{};
    std::vector<int> lepindx{};
    std::vector<int> pdgids{}; // Save to cutflow.globals
    int os=1;
    // Select electrons
    for (unsigned int elec_i = 0; elec_i < nt.Electron_pt().size(); ++elec_i)
    {
        float this_pt=nt.Electron_pt().at(elec_i);
        LorentzVector this_elec=nt.Electron_p4().at(elec_i);
        int this_elec_id=nt.Electron_jetIdx()[elec_i];
        int this_pdgid=nt.Electron_pdgId()[elec_i];
        // Use IDtight so that the requirement is for analysis
        if (this_pt > 30 && ttH_UL::electronID(elec_i, ttH::IDtight, nt.year())) { 
        // if (pt > 30) {
            n_leptons++; 
            elec.push_back(this_elec);
            lepindx.push_back(this_elec_id);
            os=os*this_pdgid;
            pdgids.push_back(this_pdgid);
        }
    }
    // Select Muons
    for (unsigned int muon_i = 0; muon_i < nt.Muon_pt().size(); ++muon_i) {
        float this_pt=nt.Muon_pt().at(muon_i);
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
    bool logic=(n_leptons==2 && (elec.size()==2 || muon.size()==2) && os<0);
    if (logic) {
        LorentzVector ld, sd;
        int ldid=-999;
        int sdid=-999;
        // Select leading and subleading leptons with their pdgids
        if (elec.size()==3) {
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
        cutflow.globals.setVal<std::vector<int>>("jetidx", lepindx);
        cutflow.globals.setVal<int>("ldid",ldid);
        cutflow.globals.setVal<int>("sdid",sdid);
    }
    
    // Require only two lepton events to pass
    return logic;
}

// Function to select events with at least one fatjet (boosted Higgs); check overlap against leptons
// Requirement on fatjet: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
bool geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_fatjet=0;
    float highesthbb=-999.999; // mass-correlated score
    float highestxbb=-999.999; // mass-decorrelated score
    int hbbindx=-999;
    int xbbindx=-999;
    double pnetmass=-999.999;
    LorentzVector sd_lep_p4=cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    LorentzVector ld_lep_p4=cutflow.globals.getVal<LorentzVector>("ld_lep_p4");
    for (unsigned int ifatjet = 0; ifatjet < nt.nFatJet(); ++ifatjet) {
        LorentzVector fatjet_p4=nt.FatJet_p4()[ifatjet];
        double m_softdrop=nt.FatJet_msoftdrop()[ifatjet];
        double jetid=nt.FatJet_jetId()[ifatjet];
        // Check against overlap with leptons by checking difference in pseudorapidity (Eta) and Phi
        bool isOverlap = false;
        if (ROOT::Math::VectorUtil::DeltaR(sd_lep_p4,fatjet_p4) < 0.8) {
            isOverlap = true;
        }
        if (ROOT::Math::VectorUtil::DeltaR(ld_lep_p4,fatjet_p4) < 0.8) {
            isOverlap = true;
        }
        if (isOverlap) {continue; }
        // Check ak8 jets with softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
        if (m_softdrop>=40 && fatjet_p4.Pt()>250 && fabs(fatjet_p4.Eta())<2.5 && jetid>0) {
            n_fatjet++;
            double this_xbb=(nt.FatJet_particleNetMD_Xbb()[ifatjet]/(nt.FatJet_particleNetMD_Xbb()[ifatjet]+nt.FatJet_particleNetMD_QCD()[ifatjet]));
            // Select boosted higgs by hbb
            if (nt.FatJet_particleNet_HbbvsQCD()[ifatjet] > highesthbb) { 
                hbbindx = ifatjet; 
                highesthbb = nt.FatJet_particleNet_HbbvsQCD()[ifatjet];
                pnetmass = nt.FatJet_particleNet_mass()[ifatjet];
            }
            // Select boosted Higgs by xbb
            if (this_xbb > highestxbb) {
                xbbindx = ifatjet;
                highestxbb = this_xbb;
            }
        }
    }
    // Select Logic: Pass events with at least one fatjets
    bool logic = (n_fatjet >= 1);
    if (logic) {
        cutflow.globals.setVal<LorentzVector>("hbbjet_p4", nt.FatJet_p4()[hbbindx]);
        cutflow.globals.setVal<LorentzVector>("xbbjet_p4", nt.FatJet_p4()[xbbindx]);
        cutflow.globals.setVal<double>("hbb",highesthbb);
        cutflow.globals.setVal<double>("xbb",highestxbb);
        cutflow.globals.setVal<int>("n_ak8",n_fatjet);
        cutflow.globals.setVal<double>("pn_mass",pnetmass);
    }
    return logic;
}


// Function to make checks on fatjet score
bool fatJetHbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore=0.5) {
    double hbbscore=cutflow.globals.getVal<double>("hbb") ;
    bool logic = (hbbscore >= setscore);
    return logic;
}

// Function to select events with at least four ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against fatjets
// Select VBF jets
bool geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    std::vector<LorentzVector> jets{};
    int n_jets=0;
    int VBF1=-999;
    int VBF2=-999;
    int w1=-999;
    int w2=-999;
    bool bvetos=false;
    bool pass=false;
    LorentzVector ld_lep_p4=cutflow.globals.getVal<LorentzVector>("ld_lep_p4");
    LorentzVector sd_lep_p4=cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    LorentzVector hbb_p4=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    int ld_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(0);
    int sd_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(1);
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
        LorentzVector jet_p4=nt.Jet_p4().at(ijet);
        if (ld_jetid == -1 || sd_jetid == -1) {
            if (ROOT::Math::VectorUtil::DeltaR(ld_lep_p4, jet_p4) < 0.4 || ROOT::Math::VectorUtil::DeltaR(sd_lep_p4, jet_p4) < 0.4)
            {
                isOverlap = true;
            }
        }
        else if (ld_jetid == (int) ijet || sd_jetid == (int) ijet)
        {
            isOverlap = true;
        }
        // Overlap check against fatjet
        if (ROOT::Math::VectorUtil::DeltaR(hbb_p4, jet_p4) <= 0.8) {
            isOverlap = true;
        }
         
        if (isOverlap) {continue; }
        if (nt.Jet_pt().at(ijet) < 30.0) {continue; }
        n_jets++;
        bvetos=(nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_medium);
        jets.push_back(jet_p4);
    }
    pass=(n_jets>=4);
    if (pass) {
        selectVBFs(VBF1, VBF2, jets);
        selectWjets(w1, w2, VBF1, VBF2, jets);
        cutflow.globals.setVal<int>("n_ak4", n_jets);
        cutflow.globals.setVal<LorentzVector>("ld_vbf_p4", nt.Jet_p4().at(VBF1));
        cutflow.globals.setVal<LorentzVector>("sd_vbf_p4", nt.Jet_p4().at(VBF2));
        cutflow.globals.setVal<LorentzVector>("ld_jet_p4", nt.Jet_p4().at(w1));
        cutflow.globals.setVal<LorentzVector>("sd_jet_p4", nt.Jet_p4().at(w2));
        cutflow.globals.setVal<bool>("medium_b_veto", bvetos);
        cutflow.globals.setVal<float>("met_pt", nt.MET_pt());
        cutflow.globals.setVal<double>("mjj", (nt.Jet_p4().at(VBF1) + nt.Jet_p4().at(VBF2)).M());
        cutflow.globals.setVal<double>("detajj", fabs(nt.Jet_eta().at(VBF1)-nt.Jet_eta().at(VBF2)));

    }
    
    return pass;
}


// Preselection functions
// Function to select VBF system mass >= 500
bool mjjGt500(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    bool logic = (cutflow.globals.getVal<double>("mjj")>=500);
    return logic;
}

// Function to select delta eta between VBF jets >= 3, then fill the selections
bool deltaEtaGt3(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    bool logic = (cutflow.globals.getVal<double>("detajj")>=3);
    if (logic) {
        LorentzVector Hbb=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
        LorentzVector wjj=cutflow.globals.getVal<LorentzVector>("ld_jet_p4")+cutflow.globals.getVal<LorentzVector>("sd_jet_p4");
        LorentzVector zll=cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
        cutflow.globals.setVal<float>("ST", (Hbb+wjj+zll).Pt());
    }
    return logic;
}


// Function to calculate and select ST
bool stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    bool logic=cutflow.globals.getVal<float>("ST")>=950;
    return logic; 
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