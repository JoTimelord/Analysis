// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"

// ====================================================================================================================================
// Namespace of functions for selecting only one fatjets
namespace ONEFATJETCUT{
    // void selectVBFS(int& i, int& j, std::vector<LorentzVector> jets);
    void initializeArbol(Arbol& arbol_); // Writes the necessary arbol branches such as mass of VVH system
    void initializeCutflow(Cutflow& cutflow_); // Holds the "unflattened" variables such as 'LorentzVector' of electrons
    void selectWjets(int& i, int& j, int VBF1, int VBF2, std::vector<LorentzVector> jets);
    void fillTree(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool fatJetHbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore);
    bool geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow);
}

// Namespace of functions for selecting two fatjets
namespace TWOFATJETSCUT{
    void initializeArbol(Arbol& arbol_);
    void initializeCutflow(Cutflow& cutflow_);
    void fillTree(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool geq2Fatjets(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool fatJetsScores(Nano& nt, Arbol& arbol, Cutflow& cutflow, double sethbbscore);
    bool geq2JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow);
    bool stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow);
}

// ====================================================================================================================================
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

// ************************************************************************************************************************************
// ====================================================================================================================================
// Function to write necessary cutflow variables (2 os leptons, 1 hbb fatjet, 4+ ak4 jets (2 VBF jets, 2 jets from W))
void ONEFATJETCUT::initializeCutflow(Cutflow& cutflow_) {
    cutflow_.globals.newVar<LorentzVector>("ld_lep_p4"); 
    cutflow_.globals.newVar<LorentzVector>("sd_lep_p4"); 
    cutflow_.globals.newVar<std::vector<int>>("jetidx");
    cutflow_.globals.newVar<int>("ld_lep_pdgid"); // leading lepton pdgid
    cutflow_.globals.newVar<int>("sd_lep_pdgid"); // subleading lepton pdgid
    cutflow_.globals.newVar<LorentzVector>("hbbjet_p4");
    cutflow_.globals.newVar<LorentzVector>("xbbjet_p4");
    // # of ak8 jets that pass the requirements: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
    cutflow_.globals.newVar<int>("n_ak8");
    cutflow_.globals.newVar<int>("n_ak4");
    cutflow_.globals.newVar<LorentzVector>("ld_wjet_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_wjet_p4");
    cutflow_.globals.newVar<LorentzVector>("ld_vbf_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_vbf_p4");
    cutflow_.globals.newVar<double>("mjj");
    cutflow_.globals.newVar<double>("detajj");
    cutflow_.globals.newVar<double>("met_pt");
    cutflow_.globals.newVar<double>("ST");
    cutflow_.globals.newVar<double>("LT");
    cutflow_.globals.newVar<bool>("medium_b_veto");
    cutflow_.globals.newVar<double>("hbb"); // highest hbb score (hbb score of hbb jet)
    cutflow_.globals.newVar<double>("xbb"); // highest xbb score (not necessarily that of hbb jet)
    cutflow_.globals.newVar<double>("pn_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<double>("sd_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<bool>("samejet"); // Check if highest xbb and highest hbb come from the same jet
}

// Function to write branches to root files
void ONEFATJETCUT::initializeArbol(Arbol& arbol_) {
    arbol_.newBranch<int>("event", -999);
    arbol_.newBranch<double>("xsec_sf", -999);
    // LHE variables
    arbol_.newBranch<double>("lhe_mvvh", -999);
    arbol_.newBranch<double>("lhe_lt", -999);
    arbol_.newBranch<double>("lhe_st", -999);
    arbol_.newBranch<double>("lhe_mjj", -999);
    arbol_.newBranch<double>("lhe_detajj", -999);
    // Reconstructed variables
    arbol_.newBranch<LorentzVector>("ld_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("mjj", -999.0);
    arbol_.newBranch<double>("detajj", -999.0);
    arbol_.newBranch<double>("dRjj", -999.0);
    arbol_.newBranch<LorentzVector>("ld_wjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_wjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("mqq", -999.0);
    arbol_.newBranch<double>("dRqq", -999.0);
    arbol_.newBranch<LorentzVector>("ld_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<int>("ld_lep_ID", -999); // pdgid for leading lepton
    arbol_.newBranch<int>("sd_lep_ID", -999);
    arbol_.newBranch<double>("mll", -999.0);
    arbol_.newBranch<double>("ptll", -999.0);
    arbol_.newBranch<double>("dRll", -999.0);
    arbol_.newBranch<double>("met", -999.0);
    arbol_.newBranch<int>("n_ak4", -999); // # of ak4 jets
    arbol_.newBranch<int>("n_ak8", -999);
    arbol_.newBranch<LorentzVector>("hbb_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("hbb_score", -999);
    arbol_.newBranch<double>("xbb_score", -999);
    arbol_.newBranch<double>("hbb_pnetmass", -999);
    arbol_.newBranch<double>("hbb_sdmass", -999);
    arbol_.newBranch<bool>("b_veto", false);
    arbol_.newBranch<double>("mvvh", -999.9);
    arbol_.newBranch<double>("st", -999.9);
    arbol_.newBranch<double>("lt", -999.9);
}

// Function to select two jets with minimum dR for W reconstruction
void ONEFATJETCUT::selectWjets(int& i, int& j, int VBF1, int VBF2, std::vector<LorentzVector> jets) {
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
void ONEFATJETCUT::fillTree(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    // Leptons
    LorentzVector ldlep_p4=cutflow.globals.getVal<LorentzVector>("ld_lep_p4");
    LorentzVector sdlep_p4=cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    arbol.setLeaf<LorentzVector>("ld_lep_LV",ldlep_p4);
    arbol.setLeaf<LorentzVector>("sd_lep_LV",sdlep_p4);
    arbol.setLeaf<int>("ld_lep_ID",cutflow.globals.getVal<int>("ld_lep_pdgid"));
    arbol.setLeaf<int>("sd_lep_ID",cutflow.globals.getVal<int>("sd_lep_pdgid"));
    arbol.setLeaf<double>("mll",(ldlep_p4+sdlep_p4).M());
    arbol.setLeaf<double>("ptll",(ldlep_p4+sdlep_p4).Pt());
    arbol.setLeaf<double>("dRll",fabs(ROOT::Math::VectorUtil::DeltaR(ldlep_p4,sdlep_p4)));
    // Fatjets
    arbol.setLeaf<int>("n_ak8", cutflow.globals.getVal<int>("n_ak8"));
    LorentzVector H=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    arbol.setLeaf<double>("hbb_score", cutflow.globals.getVal<double>("hbb"));
    arbol.setLeaf<LorentzVector>("hbb_LV", H);
    arbol.setLeaf<double>("xbb_score", cutflow.globals.getVal<double>("xbb"));
    arbol.setLeaf<double>("hbb_pnetmass", cutflow.globals.getVal<double>("pn_mass"));
    arbol.setLeaf<double>("hbb_sdmass", cutflow.globals.getVal<double>("sd_mass"));
    arbol.setLeaf<bool>("b_veto", cutflow.globals.getVal<bool>("medium_b_veto"));
    // ak4 jets
    LorentzVector ldwqqLV=cutflow.globals.getVal<LorentzVector>("ld_wjet_p4");
    LorentzVector sdwqqLV=cutflow.globals.getVal<LorentzVector>("sd_wjet_p4");
    LorentzVector ldvbfLV=cutflow.globals.getVal<LorentzVector>("ld_vbf_p4");
    LorentzVector sdvbfLV=cutflow.globals.getVal<LorentzVector>("sd_vbf_p4");
    arbol.setLeaf<LorentzVector>("ld_vbfjet_LV",ldvbfLV);
    arbol.setLeaf<LorentzVector>("sd_vbfjet_LV",sdvbfLV);
    arbol.setLeaf<LorentzVector>("ld_wjet_LV",ldwqqLV);
    arbol.setLeaf<LorentzVector>("sd_wjet_LV",sdwqqLV);
    arbol.setLeaf<double>("mjj", (ldvbfLV+sdvbfLV).M());
    arbol.setLeaf<double>("detajj", TMath::Abs(ldvbfLV.Eta()-sdvbfLV.Eta()));
    arbol.setLeaf<double>("dRjj",fabs(ROOT::Math::VectorUtil::DeltaR(ldvbfLV,sdvbfLV)));
    arbol.setLeaf<double>("mqq", (ldwqqLV+sdwqqLV).M());
    arbol.setLeaf<double>("dRqq",fabs(ROOT::Math::VectorUtil::DeltaR(ldwqqLV,sdwqqLV)));    
    arbol.setLeaf<int>("n_ak4", cutflow.globals.getVal<int>("n_ak4"));
    // Recombined 
    LorentzVector Z=ldlep_p4+sdlep_p4;
    LorentzVector W=ldwqqLV+sdwqqLV;
    arbol.setLeaf<double>("met", cutflow.globals.getVal<double>("met_pt"));
    arbol.setLeaf<double>("mvvh", (W+H+Z).M());
    arbol.setLeaf<double>("lt", Z.Pt()+cutflow.globals.getVal<double>("met_pt"));
    arbol.setLeaf<double>("st", Z.Pt()+H.Pt()+W.Pt());
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
        double this_pt=nt.Electron_pt().at(elec_i);
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
    bool logic=(n_leptons==2 && (elec.size()==2 || muon.size()==2) && os<0);
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

// Function to select events with at least one fatjet (boosted Higgs); check overlap against leptons
// Requirement on fatjet: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
bool ONEFATJETCUT::geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_fatjet=0;
    double highesthbb=-999.999; // mass-correlated score
    double highestxbb=-999.999; // mass-decorrelated score
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
        cutflow.globals.setVal<double>("sd_mass",pnetmass);
    }
    return logic;
}

// Function to make checks on fatjet score
bool ONEFATJETCUT::fatJetHbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore=0.5) {
    double hbbscore=cutflow.globals.getVal<double>("hbb") ;
    bool logic = (hbbscore >= setscore);
    return logic;
}

// Function to select events with at least four ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against fatjets
// Select VBF jets
bool ONEFATJETCUT::geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
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
    int ldlep_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(0);
    int sdlep_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(1);
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
        if (ldlep_jetid == -1 || sdlep_jetid == -1) {
            if (ROOT::Math::VectorUtil::DeltaR(ld_lep_p4, jet_p4) < 0.4 || ROOT::Math::VectorUtil::DeltaR(sd_lep_p4, jet_p4) < 0.4)
            {
                isOverlap = true;
            }
        }
        else if (ldlep_jetid == (int) ijet || sdlep_jetid == (int) ijet)
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
        cutflow.globals.setVal<LorentzVector>("ld_wjet_p4", nt.Jet_p4().at(w1));
        cutflow.globals.setVal<LorentzVector>("sd_wjet_p4", nt.Jet_p4().at(w2));
        cutflow.globals.setVal<bool>("medium_b_veto", bvetos);
        cutflow.globals.setVal<double>("met_pt", nt.MET_pt());
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
    return logic;
}

// Function to calculate and select ST
bool ONEFATJETCUT::stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    LorentzVector Hbb=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    LorentzVector wjj=cutflow.globals.getVal<LorentzVector>("ld_wjet_p4")+cutflow.globals.getVal<LorentzVector>("sd_wjet_p4");
    LorentzVector zll=cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    cutflow.globals.setVal<double>("ST", (Hbb+wjj+zll).Pt());
    bool logic=cutflow.globals.getVal<double>("ST")>=950;
    if (logic) {ONEFATJETCUT::fillTree(nt, arbol, cutflow);}
    return logic; 
}

// ************************************************************************************************************************************
// ====================================================================================================================================
// Cutflow variables: 2 os leptons, 1 hbb fatjet, 1 boosted V jet, 2+ ak4 jets (2 VBF jets)
void TWOFATJETSCUT::initializeCutflow(Cutflow& cutflow_) {
    cutflow_.globals.newVar<LorentzVector>("ld_lep_p4"); 
    cutflow_.globals.newVar<LorentzVector>("sd_lep_p4"); 
    cutflow_.globals.newVar<std::vector<int>>("jetidx");
    cutflow_.globals.newVar<int>("ld_lep_pdgid"); // leading lepton pdgid
    cutflow_.globals.newVar<int>("sd_lep_pdgid"); // subleading lepton pdgid
    cutflow_.globals.newVar<bool>("medium_b_veto");
    cutflow_.globals.newVar<LorentzVector>("hbbjet_p4");
    cutflow_.globals.newVar<LorentzVector>("xbbjet_p4");
    cutflow_.globals.newVar<LorentzVector>("xvqqjet_p4"); // p4 of boosted V jet
    cutflow_.globals.newVar<int>("n_ak8");
    cutflow_.globals.newVar<int>("n_ak4");
    cutflow_.globals.newVar<LorentzVector>("ld_vbf_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_vbf_p4");
    cutflow_.globals.newVar<double>("hbb"); // highest hbb score (hbb score of hbb jet)
    cutflow_.globals.newVar<double>("xbb"); // highest xbb score (not necessarily that of hbb jet)
    cutflow_.globals.newVar<double>("hbb_pn_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<double>("hbb_sd_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<double>("xvqq"); // highest xvqq score (hbb jet excluded)
    cutflow_.globals.newVar<double>("xvqq_pn_mass"); 
    cutflow_.globals.newVar<double>("xvqq_sd_mass"); 
    cutflow_.globals.newVar<double>("met_pt");
    cutflow_.globals.newVar<double>("mjj");
    cutflow_.globals.newVar<double>("detajj");
    cutflow_.globals.newVar<double>("ST");
    cutflow_.globals.newVar<double>("LT");
}

// Function to write branches to root files
void TWOFATJETSCUT::initializeArbol(Arbol& arbol_) {
    arbol_.newBranch<int>("event", -999);
    arbol_.newBranch<double>("xsec_sf", -999.0);
    // LHE variables
    arbol_.newBranch<double>("lhe_mvvh", -999.0);
    arbol_.newBranch<double>("lhe_lt", -999.0);
    arbol_.newBranch<double>("lhe_st", -999.0);
    arbol_.newBranch<double>("lhe_mjj", -999.0);
    arbol_.newBranch<double>("lhe_detajj", -999.0);
    // Reconstructed variables
    arbol_.newBranch<LorentzVector>("ld_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("mjj", -999.9);
    arbol_.newBranch<double>("detajj", -999.9);
    arbol_.newBranch<double>("dRjj", -999.9);
    arbol_.newBranch<LorentzVector>("ld_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<int>("ld_lep_ID", -999);
    arbol_.newBranch<int>("sd_lep_ID", -999);
    arbol_.newBranch<double>("mll", -999);
    arbol_.newBranch<double>("ptll", -999.0);
    arbol_.newBranch<double>("dRll", -999.0);
    arbol_.newBranch<double>("met", -999.0);
    arbol_.newBranch<int>("n_ak4", -999); // # of ak4 jets
    arbol_.newBranch<int>("n_ak8", -999);
    arbol_.newBranch<double>("hbb_score", -999.0);
    arbol_.newBranch<LorentzVector>("hbb_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("xbb_score", -999.0);
    arbol_.newBranch<double>("hbb_pnetmass", -999.0);
    arbol_.newBranch<double>("hbb_sdmass", -999.0);
    arbol_.newBranch<double>("xvqq_score", -999.0);
    arbol_.newBranch<double>("xvqq_sdmass", -999.0);
    arbol_.newBranch<double>("xvqq_pnetmass", -999.0);
    arbol_.newBranch<LorentzVector>("xvqq_LV", {-999,-999,-999,-999});
    arbol_.newBranch<bool>("b_veto", false);
    arbol_.newBranch<double>("st", -999.9);
    arbol_.newBranch<double>("mvvh", -999.9);
    arbol_.newBranch<double>("lt", -999.9);
}

// Function to fill arbol trees
void TWOFATJETSCUT::fillTree(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    // Leptons
    LorentzVector ldlep_p4=cutflow.globals.getVal<LorentzVector>("ld_lep_p4");
    LorentzVector sdlep_p4=cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    arbol.setLeaf<LorentzVector>("ld_lep_LV",ldlep_p4);
    arbol.setLeaf<LorentzVector>("sd_lep_LV",sdlep_p4);
    arbol.setLeaf<int>("ld_lep_ID",cutflow.globals.getVal<int>("ld_lep_pdgid"));
    arbol.setLeaf<int>("sd_lep_ID",cutflow.globals.getVal<int>("sd_lep_pdgid"));
    arbol.setLeaf<double>("mll",(ldlep_p4+sdlep_p4).M());
    arbol.setLeaf<double>("ptll",(ldlep_p4+sdlep_p4).Pt());
    arbol.setLeaf<double>("dRll",fabs(ROOT::Math::VectorUtil::DeltaR(ldlep_p4,sdlep_p4)));
    // Hbb fatjet
    arbol.setLeaf<int>("n_ak8", cutflow.globals.getVal<int>("n_ak8"));
    LorentzVector H=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    arbol.setLeaf<double>("hbb_score", cutflow.globals.getVal<double>("hbb"));
    arbol.setLeaf<LorentzVector>("hbb_LV", H);
    arbol.setLeaf<double>("xbb_score", cutflow.globals.getVal<double>("xbb"));
    arbol.setLeaf<double>("hbb_pnetmass", cutflow.globals.getVal<double>("hbb_pn_mass"));
    arbol.setLeaf<double>("hbb_sdmass", cutflow.globals.getVal<double>("hbb_sd_mass"));
    arbol.setLeaf<bool>("b_veto", cutflow.globals.getVal<bool>("medium_b_veto"));
    // Vqq fatjet
    LorentzVector V=cutflow.globals.getVal<LorentzVector>("xvqqjet_p4");
    arbol.setLeaf<LorentzVector>("xvqq_LV", V);
    arbol.setLeaf<double>("xvqq_score", cutflow.globals.getVal<double>("xvqq"));
    arbol.setLeaf<double>("xvqq_sdmass", cutflow.globals.getVal<double>("xvqq_sd_mass"));
    arbol.setLeaf<double>("xvqq_pnetmass", cutflow.globals.getVal<double>("xvqq_pn_mass"));
    // ak4 jets
    LorentzVector ldvbfLV=cutflow.globals.getVal<LorentzVector>("ld_vbf_p4");
    LorentzVector sdvbfLV=cutflow.globals.getVal<LorentzVector>("sd_vbf_p4");
    arbol.setLeaf<LorentzVector>("ld_vbfjet_LV",ldvbfLV);
    arbol.setLeaf<LorentzVector>("sd_vbfjet_LV",sdvbfLV);
    arbol.setLeaf<double>("mjj", (ldvbfLV+sdvbfLV).M());
    arbol.setLeaf<double>("detajj", TMath::Abs(ldvbfLV.Eta()-sdvbfLV.Eta()));
    arbol.setLeaf<double>("dRjj",fabs(ROOT::Math::VectorUtil::DeltaR(ldvbfLV,sdvbfLV)));
    arbol.setLeaf<int>("n_ak4", cutflow.globals.getVal<int>("n_ak4"));
    // Recombined
    double met_pt=cutflow.globals.getVal<double>("met_pt");
    LorentzVector Z=ldlep_p4+sdlep_p4;
    arbol.setLeaf<double>("met", met_pt);
    arbol.setLeaf<double>("mvvh", (V+H+Z).M());
    arbol.setLeaf<double>("lt", Z.Pt()+met_pt);
    arbol.setLeaf<double>("st", Z.Pt()+H.Pt()+V.Pt()+met_pt);
}

// Function to select at least two fatjets, one boosted Higgs and one boosted V
// check overlap against leptons
// Requirement on fatjet: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
bool TWOFATJETSCUT::geq2Fatjets(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    int n_fatjet=0;
    double highesthbb=-999.999; // mass-correlated score
    double highestxbb=-999.999; // mass-decorrelated score
    double highestxvqq=-999.999;
    int hbbindx=-999;
    int xbbindx=-999;
    int xvqqindx=-999;
    double hbb_pnetmass=-999.999;
    double hbb_sdmass=-999.9;
    double xvqq_pnetmass=-999.99;
    double xvqq_sdmass=-999.99;
    std::vector<int> fatjetindx={};
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
            fatjetindx.push_back(ifatjet);
            double this_xbb=(nt.FatJet_particleNetMD_Xbb()[ifatjet]/(nt.FatJet_particleNetMD_Xbb()[ifatjet]+nt.FatJet_particleNetMD_QCD()[ifatjet]));
            // Select boosted higgs by hbb
            if (nt.FatJet_particleNet_HbbvsQCD()[ifatjet] > highesthbb) { 
                hbbindx = ifatjet; 
                highesthbb = nt.FatJet_particleNet_HbbvsQCD()[ifatjet];
                hbb_pnetmass = nt.FatJet_particleNet_mass()[ifatjet];
                hbb_sdmass = m_softdrop;
            }
            // Select boosted Higgs by xbb
            if (this_xbb > highestxbb) {
                xbbindx = ifatjet;
                highestxbb = this_xbb;
            }
        }
    }
    // Selected boosted V by xbb+xcc+xqq
    for (int j=0; j<fatjetindx.size(); j++) {
        int this_indx=fatjetindx.at(j);
        double this_xvqq=(nt.FatJet_particleNetMD_Xbb()[this_indx]+nt.FatJet_particleNetMD_Xcc()[this_indx]+nt.FatJet_particleNetMD_Xqq()[this_indx])/(nt.FatJet_particleNetMD_Xbb()[this_indx]+nt.FatJet_particleNetMD_Xcc()[this_indx]+nt.FatJet_particleNetMD_Xqq()[this_indx]+nt.FatJet_particleNet_QCD()[this_indx]);
        if (this_xvqq > highestxvqq) {
            xvqqindx=this_indx;
            highestxvqq=this_xvqq;
            xvqq_pnetmass=nt.FatJet_particleNet_mass()[this_indx];
            xvqq_sdmass=nt.FatJet_msoftdrop()[this_indx];
        }
    }
    // Select Logic: Pass events with at least two fatjets
    bool logic = (n_fatjet >= 2);
    if (logic) {
        cutflow.globals.setVal<LorentzVector>("hbbjet_p4", nt.FatJet_p4()[hbbindx]);
        cutflow.globals.setVal<LorentzVector>("xbbjet_p4", nt.FatJet_p4()[xbbindx]);
        cutflow.globals.setVal<LorentzVector>("xvqqjet_p4", nt.FatJet_p4()[xvqqindx]);
        cutflow.globals.setVal<double>("hbb",highesthbb);
        cutflow.globals.setVal<double>("xbb",highestxbb);
        cutflow.globals.setVal<double>("xvqq",highestxvqq);
        cutflow.globals.setVal<int>("n_ak8",n_fatjet);
        cutflow.globals.setVal<double>("hbb_pn_mass",hbb_pnetmass);
        cutflow.globals.setVal<double>("hbb_sd_mass",hbb_sdmass);
        cutflow.globals.setVal<double>("xvqq_sd_mass",xvqq_sdmass);
        cutflow.globals.setVal<double>("xvqq_pn_mass",xvqq_pnetmass);
    }
    return logic;
}

// Function to make checks on fatjet score
bool TWOFATJETSCUT::fatJetsScores(Nano& nt, Arbol& arbol, Cutflow& cutflow, double sethbbscore=0.5) {
    double hbbscore=cutflow.globals.getVal<double>("hbb") ;
    bool logic = (hbbscore >= sethbbscore);
    return logic;
}

// Function to select events with at least four ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against two fatjets
// Select VBF jets
bool TWOFATJETSCUT::geq2JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    std::vector<LorentzVector> jets{};
    int n_jets=0;
    int VBF1=-999;
    int VBF2=-999;
    bool bvetos=false;
    bool pass=false;
    LorentzVector ld_lep_p4=cutflow.globals.getVal<LorentzVector>("ld_lep_p4");
    LorentzVector sd_lep_p4=cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    LorentzVector hbb_p4=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    LorentzVector xvqq_p4=cutflow.globals.getVal<LorentzVector>("xvqqjet_p4");
    int ldlep_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(0);
    int sdlep_jetid=cutflow.globals.getVal<std::vector<int>>("jetidx").at(1);
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
        if (ldlep_jetid == -1 || sdlep_jetid == -1) {
            if (ROOT::Math::VectorUtil::DeltaR(ld_lep_p4, jet_p4) < 0.4 || ROOT::Math::VectorUtil::DeltaR(sd_lep_p4, jet_p4) < 0.4)
            {
                isOverlap = true;
            }
        }
        else if (ldlep_jetid == (int) ijet || sdlep_jetid == (int) ijet)
        {
            isOverlap = true;
        }
        // Overlap check against fatjet
        if (ROOT::Math::VectorUtil::DeltaR(hbb_p4, jet_p4) <= 0.8 || ROOT::Math::VectorUtil::DeltaR(xvqq_p4, jet_p4) <= 0.8) {
            isOverlap = true;
        }
        if (isOverlap) {continue; }
        if (nt.Jet_pt().at(ijet) < 30.0) {continue; }
        n_jets++;
        bvetos=(nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_medium);
        jets.push_back(jet_p4);
    }
    pass=(n_jets>=2);
    if (pass) {
        selectVBFs(VBF1, VBF2, jets);
        cutflow.globals.setVal<int>("n_ak4", n_jets);
        cutflow.globals.setVal<LorentzVector>("ld_vbf_p4", nt.Jet_p4().at(VBF1));
        cutflow.globals.setVal<LorentzVector>("sd_vbf_p4", nt.Jet_p4().at(VBF2));
        cutflow.globals.setVal<bool>("medium_b_veto", bvetos);
        cutflow.globals.setVal<double>("met_pt", nt.MET_pt());
        cutflow.globals.setVal<double>("mjj", (nt.Jet_p4().at(VBF1) + nt.Jet_p4().at(VBF2)).M());
        cutflow.globals.setVal<double>("detajj", fabs(nt.Jet_eta().at(VBF1)-nt.Jet_eta().at(VBF2)));
    }
    return pass;
}

// Function to select events with ST >= 950
bool TWOFATJETSCUT::stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    double met=cutflow.globals.getVal<double>("met_pt");
    LorentzVector Hbb=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    LorentzVector Vqq=cutflow.globals.getVal<LorentzVector>("xvqqjet_p4");
    LorentzVector Zll=cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    cutflow.globals.setVal<double>("ST", Hbb.Pt()+Vqq.Pt()+Zll.Pt()+met);
    bool logic=cutflow.globals.getVal<double>("ST")>=950;
    if (logic) {TWOFATJETSCUT::fillTree(nt, arbol, cutflow);}
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