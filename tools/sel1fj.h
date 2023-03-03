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
    bool geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow, double sdmass, double pt, double deta);
    bool HbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore);
    bool geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow, double pt);
    bool stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow);
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
    cutflow_.globals.newVar<double>("hbb"); // highest hbb score (hbb score of hbb jet)
    cutflow_.globals.newVar<double>("xbb"); // highest xbb score (not necessarily that of hbb jet)
    cutflow_.globals.newVar<double>("pn_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<double>("sd_mass"); // particle net mass of hbb jet
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
    arbol_.newBranch<LorentzVector>("ld_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_lep_LV", {-999,-999,-999,-999});
    arbol_.newBranch<int>("ld_lep_ID", -999); // pdgid for leading lepton
    arbol_.newBranch<int>("sd_lep_ID", -999);
    arbol_.newBranch<double>("mll", -999.0);
    arbol_.newBranch<double>("ptll", -999.0);
    arbol_.newBranch<double>("dRll", -999.0);
    arbol_.newBranch<int>("n_ak8", -999);
    arbol_.newBranch<LorentzVector>("hbb_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("hbb_score", -999);
    arbol_.newBranch<double>("xbb_score", -999);
    arbol_.newBranch<double>("hbb_pnetmass", -999);
    arbol_.newBranch<double>("hbb_sdmass", -999);
    arbol_.newBranch<LorentzVector>("ld_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_vbfjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("mjj", -999.0);
    arbol_.newBranch<double>("detajj", -999.0);
    arbol_.newBranch<double>("dRjj", -999.0);
    arbol_.newBranch<LorentzVector>("ld_wjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("sd_wjet_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("mqq", -999.0);
    arbol_.newBranch<double>("dRqq", -999.0);
    arbol_.newBranch<double>("met", -999.0);
    arbol_.newBranch<int>("n_ak4", -999); // # of ak4 jets
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


// Function to select events with at least one fatjet (boosted Higgs); check overlap against leptons
// Requirement on fatjet: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
bool ONEFATJETCUT::geq1FatjetsMassGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow, double sdmass=40, double pt=250, double deta=2.5) {
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
        if (ROOT::Math::VectorUtil::DeltaR(sd_lep_p4,fatjet_p4) < 0.8 || ROOT::Math::VectorUtil::DeltaR(ld_lep_p4,fatjet_p4) < 0.8) {
            isOverlap = true;
        }
        if (isOverlap) {continue; }
        // Check ak8 jets with softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
        if (m_softdrop>=sdmass && fatjet_p4.Pt()>pt && fabs(fatjet_p4.Eta())<deta && jetid>0) {
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
bool ONEFATJETCUT::HbbScore(Nano& nt, Arbol& arbol, Cutflow& cutflow, double setscore=0.5) {
    double hbbscore=cutflow.globals.getVal<double>("hbb") ;
    bool logic = (hbbscore >= setscore);
    return logic;
}

// Function to select events with at least four ak4 jets with Pt > 30; check overlap against leptons;
// check overlap against fatjets
// Select VBF jets
bool ONEFATJETCUT::geq4JetsPtGt30(Nano& nt, Arbol& arbol, Cutflow& cutflow, double pt=30) {
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
        if (nt.Jet_pt().at(ijet) < pt) {continue; }
        n_jets++;
        bvetos=(nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_medium);
        jets.push_back(jet_p4);
    }

    pass=(n_jets>=4);
    if (pass) {
        selectVBFs(VBF1, VBF2, jets);
        selectWjets(w1, w2, VBF1, VBF2, jets);
        cutflow.globals.setVal<int>("n_ak4", n_jets);
        cutflow.globals.setVal<LorentzVector>("ld_vbf_p4", jets.at(VBF1));
        cutflow.globals.setVal<LorentzVector>("sd_vbf_p4", jets.at(VBF2));
        cutflow.globals.setVal<LorentzVector>("ld_wjet_p4", jets.at(w1));
        cutflow.globals.setVal<LorentzVector>("sd_wjet_p4", jets.at(w2));
        cutflow.globals.setVal<bool>("medium_b_veto", bvetos);
        cutflow.globals.setVal<double>("met_pt", nt.MET_pt());
        cutflow.globals.setVal<double>("mjj", (jets.at(VBF1) + jets.at(VBF2)).M());
        cutflow.globals.setVal<double>("detajj", fabs((jets.at(VBF1)).Eta()-(jets.at(VBF2)).Eta()));
    }
    return pass;
}


// Function to calculate and select ST
bool ONEFATJETCUT::stgeq950(Nano& nt, Arbol& arbol, Cutflow& cutflow) {
    LorentzVector Hbb=cutflow.globals.getVal<LorentzVector>("hbbjet_p4");
    LorentzVector wjj=cutflow.globals.getVal<LorentzVector>("ld_wjet_p4")+cutflow.globals.getVal<LorentzVector>("sd_wjet_p4");
    LorentzVector zll=cutflow.globals.getVal<LorentzVector>("ld_lep_p4")+cutflow.globals.getVal<LorentzVector>("sd_lep_p4");
    cutflow.globals.setVal<double>("ST", Hbb.Pt()+wjj.Pt()+zll.Pt());
    bool logic=cutflow.globals.getVal<double>("ST")>=950;
    if (logic) {ONEFATJETCUT::fillTree(nt, arbol, cutflow);}
    return logic; 
}

