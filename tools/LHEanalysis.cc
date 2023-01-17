#include "LHEanalysis.h"


void LHEAnalysis::Observable::selectObservable() {
    // How many particles are in each event
    unsigned int particleNo=nt.LHEPart_status().size();
    std::vector<LV> jets;
    Int_t VBFindx1, VBFindx2;
    Float_t eta1=-1000;
    Float_t eta2=1000;
    for (unsigned int i=0; i<particleNo; i++) {
        // Check if it's outgoing particles
        if (nt.LHEPart_status()[i]==1) {
            Int_t pdgid=nt.LHEPart_pdgId()[i];
            TLorentzVector p4;
            p4.SetPtEtaPhiM(nt.LHEPart_pt()[i],nt.LHEPart_eta()[i],nt.LHEPart_phi()[i],nt.LHEPart_mass()[i]);
            switch (pdgid) {
                case 24:
                    W_=RooUtil::Calc::getLV(p4);
                    break;
                case 25:
                    H_=RooUtil::Calc::getLV(p4);
                    break;
                case 23:
                    Z_=RooUtil::Calc::getLV(p4);
                    break;
                default:
                    jets.push_back(RooUtil::Calc::getLV(p4));
                    break;
            } 
        }
    }
    // Select the VBF jets with the largest etas in both directions
    for (unsigned int j=0; j<jets.size(); j++) {
        if (jets[j].Eta()>0 && jets[j].Eta()>eta1) {
            VBFindx1=j;
        }
        else if (jets[j].Eta()<0 && jets[j].Eta()<eta2) {
            VBFindx2=j;
        }
    }
    VBF_.push_back(jets[VBFindx1]);
    VBF_.push_back(jets[VBFindx2]);
}

LHEAnalysis::Histogram::Histogram () {
    MVVH_=new TH1F("MVVH", "Invariant mass of the VVH system", 1080, 0, 3500);
    LT_=new TH1F("LT", "Transverse Momentum of the Z boson", 1080, 0, 1500);
    ST_=new TH1F("ST", "Transverse Momentum of the VVH", 1080, 0, 4500);
    MJJ_=new TH1F("MJJ", "Invariant Mass of the VBF quarks system", 1080, 0, 3500);
    DETAJJ_=new TH1F("DETAJJ", "Delta Eta of the two VBF quarks", 1080, 0, 10);
}

void LHEAnalysis::Histogram::fillHistogram(LHEAnalysis::Observable& obs, Int_t wgt) {
    MVVH_->Fill((obs.W_+obs.Z_+obs.H_).M(), wgt);
    LT_->Fill(obs.Z_.Pt(), wgt);
    ST_->Fill((obs.W_+obs.Z_+obs.H_).Pt(), wgt);
    MJJ_->Fill((obs.VBF_[0]+obs.VBF_[1]).M(), wgt);
    DETAJJ_->Fill(TMath::Abs(obs.VBF_[0].Eta()-obs.VBF_[1].Eta()), wgt);
}

void LHEAnalysis::Histogram::writeHistogram(TFile* output) {
    output->cd();
    MVVH_->Write();
    LT_->Write();
    ST_->Write();
    MJJ_->Write();
    DETAJJ_->Write();
}
