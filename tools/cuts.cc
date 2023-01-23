#include "cuts.h"

bool SelectLHEVariables::evaluate() {
    unsigned int particleNo=nt.LHEPart_status().size();
    std::vector<LorentzVector> jets;
    LorentzVector W,H,Z;
    Int_t VBFindx1, VBFindx2;
    Float_t eta1=-1000;
    Float_t eta2=1000;
    
    // Selection based on pdgid
    for (unsigned int i=0; i<particleNo; i++) {
        // Check if it's outgoing particles
        if (nt.LHEPart_status()[i]==1) {
            Int_t pdgid=nt.LHEPart_pdgId()[i];
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
    
    // Select the VBF jets with the largest etas in both directions
    for (unsigned int j=0; j<jets.size(); j++) {
        if (jets[j].Eta()>0 && jets[j].Eta()>eta1) {
            VBFindx1=j;
        }
        else if (jets[j].Eta()<0 && jets[j].Eta()<eta2) {
            VBFindx2=j;
        }
    }
    
    // Set variables to arbol tree
    arbol.setLeaf<double>("mvvh", (W+H+Z).M());
    arbol.setLeaf<double>("lt", Z.Pt());
    arbol.setLeaf<double>("st", (W+H+Z).Pt());
    arbol.setLeaf<double>("mjj", (jets.at(VBFindx1)+jets.at(VBFindx2)).M());
    arbol.setLeaf<double>("detajj", (TMath::Abs(jets.at(VBFindx1).Eta()-jets.at(VBFindx2).Eta())));

    /* Cut logic */
    return true;
}

