// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"
#include "ElectronSelections.h"
#include "MuonSelections.h"


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
        double deltaEta=-999.999;
        for (unsigned int jetsi = 0; jetsi < jets.size(); jetsi++) {
            for (unsigned int jetsj = jetsi+1; jetsj < jets.size(); jetsj++) {
                double this_deltaE=fabs(jets[jetsi].Eta()-jets[jetsj].Eta());
                if (this_deltaE > deltaEta) {
                    deltaEta=this_deltaE;
                    VBFindx1=jetsi;
                    VBFindx2=jetsj;
                }
            }
        }
        
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