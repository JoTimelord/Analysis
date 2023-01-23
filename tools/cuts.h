// RAPIDO
#include "cutflow.h"
#include "arbol.h"
#include "utilities.h"
// NanoCORE
#include "Nano.h"
#include "Config.h"

bool geq2ElectronsPtGt40(Nano& nt, Arbol& arbol, Cutflow& cutflow)
{
    int n_electrons = 0;
    for (unsigned int elec_i = 0; elec_i < nt.nElectron(); ++elec_i)
    {
        if (nt.Electron_pt().at(elec_i) > 40)
        {
            cutflow.globals.setVal<LorentzVector>("elec_p4", nt.Electron_p4().at(elec_i));
            n_electrons++;
        }
    }
    return (n_electrons >= 2);
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
    bool evaluate()
    {
        return true;
    };

    double weight()
    {
        /* Event weight (applied only if Cut::evaluate() returns true) */
        return 1.;
    };

    virtual ~SelectLHEVariables() {}
};