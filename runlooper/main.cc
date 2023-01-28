#include "main.h"

void initializeArbol(Arbol& arbol_) {
    arbol_.newBranch<int>("event", -999);
    arbol_.newBranch<double>("xsec_sf", -999);
    arbol_.newBranch<double>("mvvh", -999);
    arbol_.newBranch<double>("lt", -999);
    arbol_.newBranch<double>("st", -999);
    arbol_.newBranch<double>("mjj", -999);
    arbol_.newBranch<double>("detajj", -999);
}
void initializeCutflow(Cutflow& cutflow_) {
    cutflow_.globals.newVar<LorentzVector>("ld_lep_p4"); 
    cutflow_.globals.newVar<LorentzVector>("sd_lep_p4"); 
    cutflow_.globals.newVar<LorentzVector>("fatjet_p4");
    cutflow_.globals.newVar<LorentzVector>("ld_jet_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_jet_p4");
    cutflow_.globals.newVar<LorentzVector>("ld_vbf_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_vbf_p4");
    cutflow_.globals.newVar<std::vector<int>>("jetidx");
}

int main(int argc, char** argv)
{
    // CLI
    HEPCLI cli = HEPCLI(argc, argv);

    // Initialize Looper
    Looper looper = Looper(cli);

    // Initialize Arbol 
    Arbol arbol = Arbol(cli);
    initializeArbol(arbol);

    // Initialize Cutflow
    Cutflow cutflow = Cutflow(cli.output_name + "_Cutflow");
    initializeCutflow(cutflow);

    Cut* base = new LambdaCut(
        "BaseCut",
        [&]()
        {
            /* Cut logic */
            arbol.setLeaf<int>("event", nt.event());
            arbol.setLeaf<double>("xsec_sf", cli.scale_factor*nt.genWeight());
            return true;
        },
        [&]()
        {
            /* Event weight (applied only if event passes selection above) */
            return arbol.getLeaf<double>("xsec_sf");
        }
    );
    cutflow.setRoot(base);

    // Save LHE-level information
    Cut* lhe_vars = new SelectLHEVariables("SelectLHEVariables", nt, arbol, cutflow);
    cutflow.insert(base, lhe_vars, Right);

    Cut* lep_sel = new LambdaCut(
        "eq2ElectronsPtGt30",
        [&]()
        {
            return eq2ElectronsPtGt30(nt, arbol, cutflow);
        }
    );
    cutflow.insert(lhe_vars, lep_sel, Right);

    Cut* dummycut1 = new LambdaCut("DummyCut1", [&]() { return true; });
    cutflow.insert(lep_sel, dummycut1, Right);

    // Intialize progress bar
    tqdm bar;

    // Run looper
    looper.run(
        [&](TTree* ttree)
        {
            nt.Init(ttree);
        },
        [&](int entry)
        {
            if (cli.debug && looper.n_events_processed == 10000) { looper.stop(); }
            else
            {
                // Reset branches and globals
                arbol.resetBranches();
                cutflow.globals.resetVars();
                // Run cutflow
                nt.GetEntry(entry);
                bool LHE_passed=cutflow.run("lhe_vars");
                if (LHE_passed) { arbol.fill(); }

                bool dummycut1_passed = cutflow.run("DummyCut1");
                if (dummycut1_passed) { arbol.fill(); }
                // Update progress bar
                bar.progress(looper.n_events_processed, looper.n_events_total);
            }
        }
    );

    // Wrap up
    if (!cli.is_data)
    {
        cutflow.print();
        cutflow.write(cli.output_dir);
    }
    arbol.write();
    return 0;
}