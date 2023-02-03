#include "main.h"


/* 
Per event weight
*/
// All the variables saved to preselected root files
void initializeArbol(Arbol& arbol_) {
    arbol_.newBranch<int>("event", -999);
    arbol_.newBranch<double>("xsec_sf", -999);
    // LHE variables
    arbol_.newBranch<double>("lhe_mvvh", -999);
    arbol_.newBranch<double>("lhe_lt", -999);
    arbol_.newBranch<double>("lhe_st", -999);
    arbol_.newBranch<double>("lhe_mjj", -999);
    arbol_.newBranch<double>("lhe_detajj", -999);
    // Reconstructed variables
    arbol_.newBranch<LorentzVector>("vbfjet1_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("vbfjet2_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("wjet1_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("wjet2_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("lep1_LV", {-999,-999,-999,-999});
    arbol_.newBranch<LorentzVector>("lep2_LV", {-999,-999,-999,-999});
    arbol_.newBranch<int>("lep1_ID", -999);
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
    cutflow_.globals.newVar<float>("met_pt");
    cutflow_.globals.newVar<float>("ST");
    cutflow_.globals.newVar<float>("LT");
    cutflow_.globals.newVar<bool>("medium_b_veto");
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
                bool LHE_passed=cutflow.run("SelectLHEVariables");
                if (LHE_passed) { arbol.fill(); }

                /*
                bool dummycut1_passed = cutflow.run("DummyCut1");
                if (dummycut1_passed) { arbol.fill(); }
                */
                
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