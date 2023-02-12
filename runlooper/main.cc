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
    arbol_.newBranch<int>("n_ak4", -999); // # of ak4 jets
    arbol_.newBranch<int>("n_ak8", -999);
    arbol_.newBranch<double>("hbb_score", -999);
    arbol_.newBranch<LorentzVector>("hbb_LV", {-999,-999,-999,-999});
    arbol_.newBranch<double>("xbb_score", -999);
    arbol_.newBranch<double>("hbb_pnetmass", -999);
    arbol_.newBranch<bool>("b_veto", false);
}

void initializeCutflow(Cutflow& cutflow_) {
    cutflow_.globals.newVar<LorentzVector>("ld_lep_p4"); 
    cutflow_.globals.newVar<LorentzVector>("sd_lep_p4"); 
    cutflow_.globals.newVar<int>("ldid"); // leading lepton pdgid
    cutflow_.globals.newVar<int>("sdid"); // subleading lepton pdgid
    cutflow_.globals.newVar<LorentzVector>("hbbjet_p4");
    cutflow_.globals.newVar<LorentzVector>("xbbjet_p4");
    // # of ak8 jets that pass the requirements: softdropmass larger than 40, pt greater than 250, |eta|<2.5, fatjetID>0
    cutflow_.globals.newVar<int>("n_ak8");
    cutflow_.globals.newVar<int>("n_ak4");
    cutflow_.globals.newVar<LorentzVector>("ld_jet_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_jet_p4");
    cutflow_.globals.newVar<LorentzVector>("ld_vbf_p4");
    cutflow_.globals.newVar<LorentzVector>("sd_vbf_p4");
    cutflow_.globals.newVar<std::vector<int>>("jetidx");
    cutflow_.globals.newVar<double>("mjj");
    cutflow_.globals.newVar<double>("detajj");
    cutflow_.globals.newVar<float>("met_pt");
    cutflow_.globals.newVar<float>("ST");
    cutflow_.globals.newVar<float>("LT");
    cutflow_.globals.newVar<bool>("medium_b_veto");
    cutflow_.globals.newVar<double>("hbb"); // highest hbb score (hbb score of hbb jet)
    cutflow_.globals.newVar<double>("xbb"); // highest xbb score (not necessarily that of hbb jet)
    cutflow_.globals.newVar<double>("pn_mass"); // particle net mass of hbb jet
    cutflow_.globals.newVar<bool>("samejet"); // Check if highest xbb and highest hbb come from the same jet
}

int main(int argc, char** argv)
{
    gconf.nanoAOD_ver = 9;
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
    
    // =========================================================================================
    // OBJECT SELECTION
    // Electron/Muon selection
    Cut* lep_sel = new LambdaCut(
        "eq2ElectronsPtGt30",
        [&]()
        {
            return eq2ElectronsPtGt30(nt, arbol, cutflow);
        }
    );
    cutflow.insert(lhe_vars, lep_sel, Right);

    // Fatjet Selection
    Cut* ak8_sel = new LambdaCut(
        "geq1FatjetsMassGt40",
        [&]()
        {
            return geq1FatjetsMassGt40(nt, arbol, cutflow);
        }
    );
    cutflow.insert(lep_sel, ak8_sel, Right);

    // Hbb score selection
    Cut* hbb_sel = new LambdaCut(
        "fatJetHbbScore",
        [&]()
        {
            return fatJetHbbScore(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak8_sel, hbb_sel, Right);

    // ak4 jets Selection
    Cut* ak4_sel = new LambdaCut(
        "geq4JetsPtGt30",
        [&]()
        {
            return geq4JetsPtGt30(nt, arbol, cutflow);
        }
    );
    cutflow.insert(hbb_sel, ak4_sel, Right);

    // =========================================================================================
    // PRESELECTION
    // Selection on VBF jet system mass
    Cut* mjj_sel = new LambdaCut(
        "mjjGt500",
        [&]()
        {
            return mjjGt500(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak4_sel, mjj_sel, Right);

    // Selection on VBF jet delta Eta
    Cut* deta_sel = new LambdaCut(
        "deltaEtaGt3",
        [&]()
        {
            return deltaEtaGt3(nt, arbol, cutflow);
        }
    );
    cutflow.insert(mjj_sel, deta_sel, Right);

    // Selection on ST
    Cut* ST_sel = new LambdaCut(
        "stgeq950",
        [&]()
        {
            return stgeq950(nt, arbol, cutflow);
        }
    );
    cutflow.insert(deta_sel, ST_sel, Right);

    // Cut* dummycut1 = new LambdaCut("DummyCut1", [&]() { return true; });
    // cutflow.insert(lep_sel, dummycut1, Right);

    // Intialize progress bar
    tqdm bar;

    // Run looper
    looper.run(
        [&](TTree* ttree)
        {
            nt.Init(ttree);
            TString file_name = cli.input_tchain->GetCurrentFile()->GetName();
            gconf.isAPV = (
                file_name.Contains("HIPM_UL2016")
                || file_name.Contains("NanoAODAPV")
                || file_name.Contains("UL16APV")
            );
            gconf.GetConfigs(nt.year());
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
                bool all_passed=cutflow.run("stgeq950");
                if (all_passed) { arbol.fill(); }

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