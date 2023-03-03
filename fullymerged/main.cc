#include "main.h"


// All the variables saved to preselected root files
int main(int argc, char** argv)
{
    gconf.nanoAOD_ver = 9;
    // CLI
    HEPCLI cli = HEPCLI(argc, argv);

    // Initialize Looper
    Looper looper = Looper(cli);

    // Initialize Arbol 
    Arbol arbol = Arbol(cli);
    TWOFATJETSCUT::initializeArbol(arbol);

    // Initialize Cutflow
    Cutflow cutflow = Cutflow(cli.output_name + "_Cutflow");
    TWOFATJETSCUT::initializeCutflow(cutflow);

    // Write cutflow
    // Save events and cross sections information
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
            return TWOFATJETSCUT::eq2ElectronsPtGt30(nt, arbol, cutflow);
        }
    );
    cutflow.insert(lhe_vars, lep_sel, Right);

    // Fatjet Selection
    Cut* ak8_sel = new LambdaCut(
        "geq2Fatjets",
        [&]()
        {
            return TWOFATJETSCUT::geq2Fatjets(nt, arbol, cutflow);
        }
    );
    cutflow.insert(lep_sel, ak8_sel, Right);

    // Hbb score selection
    Cut* hbb_sel = new LambdaCut(
        "fatJetsScores",
        [&]()
        {
            return TWOFATJETSCUT::fatJetsScores(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak8_sel, hbb_sel, Right);

    // ak4 jets Selection
    Cut* ak4_sel = new LambdaCut(
        "geq2JetsPtGt30",
        [&]()
        {
            return TWOFATJETSCUT::geq2JetsPtGt30(nt, arbol, cutflow);
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
            return TWOFATJETSCUT::mjjGt500(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak4_sel, mjj_sel, Right);

    // Selection on VBF jet delta Eta
    Cut* deta_sel = new LambdaCut(
        "deltaEtaGt3",
        [&]()
        {
            return TWOFATJETSCUT::deltaEtaGt3(nt, arbol, cutflow);
        }
    );
    cutflow.insert(mjj_sel, deta_sel, Right);

    // Selection on ST
    Cut* ST_sel = new LambdaCut(
        "stgeq950",
        [&]()
        {
            return TWOFATJETSCUT::stgeq950(nt, arbol, cutflow);
        }
    );
    cutflow.insert(deta_sel, ST_sel, Right); 


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