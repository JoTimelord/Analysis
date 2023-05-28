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
            arbol.setLeaf<double>("scale_fac", cli.scale_factor);
            arbol.setLeaf<double>("gen_wgt", nt.genWeight());
            arbol.setLeaf<double>("xsec_sf", cli.scale_factor*nt.genWeight());
            if (cli.is_signal && nt.nLHEReweightingWeight() > 0)
            {
                arbol.setLeaf<double>("reweight_c2v_eq_3", nt.LHEReweightingWeight().at(31));
                arbol.setLeaf<double>("reweight_c2v_eq_4", nt.LHEReweightingWeight().at(35));
                arbol.setLeaf<double>("gen_wgt",1.0);
            }
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
        "eq2OSLeptonsPtGt30",
        [&]()
        {
            return eq2OSLeptonsPtGt30(nt, arbol, cutflow);
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

    // ak4 jets Selection
    Cut* ak4_sel = new LambdaCut(
        "geq2JetsPtGt30",
        [&]()
        {
            return TWOFATJETSCUT::geq2JetsPtGt30(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak8_sel, ak4_sel, Right);

    // =========================================================================================
    // PRESELECTION
    // Hbb score selection
    /*
    Cut* hbb_sel = new LambdaCut(
        "hbbScore",
        [&]()
        {
            return TWOFATJETSCUT::hbbScore(nt, arbol, cutflow);
        }
    );
    cutflow.insert(ak4_sel, hbb_sel, Right);

    // Vxqq score selection
    Cut* vqq_sel = new LambdaCut(
        "xvqqScore",
        [&]()
        {
            return TWOFATJETSCUT::xvqqScore(nt, arbol, cutflow);
        }
    );
    cutflow.insert(hbb_sel, vqq_sel, Right);

    // Selection on VBF jet system mass
    Cut* mjj_sel = new LambdaCut(
        "mjjGt500",
        [&]()
        {
            return mjjGt500(nt, arbol, cutflow);
        }
    );
    cutflow.insert(vqq_sel, mjj_sel, Right);

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
            return TWOFATJETSCUT::stgeq950(nt, arbol, cutflow);
        }
    );
    cutflow.insert(deta_sel, ST_sel, Right); 
    */

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
                bool all_passed=cutflow.run("geq2JetsPtGt30");
                if (all_passed) { arbol.fill(); }
                
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