#include "main.h"

// Read in possible arguments and set input and output files
void readinput(int inputno, char** arguments, std::string& input_path, std::string& output_path, int& n_events) {
    // Grand option setting
    cxxopts::Options options("\n  $ parseLHEs",  "\n         **********************\n         *                    *\n         *       Looper       *\n         *                    *\n         **********************\n");

    // Read the options
    options.add_options()
        ("i,input", "Comma separated input file list OR if just a directory is provided it will glob all in the directory BUT must end with '/' for the path", cxxopts::value<std::string>())
        ("o,output", "Output file name", cxxopts::value<std::string>())
        ("n,nevents", "N events to loop over", cxxopts::value<int>()->default_value("-1"))
        ("h,help"        , "Print help")
        ;

    auto result = options.parse(argc, argv);

    //_______________________________________________________________________________
    // --help
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    } 

    //_______________________________________________________________________________
    // --input 
    if (result.count("input"))
    {
	input_path = result["input"].as<std::string>();
    }
    else 
    {
	std::cout << options.help() << std::endl;
	std::cout << "No input provided" << std::endl;
        exit(1);
    } 
    //_______________________________________________________________________________
    // --nevents
    n_events = result["nevents"].as<int>();
    
    //_______________________________________________________________________________
    // --output
    if (result.count("output"))
    {
	output_path = result["output"].as<std::string>();
    }
    else 
    {
	std::cout << options.help() << std::endl;
	std::cout << "No output provided" << std::endl;
        exit(1);
    } 
}

int main(int argc, char** argv) {
    // Files to loop through
    std::string input_path;
    std::string output_path;
    int nevents; 
    
    // Configuration for nt and NanoCORE
    nt.SetYear(2018);
    gconf.GetConfigs(nt.year());
    readinput(argc, argv, input_path, output_path, nevents);
    // Create output root file
    TFile* ofile=new TFile(output_path.c_str(),"recreate");
    // Create LHE level Histograms Class
    LHEAnalysis::Histogram H;
    // Create TChain for events looping
    TString input_tree="Events";
    TChain* events_chain=RooUtil::FileUtil::createTChain(input_tree, input_path)
    // Create a looper
    RooUtil::Looper<Nano> looper;
    // Initialize the looper
    looper.init(events_tchain, &nt, n_events);
    // Loop through events
    while (looper.nextEvent()) {
        LHEAnalysis::Observable obs;
        obs.selectObservable();
        H.fillHistogram(obs);
    }
    H.writeHistogram(ofile);
    ofile->Close();
    return 0;
}