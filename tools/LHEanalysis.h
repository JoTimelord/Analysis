#include "Nano.h"
#include "looper.h"
#include "calc.h"

namespace LHEAnalysis{
    // Define a class called Observable that extract LHE variables.
    class Observable {
    private:
        LV W_;
        LV H_;
        LV Z_;
        std::vector<LV> VBF_;
    public:
        // Used for each individual event to select the LHE level observables we care about.
        void selectObservable();
        // Destructor
        // ~Observable();
        // Allow Histogram to have access to LHE observables; 
        friend class Histogram;
    };

    /* Define a class called Histogram which is a friend class of Observable.
    This class handles the declare, define and initialize the type of histograms we care about in LHE-level Analysis. */ 
    class Histogram {
    private: 
        // Create pointer members to THistogram
        TH1F* MVVH_; 
        TH1F* LT_;
        TH1F* ST_;
        TH1F* MJJ_;
        TH1F* DETAJJ_;
    public: 
        // Constructor 
        Histogram();
        // Fill the THist objects
        void fillHistogram(LHEAnalysis::Observable& obs, Int_t wgt=1);
        // Write THist to file
        void writeHistogram(TFile* output);
    };
}