#ifndef __TIMBER_ELECTRONID_WEIGHT__
#define __TIMBER_ELECTRONID_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <algorithm>
#include <math.h>
#include "common.h"
#include "/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/py3-correctionlib/2.2.2-cb60c4327c0522c2e7ee31963c98a46f/lib/python3.9/site-packages/correctionlib/include/correction.h"

/*
C++ class to handle the application of scale factors for the electron ID.
*/

class ElectronID_weight {
    private:
        std::string _year;           // "16APV", "16", "17", "18"
        std::string _fullyear;       // "2016preVFP", "2016postVFP", "2017", "2018" - formatting for correctionlib
        std::string _wp;             // Cut-based: "Loose", "Medium", "Tight" - MVA: "wp80noiso", "wp80iso", "wp90noiso", "wp90iso"
        std::vector<std::pair<std::string,std::string>> _variations = {{"sf","ElectronIDWeight"},{"sfup","ElectronIDWeight_Up"},{"sfdown","ElectronIDWeight_Down"}}; // {variation,branchname}
        std::string _variation;      // "sf", "sfup", "sfdown"        
        std::string _filepath;       // path to file with SF's
        float GetScaleFactor(std::string variation, float eta, float pt, int type); // retreive SF using correctionlib
        int _leptonType; // 0: electron, 1: muon (must only apply SF to events where lepton is an electron)

    public:
        ElectronID_weight(std::string year, std::string wp, std::string filepath); 
        ~ElectronID_weight();

        // Return a vector of weights for each event: {nom,up,down}
        RVec<float> eval(float Lepton_eta, float Lepton_pt, int LeptonType);
};

#endif
