#include "../include/MuonID_weight.h"

// Constructor - specify file path and working point (ID) being used
MuonID_weight::MuonID_weight(std::string correctionname, std::string filepath):_correctionname(correctionname),_filepath(filepath) {}

// Destructor
MuonID_weight::~MuonID_weight() {}

RVec<float> MuonID_weight::eval(float Lepton_eta, float Lepton_pt, int LeptonType) {
    // prepares vector of weights to return
    float sf;
    float sfup;
    float sfdown;

    if (LeptonType == 0) { // Electron event, no sf
        sf = 1;
        sfup = 1;
        sfdown = 1;
    }

    // grab proper sf
    else {
        auto cs = correction::CorrectionSet::from_file(_filepath);        
        auto map = cs->at(_correctionname);
        sf = map->evaluate({std::abs(Lepton_eta), Lepton_pt, "nominal"});
        float stat = map->evaluate({std::abs(Lepton_eta), Lepton_pt, "stat"});
        float syst = map->evaluate({std::abs(Lepton_eta), Lepton_pt, "syst"});
        sfup = sf + stat + syst; 
        sfdown = sf - stat - syst;
    }

    return {sf,sfup,sfdown};
}


