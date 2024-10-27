#include "../include/MuonID_weight.h"

// Constructor - specify file path and working point (ID) being used
MuonID_weight::MuonID_weight(std::string correctionname, std::string filepath):_correctionname(correctionname),_filepath(filepath) {
    _eta_bounds[0] = 0;
    _eta_bounds[1] = 2.4;
    _pt_bounds[0] = 50; 
    if (_filepath.find("IDISO") != std::string::npos) {
        _pt_bounds[1] = 1000;
    }
    else if (_filepath.find("RECO") != std::string::npos){
        _pt_bounds[1] = 3500;
    }
}

// Destructor
MuonID_weight::~MuonID_weight() {}

bool MuonID_weight::inRange(float bounds[2], float x) {
    float low = bounds[0];
    float high = bounds[1];
    return (low <= x && x <= high);
}

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

    else if (!inRange(_eta_bounds, Lepton_eta) || !inRange(_pt_bounds, Lepton_pt)) {
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


