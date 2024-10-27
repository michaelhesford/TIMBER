#include "../include/ElectronID_weight.h"

// Constructor - set year and working point/ID being used
ElectronID_weight::ElectronID_weight(std::string year, std::string wp, std::string filepath) { 
    _year = year;
    if (_year == "16APV") {_fullyear = "2016preVFP";}
    else if (_year == "16") {_fullyear = "2016postVFP";}
    else {_fullyear = "20"+_year;}
    _wp = wp;
    _filepath = filepath;
    _eta_bounds[0] = -INFINITY;
    _eta_bounds[1] = INFINITY;
    if (_wp == "RecoBelow20") {
        _pt_bounds[0] = 10; 
	_pt_bounds[1] = 20;
    }
    else if (_wp == "RecoAbove20") {
        _pt_bounds[0] = 20;
	_pt_bounds[1] = INFINITY;
    }
    else {
        _pt_bounds[0] = 10;
	_pt_bounds[1] = INFINITY;
    }
}

// Destructor
ElectronID_weight::~ElectronID_weight() {}

bool ElectronID_weight::inRange(float bounds[2], float x) {         
    float low = bounds[0];
    float high = bounds[1];
    return (low <= x && x <= high);          
}    

float ElectronID_weight::GetScaleFactor(std::string variation, float eta, float pt, int type) {
    auto cs = correction::CorrectionSet::from_file(_filepath);
    auto map = cs->at("UL-Electron-ID-SF");
    float sf;
    if (type == 1) {sf = 1;} // muon events, no scale factor
    else if (!inRange(_eta_bounds, eta) || !inRange(_pt_bounds, pt)) {sf = 1;}
    else { // electron event
        sf = map->evaluate({_fullyear, _variation, _wp, std::abs(eta), pt});
    }
    return sf;
}

RVec<float> ElectronID_weight::eval(float Lepton_eta, float Lepton_pt, int LeptonType) {
    // prepares vector of weights to return
    RVec<float> out(3);
    std::string branchname;
    // Loop over variations ("sf","sfup","sfdown")
    for (size_t i=0; i<_variations.size(); i++) {
        _variation = _variations[i].first;
        branchname = _variations[i].second;
        float SF = GetScaleFactor(_variation, Lepton_eta, Lepton_pt, LeptonType);
        out[i] = SF;
    }
    return out;
}

 
