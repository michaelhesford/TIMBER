#include "../include/ElectronID_weight.h"

// Constructor - set year and working point/ID being used
ElectronID_weight::ElectronID_weight(std::string year, std::string wp, std::string filepath) { 
    _year = year;
    if (_year == "16APV") {_fullyear = "2016preVFP";}
    else if (_year == "16") {_fullyear = "2016postVFP";}
    else {_fullyear = "20"+_year;}
    _wp = wp;
    //_filepath = "/uscms/home/mhesford/nobackup/XHYbbWW/CMSSW_12_3_0/src/semileptonic/corrections/electron_"+_year+".json";
    _filepath = filepath;
}

// Destructor
ElectronID_weight::~ElectronID_weight() {}

float ElectronID_weight::GetScaleFactor(std::string variation, float eta, float pt, int type) {
    auto cs = correction::CorrectionSet::from_file(_filepath);
    auto map = cs->at("UL-Electron-ID-SF");
    float sf;
    if (type == 1) {sf = 1;} // muon events, no scale factor
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

 
