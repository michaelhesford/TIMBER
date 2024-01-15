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
        /*
        // get wp bin - this is kinda stupid ngl, why would they format the json like this??
        if (_wp == "NUM_HighPtID_DEN_GlobalMuonProbes") {_wpBin = 0;}
	else if (_wp == "NUM_TrkHighPtID_DEN_GlobalMuonProbes") {_wpBin = 1;}
        else if (_wp == "NUM_probe_LooseRelTkIso_DEN_HighPtProbes") {_wpBin = 2;}
        else if (_wp == "NUM_probe_TightRelTkIso_DEN_HighPtProbes") {_wpBin = 3;}
        else if (_wp == "NUM_probe_LooseRelTkIso_DEN_TrkHighPtProbes") {_wpBin = 4;}
        else if (_wp == "NUM_probe_TightRelTkIso_DEN_TrkHighPtProbes") {_wpBin = 5;}
        else if (_wp == "NUM_TightID_DEN_GlobalMuonProbes") {_wpBin = 6;}
        else if (_wp == "NUM_MediumID_DEN_GlobalMuonProbes") {_wpBin = 7;}
        else if (_wp == "NUM_probe_LooseRelTkIso_DEN_MediumIDProbes") {_wpBin = 8;}
        else if (_wp == "NUM_probe_TightRelTkIso_DEN_MediumIDProbes") {_wpBin = 9;}
        
        // get eta bin
        if (std::abs(Lepton_eta) < 0.9) {_etaBin = 0;}
        else if (std::abs(Lepton_eta) > 0.9 && std::abs(Lepton_eta) < 1.2) {_etaBin = 1;}
        else if (std::abs(Lepton_eta) > 1.2 && std::abs(Lepton_eta) < 2.1) {_etaBin = 2;}
        else if (std::abs(Lepton_eta) > 2.1 && std::abs(Lepton_eta) < 2.4) {_etaBin = 3;}
        
        // get pt bin
        if (Lepton_pt < 60) {_ptBin = 0;} // note: the actual bin is 50-60 here, so I'm cheating
        else if (Lepton_pt > 60 && Lepton_pt < 120) {_ptBin = 1;}
        else if (Lepton_pt > 120 && Lepton_pt < 200) {_ptBin = 2;}
        else {
            if (_etaBin == 3) {_ptBin = 3;} // binning goes from 200-1000 for the highest eta bin
            else {
                if (Lepton_pt > 200 && Lepton_pt < 450) {_ptBin = 3;}
                else if (Lepton_pt > 450) {_ptBin = 4;} // actual bin is 450-1000, so again I'm cheating
            }
        }
        
        // extract data
        std::ifstream file(_filepath);
        json data = json::parse(file);
        
        // grab appropriate scale factor              
        json content = data["corrections"][_wpBin]["data"]["content"][_etaBin]["content"][_ptBin]["content"];
        sf = content[4]["value"].template get<float>();
        float stat = content[1]["value"].template get<float>();
        float syst = content[2]["value"].template get<float>();
        sfup = sf + stat + syst;
        sfdown = sf - stat - syst; 
        */
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


