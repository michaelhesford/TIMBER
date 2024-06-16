#ifndef __TIMBER_MUONID_WEIGHT__
#define __TIMBER_MUONID_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <algorithm>
#include <math.h>
#include "common.h"
#include <fstream>
#include "/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/py3-correctionlib/2.2.2-cb60c4327c0522c2e7ee31963c98a46f/lib/python3.9/site-packages/correctionlib/include/correction.h"

//using json = nlohmann::json;

/*
-C++ class to handle application of scale factors for muon ID
-For muon scale factors we don't have the luxury of using correctionlib to easily access appropriate values from the centrally-produced json files, so I am creating a separate class for muons with some added bulk (otherwise, things would ideally be consolidated)
*/

class MuonID_weight {
    private:
        std::string _correctionname;  // "value to correct (mediumID, tightID etc...) - see json file for exact names
        std::string _filepath;        // "path to file with SF's
        int _leptonType;              // 0: electron, 1: muon (must only apply SF to events where lepton is a muon
        int _ptBin;                   // for indexing data by pT
        int _etaBin;                  // for indexing data by abs(eta)
        int _wpBin;                   // to access proper data in json (are you sure there's no tool for this?)
    
    public:
        MuonID_weight(std::string correctionname, std::string filepath);
        ~MuonID_weight();
 
        // Return a vector of weights for each event: {nom,up,down}
        RVec<float> eval(float Lepton_eta, float Lepton_pt, int LeptonType);
};

#endif        

 
