#ifndef __TIMBER_MUONID_WEIGHT__
#define __TIMBER_MUONID_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <algorithm>
#include <math.h>
#include <cmath>
#include "common.h"
#include <fstream>
#include "/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/py3-correctionlib/2.1.0-18f6d1dfcf4c205c38d543f9ef2a014b/lib/python3.9/site-packages/correctionlib/include/correction.h"

/*
-C++ class to handle application of scale factors for muon ID
-For muon scale factors we don't have the luxury of using correctionlib to easily access appropriate values from the centrally-produced json files, so I am creating a separate class for muons with some added bulk (otherwise, things would ideally be consolidated)
*/

class MuonID_weight {
    private:
        std::string _correctionname;  // "value to correct (mediumID, tightID etc...) - see json file for exact names
        std::string _filepath;        // "path to file with SF's
        int _leptonType;              // 0: electron, 1: muon (must only apply SF to events where lepton is a muon
	float _pt_bounds[2];
        float _eta_bounds[2];
    
    public:
        MuonID_weight(std::string correctionname, std::string filepath);
        ~MuonID_weight();

        bool inRange(float bounds[2], float x);
 
        // Return a vector of weights for each event: {nom,up,down}
        RVec<float> eval(float Lepton_eta, float Lepton_pt, int LeptonType);
};

#endif        

 
