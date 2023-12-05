#ifndef __TIMBER_PNETHBB_WEIGHT__
#define __TIMBER_PNETHBB_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <algorithm>
#include <math.h>
#include "common.h"

/**
 * C++ class to handle the application of ParticleNet MD_Hbb scale factors.
 */
class PNetHbb_weight {
    private:
	std::string _effmapname;	// name of the efficiency map ROOT file
	std::string _year;		// "16", "16APV", "17", "18"
	TFile * _effroot;		// pointer to the efficiency map file
	TH2F * _effmap;			// pointer to the efficiency map histo
	std::vector<std::pair<int, std::string>> _variations = {{0,"PNetHbbWeight"},{1,"PNetHbbWeight_Up"},{-1,"PNetHbbWeight_Down"}}; // nom, up, down
	int _variation;			// 0: nominal, 1: up, 2: down
	float _wp;			// Hbb working point delineating Faill/Pass (e.g 0.98)

    	// SF[_var][pt]
    	// pT cats are [400, 600), [600, 800), [800, +inf) across all years
    	// HP (tight) [0.98, 1.0]
    	float SF2016APV_T[3][3] = {{1.163,1.206,1.491},{1.437,1.529,2.083},{0.949,0.924,0.909}};
    	float SF2016_T[3][3]    = {{1.012,1.247,1.188},{1.180,1.509,1.424},{0.899,0.999,0.960}};
    	float SF2017_T[3][3]    = {{0.946,1.027,0.900},{1.075,1.158,1.026},{0.830,0.880,0.752}};
    	float SF2018_T[3][3]    = {{1.020,1.013,1.082},{1.146,1.110,1.240},{0.894,0.912,0.961}};
    	// MP (loose) [0.8, 0.98]
    	float SF2016APV_L[3][3] = {{1.102,1.103,0.645},{1.321,1.355,1.914},{0.918,0.871,0.955}};
    	float SF2016_L[3][3]    = {{1.032,1.173,1.145},{1.134,1.382,1.332},{0.932,0.970,0.954}};
    	float SF2017_L[3][3]    = {{0.973,1.006,1.059},{1.026,1.064,1.132},{0.904,0.931,0.982}};
    	float SF2018_L[3][3]    = {{0.904,0.921,1.087},{0.966,0.969,1.165},{0.824,0.841,0.975}};

	float GetMCEfficiency(float pt, float eta);	// get the efficiency from efficiency map based on jet's pT and eta
	float GetScaleFactor(float pt, float score);	// get the scale factor based on pT bin and internal year + tagger category variables

    public:
	PNetHbb_weight(std::string year, std::string effmapname, float wp); // pass in year and path of eff map
	~PNetHbb_weight();

	// Return a vector of weights for each event: {nom,up,down}
	RVec<float> eval(RVec<float> FatJet_pT, RVec<float> FatJet_eta, // for eff map
			 RVec<float> FatJet_PNetHbbScore);		// Hbb scores for all FatJets

};

#endif

