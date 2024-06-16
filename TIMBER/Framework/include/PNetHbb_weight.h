#ifndef __TIMBER_PNETHBB_WEIGHT__
#define __TIMBER_PNETHBB_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <set>
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
	float _wp;			// Hbb working point delineating Faill/Pass (e.g 0.98)
        int _target_flavor;             // Flavor of jet to implement SF for (0: top, 2: bq, 4: Higgs)
        std::set<int> _other;           // Flavors of jets to consider for "other" classification/SF application

    	// SF[pt][_var]
    	// WP (tight) [0.98, 1.0]
    
        // HIGGS JET SF'S: pT cats are [400, 600), [600, 800), [800, +inf) across all years
        std::vector<std::vector<float>> SF2016APV_h   = {{1.163,1.437,0.949},{1.206,1.529,0.924},{1.491,2.083,0.909}};
        std::vector<std::vector<float>> SF2016_h   = {{1.012,1.180,0.899},{1.247,1.509,0.999},{1.188,1.424,0.960}};
        std::vector<std::vector<float>> SF2017_h    = {{0.946,1.075,0.830},{1.027,1.158,0.880},{0.900,1.026,0.752}};
        std::vector<std::vector<float>> SF2018_h    = {{1.020,1.146,0.894},{1.013,1.110,0.912},{1.082,1.240,0.961}};

        // TOP JET SF'S: pT cats are [300, 450), [450, 600), [600, +inf)
        std::vector<std::vector<float>> SF2016APV_t = {{0.807,0.942,0.679},{0.763,0.871,0.657},{0.664,0.816,0.523}};
        std::vector<std::vector<float>> SF2016_t = {{0.920,1.075,0.775},{0.930,1.056,0.810},{1.045,1.242,0.864}};
        std::vector<std::vector<float>> SF2017_t = {{0.691,0.807,0.580},{1.055,1.144,0.969},{1.143,1.278,1.014}};
        std::vector<std::vector<float>> SF2018_t = {{0.738,0.819,0.660},{0.977,1.040,0.917},{1.014,1.123,0.910}};

        // BQ JET SF'S: pT cats are [300, 450), [450, 600), [600, +inf)
	std::vector<std::vector<float>> SF2016APV_bq = {{1.284,1.428,1.154},{1.261,1.442,1.097},{0.879,1.309,0.531}};
	std::vector<std::vector<float>> SF2016_bq = {{1.154,1.302,1.022},{0.962,1.152,0.797},{1.116,1.634,0.708}};
	std::vector<std::vector<float>> SF2017_bq = {{1.275,1.395,1.165},{1.333,1.475,1.202},{1.440,1.773,1.170}};
	std::vector<std::vector<float>> SF2018_bq = {{1.558,1.682,1.442},{1.349,1.476,1.233},{1.555,1.933,1.242}};


	float GetMCEfficiency(float pt, float eta, int flavor);	    // get the efficiency from efficiency map based on jet's pT and eta (for given flavor)
	std::vector<float> GetScaleFactors(float pt, int flavor);    // get the scale factor based on pT bin and internal year + tagger category variables (for given flavor)

    public:
	PNetHbb_weight(std::string year, std::string effmapname, float wp, int target_flavor, std::set<int> other = {1,3}); // pass in year and path of eff map
	~PNetHbb_weight();

	// Return a vector of weights for each event: {nom,up,down}
	RVec<float> eval(float Higgs_pt, float Higgs_eta, float Higgs_PNetHbbScore, int Higgs_jetFlavor);	

};

#endif

