#ifndef __TIMBER_PNETWQQ_WEIGHT__
#define __TIMBER_PNETWQQ_WEIGHT__

#include "TFile.h"
#include "TH2F.h"
#include <ROOT/RVec.hxx>
#include <string>
#include <set>
#include <algorithm>
#include <math.h>
#include "common.h"

/**
 * C++ class to handle the application of ParticleNet MD_Wqq scale factors.
 */
class PNetWqq_weight {
    private:
	std::string _effmapname;	// name of the efficiency map ROOT file
	std::string _year;		// "16", "16APV", "17", "18"
	TFile * _effroot;		// pointer to the efficiency map file
	TH2F * _effmap;			// pointer to the efficiency map histo
	float _wp;			// Wqq working point to distinguish pass/fail (e.g. 0.8)
        int _target_flavor;             // Flavor for SF application (0: top, 1: W, 2: bq, 3: unmerged, -1 "other" (variable))
        std::set<int> _other;           // Jet flavors which constitute the "other" category for SF application

        // https://indico.cern.ch/event/1152827/contributions/4840404/attachments/2428856/4162159/ParticleNet_SFs_ULNanoV9_JMAR_25April2022_PK.pdf
    	// SF[pt][var]

        //Top SF, pt cats: [300,450), [450,600), [600,+inf)
	float SF2016APV_t[3][3] = {{0.992,1.055,0.929},{0.972,1.029,0.915},{1.002,1.103,0.902}};
	float SF2016_t[3][3] = {{0.839,0.907,0.773},{1.001,1.062,0.941},{1.048,1.144,0.956}};
	float SF2017_t[3][3] = {{1.100,1.160,1.041},{1.027,1.070,0.985},{1.154,1.221,1.087}};
	float SF2018_t[3][3] = {{1.055,1.094,1.017},{1.015,1.046,0.984},{0.967,1.025,0.911}};

        //W SF, pt cats: [300,450), [450,600), [600,+inf)
	float SF2016APV_w[3][3] = {{1.004,1.091,0.928},{1.247,1.425,1.103},{0.723,0.952,0.546}};
	float SF2016_w[3][3] = {{1.035,1.160,0.933},{0.869,1.003,0.758},{1.301,1.696,1.014}};
	float SF2017_w[3][3] = {{1.019,1.101,0.947},{0.800,0.847,0.756},{0.877,1.018,0.750}};
	float SF2018_w[3][3] = {{0.784,0.829,0.743},{0.783,0.826,0.743},{0.755,0.855,0.670}};

        //"Other" SF, pt cats: [300,450), [450,600), [600,+inf)
	float SF2016APV_o[3][3] = {{1.025,1.060,0.990},{1.024,1.070,0.978},{1.242,1.333,1.152}};
	float SF2016_o[3][3] = {{0.988,1.026,0.950},{0.968,1.015,0.921},{1.008,1.104,0.915}};
	float SF2017_o[3][3] = {{1.016,1.054,0.979},{1.191,1.228,1.155},{1.192,1.260,1.127}};
	float SF2018_o[3][3] = {{1.080,1.111,1.049},{1.125,1.154,1.096},{1.270,1.329,1.212}};


	float GetMCEfficiency(float pt, float eta, int flavor);	// get the efficiency from efficiency map based on jet's pT and eta
	float GetScaleFactor(float pt, int flavor, int variation);	// get the scale factor based on pT bin and internal year + SF variation (0=nominal, 1=up, 2=down)

    public:
	PNetWqq_weight(std::string year, std::string effmapname, float wp, int target_flavor, std::set<int> other = {2,3}); // pass in year and path of eff map
	~PNetWqq_weight();

	// Return a vector of weights for each event: {nom,up,down}
	RVec<float> eval(float Wqq_pt, float Wqq_eta, float Wqq_PNetWqqScore, int Wqq_jetFlavor);	

};

#endif

