#include "../include/PNetHbb_weight.h"

// Constructor - set year and load eff map in memory, 
PNetHbb_weight::PNetHbb_weight(std::string year, std::string effmapname, float wp) {
    // Change this logic yourself based on how you made the efficiency map.
    // I am assuming that it is one ROOT file per year, with one histogram.
    _effmapname = effmapname;
    _year = year;
    _effroot = hardware::Open(_effmapname);
    _effmap = (TH2F*)_effroot->Get("HbbJets ratio");
    _wp = wp;
}

// Close files on destruction
PNetHbb_weight::~PNetHbb_weight() {
    _effroot->Close();
}

// Get the MC effieciency of the jets in the given (pT, eta) bin
float PNetHbb_weight::GetMCEfficiency(float pt, float eta) {
    // Efficiency map binned in pT: [60,0,3000], eta: [24,-2.4,2.4]
    int xbin = (int)(pt*60./3000.);
    int ybin = (int)((eta+2.4)*24/4.8);
    return _effmap->GetBinContent(xbin,ybin);
}

// Get the SF for the jet based on its pT and score and (internal) variation
float PNetHbb_weight::GetScaleFactor(float pt, float score) {
    // First check which pT bin we're in
    int ptBin;
    if (pt >= 400 && pt < 600) {ptBin = 0;}
    else if (pt >= 600 && pt < 800) {ptBin = 1;}
    else if (pt >= 800) {ptBin = 2;}

    // Then check if jet score is greater than WP
    bool jetPassing;
    if (score > _wp) {jetPassing = true;}
    else {jetPassing = false;}

    float SF;

    // Now check year (this can prob all be done much more elegantly)
    if (_year == "16") {
	if (jetPassing) {
	    SF = SF2016_T[_variation][ptBin];
	}
	else {
            SF = SF2016_L[_variation][ptBin];
	}	
    }
    else if (_year == "16APV") {
        if (jetPassing) {
            SF = SF2016APV_T[_variation][ptBin];
        }
        else {
            SF = SF2016APV_L[_variation][ptBin];
        }
    }
    else if (_year == "17") {
        if (jetPassing) {
            SF = SF2017_T[_variation][ptBin];
        }
        else {
            SF = SF2017_L[_variation][ptBin];
        }
    }
    else {
        if (jetPassing) {
            SF = SF2018_T[_variation][ptBin];
        }
        else {
            SF = SF2018_L[_variation][ptBin];
        }
    }

    return SF;
}

RVec<float> PNetHbb_weight::eval(RVec<float> FatJet_pt, RVec<float> FatJet_eta,
				 RVec<float> FatJet_PNetHbbScore) {
    // Prepare the vector of weights to return
    RVec<float> out(3);
    std::string branchname;
    // Loop over variations (0:nom, 1:up, 2:down)
    for (size_t i=0; i<_variations.size(); i++) {
	_variation = _variations[i].first;  // 0, 1, or 2
	branchname = _variations[i].second; // string name

	// these are for Hbb tagging SFs
	float MC_tagged = 1.0,  MC_notTagged = 1.0, data_tagged = 1.0, data_notTagged = 1.0;
   	float data_tagged_up = 1.0, data_notTagged_up = 1.0;
   	float data_tagged_down = 1.0, data_notTagged_down = 1.0;

	float HbbTagEventWeight;

	// Loop over all jets in the event
	for (size_t ijet = 0; ijet<FatJet_pt.size(); ijet++) {
            float SF, eff;    // SF(pt,score), eff(pt, eta)
            float pt = FatJet_pt[ijet];
            float eta = FatJet_eta[ijet];
            float score = FatJet_PNetHbbScore[ijet];
	    // get SF and MC efficiency for this particular jet
	    SF = GetScaleFactor(pt, score);
	    eff = GetMCEfficiency(pt, eta);
	    if (score > _wp) {	// PASS
		MC_tagged *= eff;
		data_tagged *= SF*eff;
	    }
	    else {  // FAIL
		MC_notTagged *= (1. - eff);
		data_notTagged *= (1. - SF*eff);
	    }
	} // end loop over jets

        // Calculate the event weight for this variation
        HbbTagEventWeight = (data_tagged*data_notTagged) / (MC_tagged*MC_notTagged);
        out[i] = HbbTagEventWeight;

    } // end loop over variations

    // Send it {weight_nom, weight_up, weight_down}
    return out;
}

