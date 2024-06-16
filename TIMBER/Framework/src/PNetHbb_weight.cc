#include "../include/PNetHbb_weight.h"

// Constructor - set year and load eff map in memory, 
PNetHbb_weight::PNetHbb_weight(std::string year, std::string effmapname, float wp, int target_flavor, std::set<int> other):_effmapname(effmapname),_year(year),_wp(wp),_target_flavor(target_flavor),_other(other) {
    // Change this logic yourself based on how you made the efficiency map.
    // I am assuming that it is one ROOT file per year, with one histogram.
    _effroot = hardware::Open(_effmapname, false);
    _effmap = (TH2F*)_effroot->Get("ratio");
}

// Close files on destruction
PNetHbb_weight::~PNetHbb_weight() {
    _effroot->Close();
}

// Get the MC effieciency of the jets in the given (pT, eta) bin
float PNetHbb_weight::GetMCEfficiency(float pt, float eta, int flavor) {
    // Efficiency map binned in pT: [60,0,3000], eta: [24,-2.4,2.4]
    int xbin = (int)(pt*30./3000.);
    int ybin = (int)((eta+2.4)*24/4.8);
    return _effmap->GetBinContent(xbin,ybin);
}

// Get the SF for the jet based on its pT and score and (internal) variation
std::vector<float> PNetHbb_weight::GetScaleFactors(float pt, int flavor) {
    // Output: {SF, SF_up, SF_down}
    // First check which pT bin we're in
    int ptBin;
    if (pt >= 300 && pt < 450) {ptBin = 0;}
    else if (pt >= 450 && pt < 600) {ptBin = 1;}
    else if (pt >= 600) {ptBin = 2;}

    std::vector<float> SF;

    if (flavor != _target_flavor) {
        SF = {1,1,1};
    }

    else { 
        if (flavor == 4 && pt > 400) { // SF's for merged Higgs jet
	    // separate pt bins for these SF's
	    if (pt >= 400 && pt < 600) {ptBin = 0;}
       	    else if (pt >= 600 && pt < 800) {ptBin = 1;}
	    else if (pt >= 800) {ptBin = 2;}

 	    // Now check year (this can prob all be done much more elegantly)
	    if (_year == "16") {SF = SF2016_h[ptBin];}
	    else if (_year == "16APV") {SF = SF2016APV_h[ptBin];}
	    else if (_year == "17") {SF = SF2017_h[ptBin];}
            else {SF = SF2018_h[ptBin];} // _year = "18"
        }

        else if (flavor == 0) { // SF's for merged top jet
            if (_year == "16") {SF = SF2016_t[ptBin];}
            else if (_year == "16APV") {SF = SF2016APV_t[ptBin];}
            else if (_year == "17") {SF = SF2017_t[ptBin];}
            else {SF = SF2018_t[ptBin];} // _year = "18"
        }

        else if (flavor == 2) { // SF's for merged bq jet
            if (_year == "16") {SF = SF2016_bq[ptBin];}
            else if (_year == "16APV") {SF = SF2016APV_bq[ptBin];}
            else if (_year == "17") {SF = SF2017_bq[ptBin];}
            else {SF = SF2018_bq[ptBin];} // _year = "18"
        }

        else {SF = {1,1,1};} 

    }
    return SF;
}

RVec<float> PNetHbb_weight::eval(float Higgs_pt, float Higgs_eta, float Higgs_PNetHbbScore, int Higgs_jetFlavor) {
    // Prepare the vector of weights to return
    RVec<float> out(3);

    float HbbTagEventWeight; // final multiplicative factor
    float eff;               // eff(pt, eta, flavor)
    std::vector<float> SF;       // SF(pt, flavor), {SF, SF_up, SF_down}
    // get SF's and MC efficiency for this particular jet
    SF = GetScaleFactors(Higgs_pt, Higgs_jetFlavor);
    eff = GetMCEfficiency(Higgs_pt, Higgs_eta, Higgs_jetFlavor);
    for (int i : {0,1,2}) { // {nominal, up, down}
        float MC_tagged = 1.0,  MC_notTagged = 1.0;
        float data_tagged = 1.0, data_notTagged = 1.0;
        if (Higgs_PNetHbbScore > _wp) {  // PASS
            //MC_tagged *= eff;
            //data_tagged *= SF[i]*eff;
            data_tagged *= SF[i];
        }
        else {  // FAIL
            if (eff == 1) {eff = 0.99;} // Prevent the event weight from becoming undefined
            MC_notTagged *= (1 - eff);
            data_notTagged *= (1 - SF[i]*eff);
        }
        // Calculate the event weight for this variation
        HbbTagEventWeight = (data_tagged*data_notTagged) / (MC_tagged*MC_notTagged);
        out[i] = HbbTagEventWeight;

    } // end loop over variations

    // Send it {weight_nom, weight_up, weight_down}
    return out;
}

