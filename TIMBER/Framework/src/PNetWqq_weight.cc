#include "../include/PNetWqq_weight.h"

// Constructor - set year and load eff map in memory, 
PNetWqq_weight::PNetWqq_weight(std::string year, std::string effmapname, float wp, int target_flavor, std::set<int> other):_effmapname(effmapname),_year(year),_wp(wp),_target_flavor(target_flavor),_other(other) {
    // I am assuming that it is one ROOT file per year, with one histogram.
    _effroot = hardware::Open(_effmapname, false);
    _effmap = (TH2F*)_effroot->Get("ratio");
}

// Close files on destruction
PNetWqq_weight::~PNetWqq_weight() {
    _effroot->Close();
}

// Get the MC effieciency of the jets in the given (pT, eta) bin
float PNetWqq_weight::GetMCEfficiency(float pt, float eta, int flavor) {
    // Efficiency map binned in pT: [60,0,3000], eta: [24,-2.4,2.4]
    int xbin = (int)(pt*30./3000.);
    int ybin = (int)((eta+2.4)*24/4.8);
    return _effmap->GetBinContent(xbin,ybin);
}

// Get the SF for the jet based on its pT and score and (internal) variation
float PNetWqq_weight::GetScaleFactor(float pt, int flavor, int variation) {
    // First check which pT bin we're in
    int ptBin;
    if (pt >= 300 && pt <= 450) {ptBin = 0;}
    else if (pt > 450 && pt <= 600) {ptBin = 1;}
    else if (pt > 600) {ptBin = 2;}

    float SF;

    // Check if selecting scale factor for "other" jets
    if (flavor != _target_flavor) {
        if (_target_flavor == -1 && _other.count(flavor)) {
            if (_year == "16") {SF = SF2016_o[ptBin][variation];}
            else if (_year == "16APV") {SF = SF2016APV_o[ptBin][variation];}
            else if (_year == "17") {SF = SF2017_o[ptBin][variation];}
            else {SF = SF2018_o[ptBin][variation];}
        }
        else {
            SF = 1;
        }
    }  
    else {
        if (flavor == 0) { // Get medium working point SF for merged top jet
            if (_year == "16") {SF = SF2016_t[ptBin][variation];}	
            else if (_year == "16APV") {SF = SF2016APV_t[ptBin][variation];}       
            else if (_year == "17") {SF = SF2017_t[ptBin][variation];}       
            else {SF = SF2018_t[ptBin][variation];}
        }

        else if (flavor == 1) { // Get medium working point SF for merged w jet
            if (_year == "16") {SF = SF2016_w[ptBin][variation];}
            else if (_year == "16APV") {SF = SF2016APV_w[ptBin][variation];}
            else if (_year == "17") {SF = SF2017_w[ptBin][variation];}
            else {SF = SF2018_w[ptBin][variation];}
        }
   
        else {SF = 1;}
    }
    return SF;
}

RVec<float> PNetWqq_weight::eval(float Wqq_pt, float Wqq_eta, float Wqq_PNetWqqScore, int Wqq_jetFlavor) {
    // Prepare the vector of weights to return
    RVec<float> out(3);
    // Loop over variations (0:nom, 1:up, 2:down)
    for (int var : {0,1,2}) {
	// these are for Wqq tagging SFs
	float MC_tagged = 1.0,  MC_notTagged = 1.0;
        float data_tagged = 1.0, data_notTagged = 1.0;

	float WqqTagEventWeight; // final multiplicative factor
        float SF, eff;    // SF(pt,score), eff(pt, eta)
        // get SF and MC efficiency for this particular jet
        SF = GetScaleFactor(Wqq_pt, Wqq_jetFlavor, var);
        eff = GetMCEfficiency(Wqq_pt, Wqq_eta, Wqq_jetFlavor);
        if (Wqq_PNetWqqScore > _wp) {  // PASS
            //MC_tagged *= eff;
            //data_tagged *= SF*eff;
            data_tagged *= SF;
        }
        else {  // FAIL
            if (eff == 1) {eff = 0.99;} // Prevent the event weight from becoming undefined
            MC_notTagged *= (1 - eff);
            data_notTagged *= (1 - SF*eff);
        }
        // Calculate the event weight for this variation
        WqqTagEventWeight = (data_tagged*data_notTagged) / (MC_tagged*MC_notTagged);
        out[var] = WqqTagEventWeight;

    } // end loop over variations

    // Send it {weight_nom, weight_up, weight_down}
    return out;
}

