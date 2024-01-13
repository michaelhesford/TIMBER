#include <TFile.h>
#include <TMath.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "ROOT/RVec.hxx"
#include "../include/common.h"
#include "correction.h"

using namespace ROOT::VecOps;
using namespace ROOT::VecOps;
using rvec_i  = ROOT::VecOps::RVec<int>;
using rvec_f  = ROOT::VecOps::RVec<float>;
using LVector = ROOT::Math::PtEtaPhiMVector;

Int_t leadingNonGammaAK8Idx(Int_t nFatJet, rvec_f FatJet_eta, rvec_f FatJet_phi, Float_t gammaEta, Float_t gammaPhi)
//returns index of the leading AK8 jet, with delta R > 0.8 from provided eta/phi, i.e., away from the photon
//returns -1 if no such jet is found
{
    for(Int_t i=0; i<nFatJet;i++){
        if(hardware::DeltaR(LVector(1.,FatJet_eta[i],FatJet_phi[i],0.),LVector(1.,gammaEta,gammaPhi,0.))>0.8){
            return i;
        }
    }
    return -1;
}

Int_t genGammaPt(Int_t nGenPart,rvec_i GenPart_pdgId,rvec_f GenPart_pt,rvec_i GenPart_statusFlags)
//returns gen pt of the first gen gamma from hard process
//-1 otherwise
{
    Float_t pt = -1.;
    Int_t pid;
    for(Int_t i=0; i<nGenPart;i++){
        if (!(GenPart_statusFlags[i]&(1 << 7))){
            //If not hard process particle, continue
            continue;
        }
        pid = GenPart_pdgId[i];
        if(TMath::Abs(pid)==22){//gamma
            pt = GenPart_pt[i];
            return pt;
        }
    }
    return pt;
}



rvec_i taggedJetCount(Float_t FatJet_eta, Float_t FatJet_phi, Int_t nJet, rvec_f Jet_eta, Float_t eta_cut, rvec_f Jet_phi, rvec_f Jet_pt, rvec_f Jet_bTag, Float_t bTagCut, rvec_i Jet_hadronFlavour){
    //Loop over AK4 jets with DR>0.8 from given FatJet and satisfying pt>30, |eta|<eta_cut
    // Return number of light,c,b flavour jets passing a given tagger and cut value followed by total number of light,c,b flavour jets
    rvec_i out {0,0,0,0,0,0};
    //First three elements: number of light, c and b flavoured jets passing tagger cut
    //Last three: total number of light, c and b flavored jets
    Int_t kinematicFlag;
    Int_t btagFlag;
    Int_t DRflag;

    for(Int_t i=0; i<nJet;i++){
            DRflag = hardware::DeltaR(LVector(1.,FatJet_eta,FatJet_phi,0.),LVector(1.,Jet_eta[i],Jet_phi[i],0.))>0.8;
            kinematicFlag  = TMath::Abs(Jet_eta[i])<eta_cut && Jet_pt[i]>30;
            btagFlag       = Jet_bTag[i]>bTagCut;
            if (!(kinematicFlag && DRflag)){
               continue;
            }

            if(Jet_hadronFlavour[i]==0){
                //light
                out.at(3)++;
                if(btagFlag){
                    out.at(0)++;
                }
            }

            if(Jet_hadronFlavour[i]==4){
                //c
                out.at(4)++;
                if(btagFlag){
                    out.at(1)++;
                }
            }

            if(Jet_hadronFlavour[i]==5){
                //b
                out.at(5)++;
                if(btagFlag){
                    out.at(2)++;
                }
            }

    }

/*    for (int i = 0; i < 6; i++) {
        std::cout << out[i] << " ";
        }
    std::cout<<"\n";*/
    return out;
}

//float calcBtagSF(std::unique_ptr<correction::CorrectionSet> btvjson, Int_t nJet, rvec_f Jet_eta, rvec_f Jet_phi, rvec_f Jet_pt, rvec_i Jet_hadronFlavour){
rvec_f calcBtagWeight(correction::Correction::Ref correction_light, correction::Correction::Ref correction_bc,Float_t FatJet_eta, Float_t FatJet_phi,Int_t nJet, rvec_f Jet_eta, Float_t eta_cut, rvec_f Jet_phi, rvec_f Jet_pt, rvec_i Jet_hadronFlavour, Float_t eff_l, Float_t eff_c, Float_t eff_b){
    //returns btag sf correction for the event category with 0 medium b-tagged AK4 jets (for b-tag veto)
    //SF is based on AK4 jets which are away from the provided AK8 jet, i.e., the same jets which are considered for b tag veto

    rvec_f total_weight {1.,1.,1.};//nom, down, up
    Int_t kinematicFlag;
    Int_t DRflag;
    Float_t sf, sf_up, sf_dn, w, w_up, w_dn;

    for(Int_t i=0; i<nJet;i++){
            DRflag = hardware::DeltaR(LVector(1.,FatJet_eta,FatJet_phi,0.),LVector(1.,Jet_eta[i],Jet_phi[i],0.))>0.8;
            kinematicFlag  = TMath::Abs(Jet_eta[i])<eta_cut && Jet_pt[i]>30;
            if (!(kinematicFlag && DRflag)){
               continue;
            }

            if(Jet_hadronFlavour[i]==0){
                sf      = correction_light->evaluate({"central","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_dn   = correction_light->evaluate({"down","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_up   = correction_light->evaluate({"up","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});

                w       = (1 - eff_l)/(1 - eff_l*sf);
                w_dn    = (1 - eff_l)/(1 - eff_l*sf_dn);
                w_up    = (1 - eff_l)/(1 - eff_l*sf_up);

                total_weight.at(0)*= w;
                total_weight.at(1)*= w_dn;
                total_weight.at(2)*= w_up;

            }
            else if(Jet_hadronFlavour[i]==4){
                sf      = correction_bc->evaluate({"central","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_dn   = correction_bc->evaluate({"down","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_up   = correction_bc->evaluate({"up","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});

                w       = (1 - eff_c)/(1 - eff_c*sf);
                w_dn    = (1 - eff_c)/(1 - eff_c*sf_dn);
                w_up    = (1 - eff_c)/(1 - eff_c*sf_up);

                total_weight.at(0)*= w;
                total_weight.at(1)*= w_dn;
                total_weight.at(2)*= w_up;
            }
            else if(Jet_hadronFlavour[i]==5){
                sf      = correction_bc->evaluate({"central","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_dn   = correction_bc->evaluate({"down","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});
                sf_up   = correction_bc->evaluate({"up","M",Jet_hadronFlavour[i],abs(Jet_eta[i]),Jet_pt[i]});

                w       = (1 - eff_b)/(1 - eff_b*sf);
                w_dn    = (1 - eff_b)/(1 - eff_b*sf_dn);
                w_up    = (1 - eff_b)/(1 - eff_b*sf_up);

                total_weight.at(0)*= w;
                total_weight.at(1)*= w_dn;
                total_weight.at(2)*= w_up;
            }
            else{
                continue;
            }
/*    for (int i = 0; i < 3; i++) {
        std::cout << total_weight[i] << " ";
        }
    std::cout<<"\n";*/
    }//for(Int_t i=0; i<nJet;i++){
    return total_weight;
}