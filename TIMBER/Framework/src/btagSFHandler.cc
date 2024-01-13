#include "TIMBER/Framework/include/btagSFHandler.h"

btagSFHandler::btagSFHandler(RVec<float> wps,RVec<float> effs,TString year, int var){
    //Only works for 2 wps!
    this->_wps  = wps;
    this->_effs = effs;
    this->_nWps = 2;
    this->_year = year;
    this->_var  = var;
    this->_randGen = TRandom(1234);

    if(year=="2016APV"){
        //https://coli.web.cern.ch/coli/.cms/btv/boohft-calib/20221201_bb_ULNanoV9_PNetXbbVsQCD_ak8_ext_2018/4_fit/
        //3 variations nom,dn,up
        //6 pT categories
        //our loose is actually medium
        SF_T = {{1.040,1.013,1.054,1.062,1.143,1.059},{0.979,0.969,0.986,0.98,1.058,0.917},{1.105,1.058,1.127,1.148,1.23,1.203}};
        SF_L = {{1.065,1.015,1.085,1.043,0.959,1.009},{0.993,0.964,1.024,0.962,0.882,0.898},{1.138,1.063,1.15,1.131,1.035,1.14}};
    }
    else if(year=="2016"){
        SF_T = {{1.011,1.06,1.041,1.034,1.066,1.086},{0.936,1.008,0.993,0.986,0.993,1.013},{1.084,1.113,1.091,1.086,1.143,1.164}};
        SF_L = {{0.990,0.998,1.022,1.024,1.113,1.090},{0.922,0.933,0.969,0.965,1.033,1.01},{1.056,1.064,1.079,1.092,1.199,1.17}};
    }
    else if(year=="2017"){
        SF_T = {{1.067,1.1,1.07,1.049,1.057,1.041},{0.971,1.01,1.013,0.99,0.994,0.991},{1.163,1.19,1.125,1.111,1.121,1.089}};
        SF_L = {{0.916,0.961,0.986,0.925,0.98,0.918},{0.829,0.888,0.918,0.837,0.908,0.866},{0.999,1.025,1.055,1.005,1.048,0.962}};    
    }
    else if(year=="2018"){
        SF_T = {{0.929,1.021,1.046,0.95,1.058,1.061},{0.882,0.962,0.994,0.894,1.009,1.018},{0.979,1.08,1.1,1.006,1.111,1.105}};
        SF_L = {{0.927,1.008,0.891,0.942,0.982,0.954},{0.875,0.926,0.814,0.883,0.927,0.871},{0.978,1.076,0.971,1.001,1.032,1.018}};
    }
}

int btagSFHandler::createTaggingCategories(float taggerVal){
    int cat = 0;
    for(int i=0; i<this->_nWps;i++){
        if(taggerVal>this->_wps[i]){
            cat = i+1;//cat 0 is fail, cat1 is wp1<tag<wp2
        }
    }
    return cat;
}



int btagSFHandler::updateTaggingCategories(int jetCat,float jetPt){
    //Method 2a) https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#2a_Jet_by_jet_updating_of_the_b
    float sfT   = this->getSF(jetPt,2);
    float sfL   = this->getSF(jetPt,1);
    int cat;
    //Case 1
    if(sfT<1 && sfL<1){
        //For L or T check if needs to be untagged
        cat = this->bothLessThanOne(jetCat,sfT,sfL);
        return cat;
    }
    //Case 2
    else if(sfT>1 && sfL>1){
        //Check if untagged gets promoted to T, if not, check if untagged needs to be promoted to L
        cat = this->bothGreaterThanOne(jetCat,sfT,sfL);
        return cat;
    }
    //Case 3
    else if(sfT<1 && sfL>1){
        //Check if T needs to be untagged, from *originally* untagged, check if need to be promoted to L
        cat = this->TLowerLGreaterThanOne(jetCat,sfT,sfL);
        return cat;         
        }
    //Case 4
    else if(sfT>1 && sfL<1){
        //If untagged, check if promoted to T. If loose, check if need to be demoted to untagged
        cat = this->LLowerTGreaterThanOne(jetCat,sfT,sfL);
        return cat; 
    }
    else{
    return jetCat;
    }
}

int btagSFHandler::LLowerTGreaterThanOne(int jetCat, float sfT, float sfL){
    float freq;
    float rand  = this->_randGen.Rndm();

    float effL  = this->_effs[0];
    float effT  = this->_effs[1];

    if(jetCat==0){
        float nom = effT*(sfT-1.);
        float denom = 1.-(effT+effL);
        freq = nom/denom;
        if(rand<freq){
            return 2;//upgrade to T
        }
        else{
            return 0;
        }
    }
    else if(jetCat==1){
        freq = 1.-sfL;
        if(rand<freq){
            return 0;
        }
        else{
            return 1;
        }
    }
    else{//don't change tight
        return jetCat;
    }
}

int btagSFHandler::TLowerLGreaterThanOne(int jetCat, float sfT, float sfL){
    float freq;
    float rand  = this->_randGen.Rndm();
    float effL  = this->_effs[0];
    float effT  = this->_effs[1];
    if(jetCat==2){
        freq = 1.-sfT;
    if(rand<freq){
        return 0;//demote to untagged
        }
    else{
        return 2;
        }
    }
    else if(jetCat==0){
        float effL  = this->_effs[0];
        float nom = effL*(sfL-1.);
        float denom = 1.-(effT+effL);
        freq = nom/denom;
        if(rand<freq){
            return 1;//update to L
        }
        else{
            return 0;//keep fail
        }
    }
    else{
        return jetCat;
    }
}


int btagSFHandler::bothLessThanOne(int jetCat, float sfT, float sfL){
    float freq;
    float rand  = this->_randGen.Rndm();
    if(jetCat==2){
        freq = 1.-sfT;
        if(rand<freq){
            return 0;//demote to untagged
        }
        else{
            return 2;
        }       
    }
    else if(jetCat==1){
        freq = 1.-sfL;
        if(rand<freq){
            return 0;//demote to untagged
        }
        else{
            return 1;
        }
    }       
    else{
        return jetCat;
    }
}


int btagSFHandler::bothGreaterThanOne(int jetCat, float sfT, float sfL){
    float freq_t, freq_l;
    float rand  = this->_randGen.Rndm();
    float effT  = this->_effs[1];
    float effL  = this->_effs[0];
    if(jetCat==0){
        float nom = effT*(sfT-1.);
        float denom = 1.-(effT+effL);
        freq_t = nom/denom;
        if(rand<freq_t){
            return 2;//update to T
        }
        else{
            rand  = this->_randGen.Rndm();
            nom   = effL*(sfL-1);
            denom = (1.-effL-effT)*(1-freq_t);
            freq_l= nom/denom;
            if(rand<freq_l){
                return 1;//promote to L
            }
            else{           
            return 0;//keep untagged
            }
        }
    }
    else{
        return jetCat;
    }
}

float btagSFHandler::getSF(float jetPt,int wpCat){
    int ptCat;
    int var = this->_var;
    //Deduce pT bin
    if(jetPt>600.){
        ptCat = 5;
        }
    else if(jetPt>500.){
        ptCat = 4;
    }
    else if(jetPt>450.){
        ptCat = 3;
    }
    else if(jetPt>400.){
        ptCat = 2;
    }
    else if(jetPt>350.){
        ptCat = 1;
    }
    else{
        ptCat = 0;
    }

    //Get SF
    if(wpCat==1){
        return this->SF_L[var][ptCat];          
    }
    else if(wpCat==2){
        return this->SF_T[var][ptCat];          
    }
    else{
        std::cout<<"Error getting SF\n";
        return -1;
    }

}