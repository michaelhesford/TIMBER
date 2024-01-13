#pragma once
#include <ROOT/RVec.hxx>
#include <TRandom.h>
#include <string.h>
#include <vector>

using namespace ROOT::VecOps;

class btagSFHandler {
    public:
        int _nWps;
        RVec<float> _wps;
        RVec<float> _effs;
        TString _year;
        int _var;
        TRandom _randGen;
        std::vector<std::vector<float>> SF_L;
        std::vector<std::vector<float>> SF_T;

        btagSFHandler(RVec<float> wps,RVec<float> effs,TString year, int var);
        ~btagSFHandler(){};
        int createTaggingCategories(float taggerVal);
        int updateTaggingCategories(int jetCat,float jetPt);
        int bothLessThanOne(int jetCat, float sfT, float sfL);
        int bothGreaterThanOne(int jetCat, float sfT, float sfL);
        int LLowerTGreaterThanOne(int jetCat, float sfT, float sfL);
        int TLowerLGreaterThanOne(int jetCat, float sfT, float sfL);
        float getSF(float jetPt,int wpCat);
};

