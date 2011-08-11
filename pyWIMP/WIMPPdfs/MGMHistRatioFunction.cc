 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * This code was autogenerated by RooClassFactory                            * 
  *****************************************************************************/ 

 // Your description goes here... 

 #include "MGMHistRatioFunction.hh" 
 #include "RooRealVar.h" 
 #include "TH1.h" 

 ClassImp(MGMHistRatioFunction) 

 MGMHistRatioFunction::MGMHistRatioFunction(const char *name, const char *title,
	              const RooRealVar& _depVar,
	              RooAbsReal& _offsetMS,
	              RooAbsReal& _slopeMS,
	              RooAbsReal& _offsetSS,
	              RooAbsReal& _slopeSS,
	              RooAbsReal& _relative,
	              TH1&  _msHist,
	              TH1&  _ssHist
              ) :
   RooAbsReal(name,title), 
   fOffsetMS("fOffsetMS","fOffsetMS",this,_offsetMS),
   fSlopeMS("fSlopeMS","fSlopeMS",this,_slopeMS),
   fOffsetSS("fOffsetSS","fOffsetSS",this,_offsetSS),
   fSlopeSS("fSlopeSS","fSlopeSS",this,_slopeSS),
   fRelative("fRelative","fRelative",this,_relative),
   fMSHist(&_msHist),
   fSSHist(&_ssHist),
   fLowerLimit(_depVar.getMin()),
   fUpperLimit(_depVar.getMax())
 { 
 } 


 MGMHistRatioFunction::MGMHistRatioFunction(const MGMHistRatioFunction& other, const char* name) :  
   RooAbsReal(other, name), 
   fOffsetMS("fOffsetMS",this,other.fOffsetMS),
   fSlopeMS("fSlopeMS",this,other.fSlopeMS),
   fOffsetSS("fOffsetSS",this,other.fOffsetSS),
   fSlopeSS("fSlopeSS",this,other.fSlopeSS),
   fRelative("fRelative",this,other.fRelative),
   fMSHist(other.fMSHist),
   fSSHist(other.fSSHist),
   fLowerLimit(other.fLowerLimit),
   fUpperLimit(other.fUpperLimit)
 { 
 } 



 Double_t MGMHistRatioFunction::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   Double_t minEnergyMS = fLowerLimit*fSlopeMS + fOffsetMS;
   Double_t maxEnergyMS = fUpperLimit*fSlopeMS + fOffsetMS;
   Double_t minEnergySS = fLowerLimit*fSlopeSS + fOffsetSS;
   Double_t maxEnergySS = fUpperLimit*fSlopeSS + fOffsetSS;
        
   Int_t minEnergyBinMS = fMSHist->FindBin(minEnergyMS);
   Int_t maxEnergyBinMS = fMSHist->FindBin(maxEnergyMS);
   Int_t minEnergyBinSS = fSSHist->FindBin(minEnergySS);
   Int_t maxEnergyBinSS = fSSHist->FindBin(maxEnergySS);
   Double_t ratio = fMSHist->Integral(minEnergyBinMS, maxEnergyBinMS); 
   if (ratio == 0.0) return ratio;
   Double_t temp = fSSHist->Integral(minEnergyBinSS, maxEnergyBinSS);
   if (temp == 0.0) return 0.0;
   ratio /= temp;
   return ratio*fRelative/(1+ratio*fRelative);
 } 



