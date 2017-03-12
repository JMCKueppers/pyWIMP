 /***************************************************************************** 
  * Project: MGMWimpDiffRateEscapeVelPdf                                      * 
  *                                                                           * 
  *****************************************************************************/ 

 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.13)
       Rate with escape velocity v_esc, with finite annual modulated velocity v_E of earth*/

 #include "MGMWimpDiffRateEscapeVelPdf.hh" 
 #include "RooAbsReal.h" 
 #include "TMath.h" 

 ClassImp(MGMWimpDiffRateEscapeVelPdf) 

 MGMWimpDiffRateEscapeVelPdf::MGMWimpDiffRateEscapeVelPdf(const char *name, const char *title, 
                        RooAbsReal& _v_sub_0,                   // Base velocity of velocity-distribution
                        RooAbsReal& _v_sub_min,                 // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
                        RooAbsReal& _v_sub_E,                   // Earth velocity (with annual modulation)
                        RooAbsReal& _R_sub_0,                   // Base rate
                        RooAbsReal& _E_sub_0,                   // Lewin/Smith E_0
                        RooAbsReal& _r,                         // Lewin/Smith r
	                RooAbsReal& _v_sub_esc,                       // Escape velocity
	                MGMVWimpFormFactor& _form_factor) :           // Form factor != 1 can be included
   MGMWimpDiffRatePdf(name, title,
                      _v_sub_0, _v_sub_min, _v_sub_E,
                      _R_sub_0, _E_sub_0, _r, _form_factor), 
   v_sub_esc("v_sub_esc","v_sub_esc",this,_v_sub_esc)
   
 { 
 } 


 // Constructor for cloning
 MGMWimpDiffRateEscapeVelPdf::MGMWimpDiffRateEscapeVelPdf(const MGMWimpDiffRateEscapeVelPdf& other, const char* name) :  
   MGMWimpDiffRatePdf(other,name), 
   v_sub_esc("v_sub_esc",this,other.v_sub_esc)
 { 
 } 



 Double_t MGMWimpDiffRateEscapeVelPdf::EvaluatePDF() const 
 { 
   Double_t temp = 
          getK0OverK1()*( MGMWimpDiffRatePdf::EvaluatePDF() - 
          (R_sub_0/(E_sub_0*r))*
          TMath::Exp(-TMath::Power(v_sub_esc/v_sub_0,2)) );
   return temp; 
 } 

 Double_t MGMWimpDiffRateEscapeVelPdf::getK0OverK1() const 
 { 
   // Calculate k_0/k_1 according to eqn. (2.2)
   return 1./(TMath::Erf(v_sub_esc/v_sub_0) - (2./TMath::Sqrt(TMath::Pi()))*
          (v_sub_esc/v_sub_0)*TMath::Exp(-TMath::Power(v_sub_esc/v_sub_0,2))); 
 } 



