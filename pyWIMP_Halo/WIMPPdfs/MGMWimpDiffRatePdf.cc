 /***************************************************************************** 
  * Project: MGMWimpDiffRatePdf                                               * 
  *                                                                           * 
  *****************************************************************************/ 

 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.12)
       Rate without escape velocity v_esc, with finite annual modulated velocity v_E of earth*/

 #include "MGMWimpDiffRatePdf.hh" 
 #include "RooAbsReal.h" 
 #include "TMath.h" 

 using namespace std;

 ClassImp(MGMWimpDiffRatePdf) 

 MGMWimpDiffRatePdf::MGMWimpDiffRatePdf(const char *name, const char *title, 
                        RooAbsReal& _v_sub_0,                   // Base velocity of velocity-distribution
                        RooAbsReal& _v_sub_min,                 // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
                        RooAbsReal& _v_sub_E,                   // Earth velocity (with annual modulation)
                        RooAbsReal& _R_sub_0,                   // Base rate
                        RooAbsReal& _E_sub_0,                   // Lewin/Smith E_0
                        RooAbsReal& _r,                         // Lewin/Smith r
                        MGMVWimpFormFactor& _form_factor) :     // Form factor != 1 can be included
   RooAbsPdf(name,title), 
   v_sub_0("v_sub_0","v_sub_0",this,_v_sub_0),
   v_sub_min("v_sub_min","v_sub_min",this,_v_sub_min),
   v_sub_E("v_sub_E","v_sub_E",this,_v_sub_E),
   R_sub_0("R_sub_0","R_sub_0",this,_R_sub_0),
   E_sub_0("E_sub_0","E_sub_0",this,_E_sub_0),
   r("r","r",this,_r),
   form_factor("form_factor", "form_factor", this, _form_factor)
 { 
   fPrecisionCutOff = FindPrecisionCutoff();
 } 


 // Constructor for cloning
 MGMWimpDiffRatePdf::MGMWimpDiffRatePdf(const MGMWimpDiffRatePdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   v_sub_0("v_sub_0",this,other.v_sub_0),
   v_sub_min("v_sub_min",this,other.v_sub_min),
   v_sub_E("v_sub_E",this,other.v_sub_E),
   R_sub_0("R_sub_0",this,other.R_sub_0),
   E_sub_0("E_sub_0",this,other.E_sub_0),
   r("r",this,other.r),
   form_factor("form_factor", this, other.form_factor)
 { 
   fPrecisionCutOff = FindPrecisionCutoff();
 } 

Double_t MGMWimpDiffRatePdf::FindPrecisionCutoff() 
{

  /* This functions determines when we switch from using erfs to erfcs. */
  Double_t lower_limit = 0;
  Double_t upper_limit = 0;
  // Lower limit will be almost zero, but for posterity:
  for (Double_t i=1; i>1e-100;i*=0.1) {
    Double_t diff = TMath::Erfc(i);
    if (diff == 1.0) {
      lower_limit = i;
      break;
    }
  }
  // The upper limit is really what counts 
  for (Double_t i=3; i<10;i+=0.000001) {
    Double_t diff = TMath::Erf(i);
    if (diff == 1.0) {
      upper_limit = i;
      break;
    }
  }

  if ( upper_limit == 0 ) {
    cout << "Failed to find precision!" << endl;
    upper_limit = 100;
  } 
  // Naively choose a point right in the middle:
  return (upper_limit - lower_limit)/2.;

}


 Double_t MGMWimpDiffRatePdf::EvaluatePDF() const 
 { 
   Double_t erf_value_one = ( v_sub_min + v_sub_E )/v_sub_0;
   Double_t erf_value_two = ( v_sub_min - v_sub_E )/v_sub_0;
   Double_t erf_val = (erf_value_one >= fPrecisionCutOff) ?
       (TMath::Erfc(erf_value_two) - TMath::Erfc(erf_value_one)) : // Two high, use erfcs. erfc(a)-erfc(b)=erf(b)-erf(a)
       (TMath::Erf(erf_value_one) - TMath::Erf(erf_value_two)); // Default, use erf

   // Rate without form factor
   Double_t temp = (R_sub_0/(E_sub_0*r))*(TMath::Sqrt(TMath::Pi())/4.)*
          (v_sub_0/v_sub_E)*erf_val;
   return temp;
 } 

 Double_t MGMWimpDiffRatePdf::EvaluateFormFactor() const
{
  return form_factor;
}

 Double_t MGMWimpDiffRatePdf::evaluate() const
{
  // Form factor correction
  Double_t func = EvaluateFormFactor()*EvaluatePDF();
  return ( func < 1e-15 ) ? 1e-15 : func;
}

