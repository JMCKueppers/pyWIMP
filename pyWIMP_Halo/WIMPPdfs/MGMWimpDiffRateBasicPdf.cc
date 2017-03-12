 /***************************************************************************** 
  * Project: MGMWimpDiffRateBasicPdf                                          * 
  *                                                                           * 
  *****************************************************************************/ 

 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.14)
       Approximation for Rate without escape velocity v_esc, with finite annual modulated velocity v_E of earth. Approximation is independent of v_E*/

 #include "MGMWimpDiffRateBasicPdf.hh" 
 #include "RooAbsReal.h" 
 #include "TMath.h" 

 ClassImp(MGMWimpDiffRateBasicPdf) 

 MGMWimpDiffRateBasicPdf::MGMWimpDiffRateBasicPdf(const char *name, const char *title, 
                        RooAbsReal& _R_sub_0,                // Base rate
                        RooAbsReal& _E_sub_0,           	 // Lewin/Smith E_0    
                        RooAbsReal& _Q,                      // Recoil energy
                        RooAbsReal& _r,						 // Lewin/Smith r
                        MGMVWimpFormFactor& _form_factor) :  // Formfactor != 1 can be included
   RooAbsPdf(name,title), 
   R_sub_0("R_sub_0","R_sub_0",this,_R_sub_0),
   E_sub_0("E_sub_0","E_sub_0",this,_E_sub_0),
   Q("Q","Q",this,_Q),
   r("r","r",this,_r),
   form_factor("form_factor", "form_factor", this, _form_factor)
 { 
 } 

 // Constructor for cloning
 MGMWimpDiffRateBasicPdf::MGMWimpDiffRateBasicPdf(const MGMWimpDiffRateBasicPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   R_sub_0("R_sub_0",this,other.R_sub_0),
   E_sub_0("E_sub_0",this,other.E_sub_0),
   Q("Q",this,other.Q),
   r("r",this,other.r),
   form_factor("form_factor", this, other.form_factor)
 { 
 } 

 Double_t MGMWimpDiffRateBasicPdf::EvaluatePDF() const 
 { 
   // Rate without form factor 
   Double_t temp = 0.751*(R_sub_0/(E_sub_0*r))*
          TMath::Exp(-0.561*Q/(E_sub_0*r));
   return temp;
 } 

 Double_t MGMWimpDiffRateBasicPdf::EvaluateFormFactor() const
{
  return form_factor;
}

 Double_t MGMWimpDiffRateBasicPdf::evaluate() const
{
  // Form factor correction
  return EvaluateFormFactor()*EvaluatePDF();
}

