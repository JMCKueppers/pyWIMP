/*****************************************************************************
 * Project: JKWimpDiffRateSFGPdf                                             *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to Savage/Freese/Gondolo 'Annual Modulation of Dark Matter in the Presence of Streams' eqn. (11) in earth's rest system*/

#include "Riostream.h" 

#include "JKWimpDiffRateSFGPdf.hh" 
#include "RooAbsReal.h"
#include "TMath.h" 
#include "TH1D.h"
#include "TF1.h"

ClassImp(JKWimpDiffRateSFGPdf) 

 //JKWimpDiffRateSFGPdf::JKWimpDiffRateSFGPdf() {} ;
 JKWimpDiffRateSFGPdf::JKWimpDiffRateSFGPdf(const char *name, const char *title, 
                        RooAbsReal& _constant_factors,            // Prefactor containing all factors independent of velocities
                        RooAbsReal& _v_sub_E,                     // Earth velocity
                        RooAbsReal& _v_sub_min,                   // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
                        RooAbsReal& _v_sub_esc,                   // Escape velocity
                        RooAbsReal& _v_sub_0,                     // Galactic base velocity
                        MGMVWimpFormFactor& _form_factor_squared) :       // Form factor != 1 can be included
   RooAbsPdf(name,title), 
   constant_factors("constant_factors","constant_factors",this,_constant_factors),
   v_sub_E("v_sub_E","v_sub_E",this,_v_sub_E),
   v_sub_min("v_sub_min","v_sub_min",this,_v_sub_min),
   v_sub_esc("v_sub_esc","v_sub_esc",this,_v_sub_esc),
   v_sub_0("v_sub_0","v_sub_0",this,_v_sub_0),
   form_factor_squared("form_factor_squared","form_factor_squared",this,_form_factor_squared)
 { 
 } 


 // Constructor for cloning
 JKWimpDiffRateSFGPdf::JKWimpDiffRateSFGPdf(const JKWimpDiffRateSFGPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   constant_factors("constant_factors",this,other.constant_factors),
   v_sub_E("v_sub_E",this,other.v_sub_E),
   v_sub_min("v_sub_min",this,other.v_sub_min),
   v_sub_esc("v_sub_esc",this,other.v_sub_esc),
   v_sub_0("v_sub_0",this,other.v_sub_0),
   form_factor_squared("form_factor_squared",this,other.form_factor_squared)
 { 
 } 


 // Calculate integral form v_min to v_esc over (1/v)*velocity_distribution
 Double_t JKWimpDiffRateSFGPdf::GetVelocityIntegral() const
 {
   Double_t x = v_sub_min/v_sub_0;
   Double_t y = v_sub_E/v_sub_0;
   Double_t z = v_sub_esc/v_sub_0;
   if (z<y and x<TMath::Abs(y-z))
   {
    return 1./v_sub_E;
   }
   else if (z>y and x<TMath::Abs(y-z))
   {
    return (TMath::Erf(x+y)-TMath::Erf(x-y)-4*TMath::Power(TMath::Pi(),-0.5)*y*TMath::Exp(-TMath::Power(z,2.)))/(2*v_sub_E*(TMath::Erf(z)-2*TMath::Power(TMath::Pi(),-0.5)*z*TMath::Exp(-TMath::Power(z,2.))));
   }
   else if (TMath::Abs(y-z)<x and x<y+z)
   {
    return (TMath::Erf(z)-TMath::Erf(x-y)-2*TMath::Power(TMath::Pi(),-0.5)*(y+z-x)*TMath::Exp(-TMath::Power(z,2.)))/(2*v_sub_E*(TMath::Erf(z)-2*TMath::Power(TMath::Pi(),-0.5)*z*TMath::Exp(-TMath::Power(z,2.))));
   }
   else if (y+z<x)
   {
    return 0;
   }
   return 0;
 }

 // Calculate differential rate
 Double_t JKWimpDiffRateSFGPdf::EvaluatePDF() const 
 {
   return constant_factors*GetVelocityIntegral(); 
 }

 Double_t JKWimpDiffRateSFGPdf::EvaluateFormFactor() const
 {
   return form_factor_squared;
 }

 Double_t JKWimpDiffRateSFGPdf::evaluate() const
 {
   return EvaluatePDF()*EvaluateFormFactor();
 }



