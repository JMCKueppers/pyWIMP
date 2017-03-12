/*****************************************************************************
 * Project: JKWimpDiffRateGalaxyPdf                                          *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to A. Green 'Dependence of direct detection signals on the WIMP velocity distribution' eqn. (2.1) with finite escape velocity and arbitrary velocity distribution in galactic rest system*/

#include "Riostream.h" 

#include "JKWimpDiffRateGalaxyPdf.hh" 
#include "RooAbsReal.h"
#include "TMath.h" 
#include "TH1D.h"
#include "TF1.h"

ClassImp(JKWimpDiffRateGalaxyPdf) 

 //JKWimpDiffRateGalaxyPdf::JKWimpDiffRateGalaxyPdf() {} ;
 JKWimpDiffRateGalaxyPdf::JKWimpDiffRateGalaxyPdf(const char *name, const char *title, 
                        RooAbsReal& _constant_factors,            // Prefactor containing all factors independent of velocities
                        RooAbsReal& _v_sub_min,                   // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
                        RooAbsReal& _v_sub_esc,                   // Escape velocity
                        TH1D& _velocity_distribution,             // Arbitrary velocity distribution as histogram
                        MGMVWimpFormFactor& _form_factor_squared) :       // Form factor != 1 can be included
   RooAbsPdf(name,title), 
   constant_factors("constant_factors","constant_factors",this,_constant_factors),
   v_sub_min("v_sub_min","v_sub_min",this,_v_sub_min),
   v_sub_esc("v_sub_esc","v_sub_esc",this,_v_sub_esc),
   velocity_distribution(_velocity_distribution),
   form_factor_squared("form_factor_squared","form_factor_squared",this,_form_factor_squared)
 { 
   SetDistribution();
 } 


 // Constructor for cloning
 JKWimpDiffRateGalaxyPdf::JKWimpDiffRateGalaxyPdf(const JKWimpDiffRateGalaxyPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   constant_factors("constant_factors",this,other.constant_factors),
   v_sub_min("v_sub_min",this,other.v_sub_min),
   v_sub_esc("v_sub_esc",this,other.v_sub_esc),
   velocity_distribution(other.velocity_distribution),
   form_factor_squared("form_factor_squared",this,other.form_factor_squared)
 { 
   SetDistribution();   
 } 


 void JKWimpDiffRateGalaxyPdf::SetDistribution()
 {
   // Inverse velocity function is needed for calculating velocity integral
   inverse_velocity = new TF1("v^-1", "TMath::Power(x, -1.)", velocity_distribution.GetXaxis()->GetXmin(), velocity_distribution.GetXaxis()->GetXmax());
   velocity_distribution_temp = velocity_distribution;
   velocity_distribution_temp.Multiply(inverse_velocity);
 }


 // Calculate integral form v_min to v_esc over (1/v)*velocity_distribution
 Double_t JKWimpDiffRateGalaxyPdf::GetVelocityIntegral() const
 {
   if (v_sub_min >= v_sub_esc)
   {
    return 0;
   }
   Int_t lower_limit = velocity_distribution.GetXaxis()->FindBin(v_sub_min);
   Int_t upper_limit = velocity_distribution.GetXaxis()->FindBin(v_sub_esc);
   Double_t velocity_integral = velocity_distribution_temp.Integral(lower_limit, upper_limit, "width");
   return velocity_integral;
 }

 // Calculate differential rate
 Double_t JKWimpDiffRateGalaxyPdf::EvaluatePDF() const 
 {
   return constant_factors*GetVelocityIntegral(); 
 }

 Double_t JKWimpDiffRateGalaxyPdf::EvaluateFormFactor() const
 {
   return form_factor_squared;
 }

 Double_t JKWimpDiffRateGalaxyPdf::evaluate() const
 {
   return EvaluatePDF()*EvaluateFormFactor();
 }



