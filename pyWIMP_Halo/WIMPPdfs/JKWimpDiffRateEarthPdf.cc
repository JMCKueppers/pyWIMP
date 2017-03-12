/*****************************************************************************
 * Project: JKWimpDiffRateEarthPdf                                           *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to A. Green 'Dependence of direct detection signals on the WIMP velocity distribution' eqn. (2.1) with finite escape velocity and arbitrary velocity distribution in earth's rest system*/

#include "Riostream.h" 

#include "JKWimpDiffRateEarthPdf.hh" 
#include "RooAbsReal.h"
#include "TMath.h" 
#include "TH1D.h"
#include "TH3D.h"
#include "TF1.h"

ClassImp(JKWimpDiffRateEarthPdf) 

 //JKWimpDiffRateEarthPdf::JKWimpDiffRateEarthPdf() {} ;
 JKWimpDiffRateEarthPdf::JKWimpDiffRateEarthPdf(const char *name, const char *title, 
                        RooAbsReal& _constant_factors,            // Prefactor containing all factors independent of velocities
                        RooAbsReal& _v_sub_E_sub_x,               // Earth velocity
                        RooAbsReal& _v_sub_E_sub_y,
                        RooAbsReal& _v_sub_E_sub_z,
                        RooAbsReal& _v_sub_min,                   // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
                        RooAbsReal& _v_sub_esc,                   // Escape velocity
                        TH1D& _velocity_distribution,             // Arbitrary velocity distribution as histogram
                        MGMVWimpFormFactor& _form_factor_squared) :       // Form factor != 1 can be included
   RooAbsPdf(name,title), 
   constant_factors("constant_factors","constant_factors",this,_constant_factors),
   v_sub_E_sub_x("v_sub_E_sub_x","v_sub_E_sub_x",this,_v_sub_E_sub_x),
   v_sub_E_sub_y("v_sub_E_sub_y","v_sub_E_sub_y",this,_v_sub_E_sub_y),
   v_sub_E_sub_z("v_sub_E_sub_z","v_sub_E_sub_z",this,_v_sub_E_sub_z),
   v_sub_min("v_sub_min","v_sub_min",this,_v_sub_min),
   v_sub_esc("v_sub_esc","v_sub_esc",this,_v_sub_esc),
   velocity_distribution(_velocity_distribution),
   form_factor_squared("form_factor_squared","form_factor_squared",this,_form_factor_squared)
 { 
   SetDistribution();
 } 


 // Constructor for cloning
 JKWimpDiffRateEarthPdf::JKWimpDiffRateEarthPdf(const JKWimpDiffRateEarthPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   constant_factors("constant_factors",this,other.constant_factors),
   v_sub_E_sub_x("v_sub_E_sub_x",this,other.v_sub_E_sub_x),
   v_sub_E_sub_y("v_sub_E_sub_y",this,other.v_sub_E_sub_y),
   v_sub_E_sub_z("v_sub_E_sub_z",this,other.v_sub_E_sub_z),
   v_sub_min("v_sub_min",this,other.v_sub_min),
   v_sub_esc("v_sub_esc",this,other.v_sub_esc),
   velocity_distribution(other.velocity_distribution),
   form_factor_squared("form_factor_squared",this,other.form_factor_squared)
 { 
   SetDistribution();
 }


 void JKWimpDiffRateEarthPdf::SetDistribution()
 {
   // Divide by angle_element
   angle_element = new TF1("4pi*v^2", "4.*TMath::Power(x, 2.)*TMath::Pi()", velocity_distribution.GetXaxis()->GetXmin(), velocity_distribution.GetXaxis()->GetXmax());
   velocity_distribution_3d = velocity_distribution;
   velocity_distribution_3d.Divide(angle_element);
   
   NbinsV = velocity_distribution_3d.GetNbinsX();
   
   velocity_distribution_prime = TH3D("f^prime", "f in earth rest frame", NbinsV, velocity_distribution_3d.GetXaxis()->GetXmin(), velocity_distribution_3d.GetXaxis()->GetXmax(), NbinsPhi, 0., TMath::TwoPi(), NbinsTheta, 0., TMath::Pi());
   for (int i = 0; i < NbinsV; i++)
   {
    Double_t v = velocity_distribution_3d.GetBinLowEdge(i);
    for (int j = 0; j < NbinsPhi; j++)
    {
    Double_t phi = double(j)*TMath::TwoPi()/double(NbinsPhi-1);
     for (int k = 0; k < NbinsTheta; k++)
     {
      Double_t theta = double(k)*TMath::Pi()/double(NbinsTheta-1);
      
      Double_t v_prime = TMath::Sqrt(TMath::Power(v, 2.)+TMath::Power(v_sub_E_sub_x,2)+TMath::Power(v_sub_E_sub_y,2)+TMath::Power(v_sub_E_sub_z,2)+2.*v*((v_sub_E_sub_x*TMath::Cos(phi)+v_sub_E_sub_y*TMath::Sin(phi))*TMath::Sin(theta)+v_sub_E_sub_z*TMath::Cos(theta)));
      
      Double_t bin_content = velocity_distribution_3d.GetBinContent(velocity_distribution_3d.GetXaxis()->FindFixBin(v_prime));
      // Multiply with v*sin(theta) for integration: f'=f^~*v*sin(theta)
      bin_content *= v*TMath::Sin(theta);
      velocity_distribution_prime.SetBinContent(i, j, k, bin_content);
     }
    }
   }
 }


 // Calculate integral form v_min to v_max over (1/v)*velocity_distribution
 Double_t JKWimpDiffRateEarthPdf::GetVelocityIntegral() const
 {
   Int_t lower_limit = velocity_distribution_prime.GetXaxis()->FindBin(v_sub_min);
   Double_t velocity_integral = 0.;
   for (int i = 0; i < NbinsPhi; i++)
    {
     Double_t phi = double(i)*TMath::TwoPi()/double(NbinsPhi-1);
     for (int j = 0; j < NbinsTheta; j++)
      {
       Double_t theta = double(j)*TMath::Pi()/double(NbinsTheta-1);
       Double_t v_sub_max;
       if (TMath::Sqrt(TMath::Power(v_sub_E_sub_x,2)+TMath::Power(v_sub_E_sub_y,2)+TMath::Power(v_sub_E_sub_z,2)) == 0)
       {
        v_sub_max = v_sub_esc;
       }
       else
       {
        //v_sub_max = v_sub_esc+v_sub_E_sub_x*TMath::Cos(phi)*TMath::Sin(theta)+v_sub_E_sub_y*TMath::Sin(phi)*TMath::Sin(theta)+v_sub_E_sub_z*TMath::Cos(theta);
        v_sub_max = v_sub_esc+TMath::Abs(v_sub_E_sub_x*TMath::Cos(phi)*TMath::Sin(theta)+v_sub_E_sub_y*TMath::Sin(phi)*TMath::Sin(theta)+v_sub_E_sub_z*TMath::Cos(theta));
       }
       Int_t upper_limit = velocity_distribution_prime.GetXaxis()->FindBin(v_sub_max);
       if (lower_limit < upper_limit)
       {
	    for (int k = lower_limit; k <= upper_limit; k++)
         {
          Double_t bin_size = velocity_distribution_prime.GetXaxis()->GetBinWidth(k)*velocity_distribution_prime.GetYaxis()->GetBinWidth(i)*velocity_distribution_prime.GetZaxis()->GetBinWidth(j);
          velocity_integral += velocity_distribution_prime.GetBinContent(k,i,j)*bin_size;
         }
	   }
      }
    }
   return velocity_integral;
 }

 // Calculate differential rate
 Double_t JKWimpDiffRateEarthPdf::EvaluatePDF() const
 {
   return constant_factors*GetVelocityIntegral();
 }

 Double_t JKWimpDiffRateEarthPdf::EvaluateFormFactor() const
 {
   return form_factor_squared;
 }

 Double_t JKWimpDiffRateEarthPdf::evaluate() const
 {
   return EvaluatePDF()*EvaluateFormFactor();
 }



