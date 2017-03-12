/*****************************************************************************
 * Project: JKWimpDiffRateEarthPdf                                           *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to A. Green 'Dependence of direct detection signals on the WIMP velocity distribution' eqn. (2.1) with finite escape velocity and arbitrary velocity distribution in earth's rest system*/

#ifndef _JKWimpDiffRateEarthPdf_hh_
#define _JKWimpDiffRateEarthPdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TMath.h"
#include "MGMVWimpFormFactor.hh"
 
class JKWimpDiffRateEarthPdf : public RooAbsPdf {
public:
  //JKWimpDiffRateEarthPdf() ; 
  JKWimpDiffRateEarthPdf(const char *name, const char *title,
        RooAbsReal& _constant_factors,           // Prefactor containing all factors independent of velocities
        RooAbsReal& _v_sub_E_sub_x,              // Earth velocity
        RooAbsReal& _v_sub_E_sub_y,
        RooAbsReal& _v_sub_E_sub_z,
        RooAbsReal& _v_sub_min,                  // Velocity corresponding to minimum energy which can give recoil energy E_R (recoil energy dependent)
        RooAbsReal& _v_sub_esc,                  // Escape velocity
        TH1D& _velocity_distribution,            // Arbitrary velocity distribution as histogram
        MGMVWimpFormFactor& _form_factor_squared = MGMVWimpFormFactor::DefaultFormFactor());         // Form factor != 1 can be included
  JKWimpDiffRateEarthPdf(const JKWimpDiffRateEarthPdf& other, const char* name=0) ;         // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new JKWimpDiffRateEarthPdf(*this,newname); }
  inline virtual ~JKWimpDiffRateEarthPdf() { }
  virtual ExtendMode extendMode() const { return MustBeExtended; }     // PDF is definded as extended PDF. Extended PDF is not normalized to 1.
  using RooAbsPdf::expectedEvents;                                     // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const 
    { Double_t temp = getNorm(nset);                                   // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}

protected:

  RooRealProxy constant_factors;
  RooRealProxy v_sub_E_sub_x ;
  RooRealProxy v_sub_E_sub_y ;
  RooRealProxy v_sub_E_sub_z ;
  RooRealProxy v_sub_min ;
  RooRealProxy v_sub_esc ;
  TH1D velocity_distribution ;
  RooRealProxy form_factor_squared ;
  
  Int_t NbinsV ;
  static const Int_t NbinsPhi = 21;
  static const Int_t NbinsTheta = 11;
  
  TF1* angle_element ;
  TH1D velocity_distribution_3d ;
  TH3D velocity_distribution_prime ;
  
  virtual void SetDistribution() ;
  virtual Double_t GetVelocityIntegral() const ;
  virtual Double_t EvaluatePDF() const ;
  virtual Double_t EvaluateFormFactor() const ;
  virtual Double_t evaluate() const ;

private:

  ClassDef(JKWimpDiffRateEarthPdf,1)
};
 
#endif
