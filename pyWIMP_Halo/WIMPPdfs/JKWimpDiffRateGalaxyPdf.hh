/*****************************************************************************
 * Project: JKWimpDiffRateGalaxyPdf                                          *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to A. Green 'Dependence of direct detection signals on the WIMP velocity distribution' eqn. (2.1) with finite escape velocity and arbitrary velocity distribution in galactic rest system*/

#ifndef _JKWimpDiffRateGalaxyPdf_hh_
#define _JKWimpDiffRateGalaxyPdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TF1.h"
#include "MGMVWimpFormFactor.hh"
 
class JKWimpDiffRateGalaxyPdf : public RooAbsPdf {
public:
  //JKWimpDiffRateGalaxyPdf() ; 
  JKWimpDiffRateGalaxyPdf(const char *name, const char *title,
	      RooAbsReal& _constant_factors,           // Prefactor containing all factors independent of velocities
	      RooAbsReal& _v_sub_min,                  // Velocity corresponding to minimum energy which can give recoil energy E_R (recoil energy dependent)
	      RooAbsReal& _v_sub_esc,                  // Escape velocity
        TH1D& _velocity_distribution,            // Arbitrary velocity distribution as histogram
        MGMVWimpFormFactor& _form_factor_squared = MGMVWimpFormFactor::DefaultFormFactor());         // Form factor != 1 can be included
  JKWimpDiffRateGalaxyPdf(const JKWimpDiffRateGalaxyPdf& other, const char* name=0) ;         // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new JKWimpDiffRateGalaxyPdf(*this,newname); }
  inline virtual ~JKWimpDiffRateGalaxyPdf() { }
  virtual ExtendMode extendMode() const { return MustBeExtended; }     // PDF is definded as extended PDF. Extended PDF is not normalized to 1.
  using RooAbsPdf::expectedEvents;                                     // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const 
    { Double_t temp = getNorm(nset);                                   // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}

protected:

  RooRealProxy constant_factors;
  RooRealProxy v_sub_min ;
  RooRealProxy v_sub_esc ;
  TH1D velocity_distribution ;
  RooRealProxy form_factor_squared ;
  
  TF1* inverse_velocity ;
  TH1D velocity_distribution_temp ;
  
  
  virtual void SetDistribution() ;
  virtual Double_t GetVelocityIntegral() const ;
  virtual Double_t EvaluatePDF() const ;
  virtual Double_t EvaluateFormFactor() const;
  virtual Double_t evaluate() const ;

private:

  ClassDef(JKWimpDiffRateGalaxyPdf,1)
};
 
#endif
