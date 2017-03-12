/*****************************************************************************
 * Project: JKWimpDiffRateLSPdf                                              *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections dark matter experiments based on elastic nuclear recoil' eqn. (3.13) in earth's rest system*/

#ifndef _JKWimpDiffRateLSPdf_hh_
#define _JKWimpDiffRateLSPdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TF1.h"
#include "MGMVWimpFormFactor.hh"
 
class JKWimpDiffRateLSPdf : public RooAbsPdf {
public:
  //JKWimpDiffRateLSPdf() ; 
  JKWimpDiffRateLSPdf(const char *name, const char *title,
        RooAbsReal& _constant_factors,           // Prefactor containing all factors independent of velocities
        RooAbsReal& _v_sub_E,                    // Earth velocity
        RooAbsReal& _v_sub_min,                  // Velocity corresponding to minimum energy which can give recoil energy E_R (recoil energy dependent)
        RooAbsReal& _v_sub_esc,                  // Escape velocity
        RooAbsReal& _v_sub_0,                    // Galactic base velocity
        MGMVWimpFormFactor& _form_factor_squared = MGMVWimpFormFactor::DefaultFormFactor());         // Form factor != 1 can be included
  JKWimpDiffRateLSPdf(const JKWimpDiffRateLSPdf& other, const char* name=0) ;         // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new JKWimpDiffRateLSPdf(*this,newname); }
  inline virtual ~JKWimpDiffRateLSPdf() { }
  virtual ExtendMode extendMode() const { return MustBeExtended; }     // PDF is definded as extended PDF. Extended PDF is not normalized to 1.
  using RooAbsPdf::expectedEvents;                                     // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const 
    { Double_t temp = getNorm(nset);                                   // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}

protected:

  RooRealProxy constant_factors;
  RooRealProxy v_sub_E ;
  RooRealProxy v_sub_min ;
  RooRealProxy v_sub_esc ;
  RooRealProxy v_sub_0 ;
  RooRealProxy form_factor_squared ;
  
  
  virtual Double_t GetVelocityIntegral() const ;
  virtual Double_t EvaluatePDF() const ;
  virtual Double_t EvaluateFormFactor() const;
  virtual Double_t evaluate() const ;

private:

  ClassDef(JKWimpDiffRateLSPdf,1)
};
 
#endif
