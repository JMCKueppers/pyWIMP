/*****************************************************************************
 * Project: JKWimpDiffRateSFGPdf                                             *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to Savage/Freese/Gondolo 'Annual Modulation of Dark Matter in the Presence of Streams' eqn. (11) in earth's rest system*/

#ifndef _JKWimpDiffRateSFGPdf_hh_
#define _JKWimpDiffRateSFGPdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "TH1D.h"
#include "TF1.h"
#include "MGMVWimpFormFactor.hh"
 
class JKWimpDiffRateSFGPdf : public RooAbsPdf {
public:
  //JKWimpDiffRateSFGPdf() ; 
  JKWimpDiffRateSFGPdf(const char *name, const char *title,
        RooAbsReal& _constant_factors,           // Prefactor containing all factors independent of velocities
        RooAbsReal& _v_sub_E,                    // Earth velocity
        RooAbsReal& _v_sub_min,                  // Velocity corresponding to minimum energy which can give recoil energy E_R (recoil energy dependent)
        RooAbsReal& _v_sub_esc,                  // Escape velocity
        RooAbsReal& _v_sub_0,                    // Galactic base velocity
        MGMVWimpFormFactor& _form_factor_squared = MGMVWimpFormFactor::DefaultFormFactor());         // Form factor != 1 can be included
  JKWimpDiffRateSFGPdf(const JKWimpDiffRateSFGPdf& other, const char* name=0) ;         // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new JKWimpDiffRateSFGPdf(*this,newname); }
  inline virtual ~JKWimpDiffRateSFGPdf() { }
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

  ClassDef(JKWimpDiffRateSFGPdf,1)
};
 
#endif
