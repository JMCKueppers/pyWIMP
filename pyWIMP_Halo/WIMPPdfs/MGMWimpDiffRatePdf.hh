/*****************************************************************************
 * Project: MGMWimpDiffRatePdf                                               *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.12)
       Rate without escape velocity v_esc, with finite annual modulated velocity v_E of earth*/

#ifndef _MGMWimpDiffRatePdf_hh_
#define _MGMWimpDiffRatePdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "MGMVWimpFormFactor.hh"
 
class MGMWimpDiffRatePdf : public RooAbsPdf {
public:
  MGMWimpDiffRatePdf() {} ; 
  MGMWimpDiffRatePdf(const char *name, const char *title,
	      RooAbsReal& _v_sub_0,          // Base velocity of velocity distribution
	      RooAbsReal& _v_sub_min,        // Velocity corresponding to minimum energy which can give regoil energy E_R (recoil energy dependent)
	      RooAbsReal& _v_sub_E,          // Earth velocity (with annual modulation)
	      RooAbsReal& _R_sub_0,          // Base rate
	      RooAbsReal& _E_sub_0,          // Lewin/Smith E_0
	      RooAbsReal& _r,                // Lewin/Smith r
              MGMVWimpFormFactor& _form_factor = MGMVWimpFormFactor::DefaultFormFactor());   // Form factor != 1 can be included
  MGMWimpDiffRatePdf(const MGMWimpDiffRatePdf& other, const char* name=0) ;                  // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new MGMWimpDiffRatePdf(*this,newname); }       // Cloning with different name
  inline virtual ~MGMWimpDiffRatePdf() { }
  virtual ExtendMode extendMode() const { return MustBeExtended; }                                  // PDF is definded as extended PDF. Extended PDF is not normalized to 1.
  using RooAbsPdf::expectedEvents;                                                                  // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const 
    { Double_t temp = getNorm(nset); //cout << temp << endl;                                        // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}

protected:

  RooRealProxy v_sub_0 ;
  RooRealProxy v_sub_min ;
  RooRealProxy v_sub_E ;
  RooRealProxy R_sub_0 ;
  RooRealProxy E_sub_0 ;
  RooRealProxy r ;
  RooRealProxy form_factor ;
  Double_t fPrecisionCutOff ;
  
  //! Overload the Evaluate PDF when deriving from this class
  virtual Double_t EvaluatePDF() const;
  virtual Double_t EvaluateFormFactor() const;
  virtual Double_t evaluate() const ;

  virtual Double_t FindPrecisionCutoff();
private:

  ClassDef(MGMWimpDiffRatePdf,1)
};
 
#endif
