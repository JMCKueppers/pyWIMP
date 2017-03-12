/*****************************************************************************
 * Project: MGMWimpDiffRateBasicPdf                                          *
 *                                                                           *
 *****************************************************************************/
 
 /*Pdf according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.14)
       Approximation for Rate without escape velocity v_esc, with finite annual modulated velocity v_E of earth. Approximation is independent of v_E*/

#ifndef _MGMWimpDiffRateBasicPdf_hh_
#define _MGMWimpDiffRateBasicPdf_hh_

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "MGMVWimpFormFactor.hh"
 
class MGMWimpDiffRateBasicPdf : public RooAbsPdf {
public:
  MGMWimpDiffRateBasicPdf() {} ; 
  MGMWimpDiffRateBasicPdf(const char *name, const char *title,
	      RooAbsReal& _R_sub_0,			// Base rate
	      RooAbsReal& _E_sub_0,			// Lewin/Smith E_0
          RooAbsReal& _Q,				// Recoil energy
	      RooAbsReal& _r,				// Lewin/Smith r
          MGMVWimpFormFactor& _form_factor = MGMVWimpFormFactor::DefaultFormFactor());		// Formfactor != 1 can be included
  MGMWimpDiffRateBasicPdf(const MGMWimpDiffRateBasicPdf& other, const char* name=0) ;     // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new MGMWimpDiffRateBasicPdf(*this,newname); }   // Cloning with different name
  inline virtual ~MGMWimpDiffRateBasicPdf() { }
  virtual ExtendMode extendMode() const { return MustBeExtended; }                                       // PDF is definded as extended PDF. Extended PDF is not normalized to 1.
  using RooAbsPdf::expectedEvents;                                                                       // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const                                           
    { Double_t temp = getNorm(nset); //cout << temp << endl;                                             // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}

protected:

  RooRealProxy R_sub_0 ;
  RooRealProxy E_sub_0 ;
  RooRealProxy Q ;
  RooRealProxy r ;
  RooRealProxy form_factor ;
  
  //! Overload the Evaluate PDF when deriving from this class
  virtual Double_t EvaluatePDF() const;
  virtual Double_t EvaluateFormFactor() const;
  virtual Double_t evaluate() const ;

private:

  ClassDef(MGMWimpDiffRateBasicPdf,1)
};
 
#endif
