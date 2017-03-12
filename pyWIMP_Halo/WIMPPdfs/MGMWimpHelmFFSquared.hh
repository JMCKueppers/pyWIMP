/*****************************************************************************
 * Project: MGMWimpHelmFFSquared                                             *
 *                                                                           *
 *****************************************************************************/

/*Form factor according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (4.7)*/

#ifndef _MGMWimpHelmFFSquared_hh_
#define _MGMWimpHelmFFSquared_hh_ 1

#include "MGMVWimpFormFactor.hh"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
 
class MGMWimpHelmFFSquared : public MGMVWimpFormFactor {
public:
  MGMWimpHelmFFSquared() {} ; 
  MGMWimpHelmFFSquared(const char *name, const char *title,
	      RooAbsReal& _q,			// Momentum transfer
	      RooAbsReal& _r_sub_n,		// Effective nuclear radius
	      RooAbsReal& _s);			// Nuclear skin thickness
  MGMWimpHelmFFSquared(const MGMWimpHelmFFSquared& other, const char* name=0) ;   // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new MGMWimpHelmFFSquared(*this,newname); }     // Cloning with different name
  inline virtual ~MGMWimpHelmFFSquared() { }

protected:

  RooRealProxy q ;
  RooRealProxy r_sub_n ;
  RooRealProxy s ;
  
  Double_t evaluate() const ;

private:

  ClassDef(MGMWimpHelmFFSquared,1)
};
 
#endif
