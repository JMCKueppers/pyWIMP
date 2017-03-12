/***************************************************************************** 
 * Project: MGMWimpTimeFunction                                              * 
 *                                                                           * 
 *****************************************************************************/ 
 
/*Velocity modulation according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.6)
   Class implementing a generic time pdf for a WIMP search with the functional form: v_0+v_1*sin(2pi*t).
   This is useful for deriving sensitivity or looking for limits on a generic oscillation signal.*/

#ifndef _MGMWimpTimeFunction_hh_
#define _MGMWimpTimeFunction_hh_

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class MGMWimpTimeFunction : public RooAbsReal {
public:
  MGMWimpTimeFunction() {} ; 
  MGMWimpTimeFunction(const char *name, const char *title,
	      RooAbsReal& _velocity_0,		// Base component of modulating earth velocity
	      RooAbsReal& _velocity_1,		// Amplitude component of modulating earth velocity
	      RooAbsReal& _time,          // Time
	      RooAbsReal& _offset);       // Time Offset
  MGMWimpTimeFunction(const MGMWimpTimeFunction& other, const char* name=0) ;  // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new MGMWimpTimeFunction(*this,newname); }   // Cloning with different name
  inline virtual ~MGMWimpTimeFunction() { }

protected:

  RooRealProxy velocity_0 ;
  RooRealProxy velocity_1 ;
  RooRealProxy time ;
  RooRealProxy offset ;
  
  Double_t evaluate() const ;

private:

  ClassDef(MGMWimpTimeFunction,1)
};
 
#endif
