/*****************************************************************************
 * Project: MGMWimpClasses                                                   *
 *                                                                           *
 *****************************************************************************/
 
/*Base class for all classes describing a form factor.*/

#ifndef _MGMVWimpFormFactor_hh_
#define _MGMVWimpFormFactor_hh_ 1

#include "RooAbsPdf.h"
 
class MGMVWimpFormFactor : public RooAbsPdf {
public:
  MGMVWimpFormFactor(const char *name = "", const char *title = "");
  inline virtual ~MGMVWimpFormFactor() {}
  MGMVWimpFormFactor(const MGMVWimpFormFactor& other, const char* name=0) ;    // Constructor for cloning
  virtual TObject* clone(const char* newname) const { return new MGMVWimpFormFactor(*this,newname); }    // Cloning with different name
  
  using RooAbsPdf::expectedEvents;                                     // Get total number of events
  virtual Double_t expectedEvents(const RooArgSet* nset) const 
    { Double_t temp = getNorm(nset); //cout << temp << endl;           // Compute normalization factor from integrating PDF over arguments in ArgSet
      return temp;}
  // Force no normalization explicitly
  virtual Double_t getVal(const RooArgSet* /*nset*/ = 0) const 
    { return RooAbsPdf::getVal(0); }
  
  static MGMVWimpFormFactor& DefaultFormFactor();

protected:

  static MGMVWimpFormFactor fBasicFormFactor;
  Double_t evaluate() const {return 1;}

private:

  ClassDef(MGMVWimpFormFactor,1) 
};
 
#endif
