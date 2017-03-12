/*Base class for all classes describing a form factor.*/

#include "MGMVWimpFormFactor.hh"

ClassImp(MGMVWimpFormFactor)

MGMVWimpFormFactor MGMVWimpFormFactor::fBasicFormFactor;

MGMVWimpFormFactor::MGMVWimpFormFactor(const char *name, const char *title) : 
  RooAbsPdf(name, title) {}

// Constructor for cloning
MGMVWimpFormFactor::MGMVWimpFormFactor(const MGMVWimpFormFactor& other, const char* name) : 
  RooAbsPdf(other, name) {}

MGMVWimpFormFactor& MGMVWimpFormFactor::DefaultFormFactor()
{ return fBasicFormFactor; }
