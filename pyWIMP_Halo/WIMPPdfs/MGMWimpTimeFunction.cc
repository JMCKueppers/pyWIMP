 /***************************************************************************** 
  * Project: MGMWimpTimeFunction                                              * 
  *                                                                           * 
  *****************************************************************************/ 
  
 /*Velocity modulation according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (3.6)
    Class implementing a generic time pdf for a WIMP search with the functional form: v_0+v_1*sin(2pi*t).
    This is useful for deriving sensitivity or looking for limits on a generic oscillation signal.*/

 #include "MGMWimpTimeFunction.hh" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 
 #include <math.h> 
 #include "TMath.h" 

 ClassImp(MGMWimpTimeFunction) 

 MGMWimpTimeFunction::MGMWimpTimeFunction(const char *name, const char *title, 
                        RooAbsReal& _velocity_0,	// Base component of modulating earth velocity
                        RooAbsReal& _velocity_1,	// Amplitude component of modulating earth velocity
                        RooAbsReal& _time,        // Time
                        RooAbsReal& _offset) :    // Time Offset
   RooAbsReal(name,title), 
   velocity_0("velocity_0","velocity_0",this,_velocity_0),
   velocity_1("velocity_1","velocity_1",this,_velocity_1),
   time("time","time",this,_time),
   offset("offset","offset",this,_offset)
 { 
 } 

 // Constructor for Cloning
 MGMWimpTimeFunction::MGMWimpTimeFunction(const MGMWimpTimeFunction& other, const char* name) :  
   RooAbsReal(other,name), 
   velocity_0("velocity_0",this,other.velocity_0),
   velocity_1("velocity_1",this,other.velocity_1),
   time("time",this,other.time),
   offset("offset",this,other.offset)
 { 
 } 



 Double_t MGMWimpTimeFunction::evaluate() const 
 { 
   // Calculate modulating earth velocity
   return velocity_0 + velocity_1*TMath::Sin(TMath::TwoPi()*(time + offset)); //set offset to -61./365.25 in order for this to start on Jan 1
 } 



