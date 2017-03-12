 /***************************************************************************** 
  * Project: MGMWimpHelmFFSquared                                             * 
  *                                                                           * 
  *****************************************************************************/ 

 /*Form factor according to Lewin/Smith 'Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil' eqn. (4.7)*/

 #include "MGMWimpHelmFFSquared.hh" 
 #include "RooAbsReal.h" 
 #include "RooAbsCategory.h" 
 #include <math.h> 
 #include "TMath.h" 

 ClassImp(MGMWimpHelmFFSquared) 

 MGMWimpHelmFFSquared::MGMWimpHelmFFSquared(const char *name, const char *title, 
                        RooAbsReal& _q,			// Momentum transfer
                        RooAbsReal& _r_sub_n,	// Effective nuclear radius
                        RooAbsReal& _s) :		// Nuclear skin thickness
   MGMVWimpFormFactor(name,title), 
   q("q","q",this,_q),
   r_sub_n("r_sub_n","r_sub_n",this,_r_sub_n),
   s("s","s",this,_s)
 { 
 } 

 // Constructor for cloning
 MGMWimpHelmFFSquared::MGMWimpHelmFFSquared(const MGMWimpHelmFFSquared& other, const char* name) :  
   MGMVWimpFormFactor(other,name), 
   q("q",this,other.q),
   r_sub_n("r_sub_n",this,other.r_sub_n),
   s("s",this,other.s)
 { 
 } 



 Double_t MGMWimpHelmFFSquared::evaluate() const 
 { 
   // Calculate form factor squared
   Double_t q_times_r = q*r_sub_n;
   Double_t temp = (q_times_r <= 2) ? 
            (1-TMath::Power(q_times_r,2)/5. + 3*TMath::Power(q_times_r,4)/175. // Series expansion of eqn. (4.3)^2 for q*r_n <= 2
             - 4*TMath::Power(q_times_r, 6)/4725.) :
            9*TMath::Power(                                                    // Eqn. (4.3)^2 for q*r_n > 2
              (TMath::Sin(q_times_r) - q_times_r*TMath::Cos(q_times_r))/
              (TMath::Power(q_times_r, 3)),2);
   return temp*TMath::Exp(-TMath::Power(q*s,2));                               // Eqn. (4.7)^2
 } 

