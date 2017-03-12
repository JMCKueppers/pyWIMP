import ROOT
from base_model import BaseModel
import pyWIMP_Halo.WIMPPdfs as pdfs

#For understanding all parameters, read
#Lewin/Smith 'Review of mathematics, numerical factors, and corrections
#for dark matter experiments based on elastic nuclear recoil'

class WIMPModel(BaseModel):     #A class wrapping several different WIMP PDFs
  def __init__(self, 
               basevars, 
               mass_of_wimp=20,                  #Mass of the WIMP (GeV)
               kilograms=1,                      #Mass of the detector (kg)
               constant_quenching=True,          #If true, use constant, non-energy-dependent quenching.
                                                 #Else use energy-dependent quenching.
               nucl_recoil = True):

    BaseModel.__init__(self, basevars)

    # constant quenching
    if nucl_recoil == False:
      if constant_quenching:
        #CoGeNT parameterization
        #self.quenching = ROOT.RooRealVar("quenching", "quenching", 0.2)
        #self.dQ_over_dE = ROOT.RooFormulaVar("dQ_over_dE", "#frac{dQ}{dE}",\
        #                                     "1./@0", ROOT.RooArgList(self.quenching))  #Q=quenching, E=recoil energy
        #self.recoil_energy = ROOT.RooFormulaVar("energy", "Energy", \
        #                                        "@0/@1", ROOT.RooArgList(basevars.get_energy(), \
        #                                        self.quenching))


        #Edelweiss parameterization of quenching
        self.quenching = ROOT.RooRealVar("quenching", "quenching", 0.16)
        self.dQ_over_dE = ROOT.RooFormulaVar("dQ_over_dE", "#frac{dQ}{dE}",\
                                             "1./@0", ROOT.RooArgList(self.quenching))  #Q=quenching, E=recoil energy
        self.recoil_energy = ROOT.RooFormulaVar("energy", "Energy", \
                                                "@0/@1", ROOT.RooArgList(basevars.get_energy(), \
                                                self.quenching))

      else:
        #CoGeNT parameterization
        #self.recoil_energy = ROOT.RooFormulaVar("energy", "Energy", \
        #                                        "4.03482*TMath::Power(@0,0.880165)", \
        #                                        ROOT.RooArgList(basevars.get_energy()))
        #self.dQ_over_dE = ROOT.RooFormulaVar("dQ_over_dE", "#frac{dQ}{dE}",\
        #                                     "3.55131*TMath::Power(@0, -0.119835)", \
        #                                     ROOT.RooArgList(basevars.get_energy()))  #Q=quenching, E=recoil energy

        #reparameterized for Edelweiss quenching factor

        self.recoil_energy = ROOT.RooFormulaVar("energy", "Energy", \
                                                "4.57458*TMath::Power(@0,0.84389)", \
                                                ROOT.RooArgList(basevars.get_energy()))
        self.dQ_over_dE = ROOT.RooFormulaVar("dQ_over_dE", "#frac{dQ}{dE}",\
                                             "3.86044*TMath::Power(@0, -0.11561)", \
                                             ROOT.RooArgList(basevars.get_energy()))  #Q=quenching, E=recoil energy

    #nuclear recoil energy scale - no quenching.
    else:
      self.quenching = ROOT.RooRealVar("quenching", "quenching", 1.0)
      self.dQ_over_dE = ROOT.RooFormulaVar("dQ_over_dE", "#frac{dQ}{dE}",\
                                           "1./@0", ROOT.RooArgList(self.quenching))  #Q=quenching, E=recoil energy
      self.recoil_energy = ROOT.RooFormulaVar("energy", "Energy", \
                                              "@0/@1", ROOT.RooArgList(basevars.get_energy(), \
                                              self.quenching))

    self.kilograms = ROOT.RooRealVar("kilograms", "kilograms", \
                                     kilograms)

    # Components of modulating earth velocity
    self.v_sub_E_sub_0 = ROOT.RooRealVar("v_sub_E_sub_0", \
                                         "Constant in Velocity Function", 232.59, "km s^-1")
    self.v_sub_E_sub_1 = ROOT.RooRealVar("v_sub_E_sub_1", \
                                         "Modulation Amplitude in Velocity Function", 0.001, \
                                         "km s^-1") 

    self.atomic_mass_of_target = ROOT.RooRealVar("atomic_mass_of_target", \
                                                 "Atomic Mass of Target", 72.6, "amu")  #Atomic mass A=68/0.932
    self.density_of_dark_matter = ROOT.RooRealVar("density_of_dark_matter", \
                                                  "Density of Dark Matter", 0.3, "Gev c^-2 cm^-3") 
    self.speed_of_light = ROOT.RooRealVar("speed_of_light", \
                                          "Speed of Light", 299792.458, "km s^-1")
    self.v_sub_0 = ROOT.RooRealVar("v_sub_0", \
                                   "Base Velocity", 220, "km s^-1")  #Base velocity of velocity-distribution
    self.v_sub_esc = ROOT.RooRealVar("v_sub_esc", \
                                     "Escape Velocity", 544, "km s^-1")
    self.mass_of_target = ROOT.RooFormulaVar("mass_of_target", \
                                             "Mass of Target", "0.932*@0", \
                                             ROOT.RooArgList(self.atomic_mass_of_target))  #Mass of target M_T=68 
    self.mass_of_target.setUnit("GeV c^02")


    # Following is for the Form Factors

    # Momentum transfer as calculated in Lewin/Smith sec. 4
    self.q = ROOT.RooFormulaVar("q", "Momentum Transfer",\
                                "sqrt(2*@0*@1)/197.3", ROOT.RooArgList(\
                                self.recoil_energy, self.mass_of_target))  #197.3=hbar
    self.q.setUnit("fm^-1")

    # Effective nuclear radius as calculated in Lewin/Smith sec. 4
    self.r_sub_n = ROOT.RooFormulaVar("r_sub_n", "Effective Nuclear Radius",\
                                      "1.14*TMath::Power(@0, 1./3.)", ROOT.RooArgList(\
                                      self.atomic_mass_of_target))
    self.r_sub_n.setUnit("fm")

    self.s = ROOT.RooRealVar("s", "Nuclear Skin Thickness",0.9)
    self.s.setUnit("fm")
    
    #self.r_sub_0 = ROOT.RooFormulaVar("r_sub_0", "Nuclear Radius",\
    #                                  "(0.3 + 0.91*TMath::Power(@0, 1./3.))", \
    #                                  ROOT.RooArgList(self.mass_of_target))
    #self.r_sub_0.setUnit("fm")
    #self.q_sub_0 = ROOT.RooFormulaVar("q_sub_0", "Coherence Energy",\
    #                                  "1.5*(197.3*197.3)/(@0*@1*@1)", \
    #                                  ROOT.RooArgList(self.mass_of_target,\
    #                                  self.r_sub_0))
    #self.q_sub_0.setUnit("keV")

    self.mass_of_wimp = ROOT.RooRealVar("mass_of_wimp", \
                                        "Mass of Wimp", mass_of_wimp, "GeV c^{-2}") 


    # The following takes into account the rate with days vs.
    # years and the kilogram mass of the detector
    # Be careful here, if time is constant be sure to take that into account:
    # Taking into account, that eq. (3.7) is given in tru. xs normalized to 1 pb
    if basevars.get_time().isConstant():
      time_dif = basevars.get_time().getMax() - basevars.get_time().getMin()
      # This is the time in units of years
      self.R_sub_0 = ROOT.RooFormulaVar("R_sub_0", "Base Rate",\
                                        "365*%f*@4*@5*503.4/(@0*@1)*(@2/0.4)*(@3/230.)" % time_dif, \
                                        #"503.4/(@0*@1)*(@2/0.4)*(@3/230.)", \
                                        ROOT.RooArgList(self.mass_of_target, self.mass_of_wimp,\
                                        self.density_of_dark_matter, self.v_sub_0,\
                                        self.kilograms, self.dQ_over_dE))
      self.R_sub_0.setUnit("pb^{-1}") 

    else:
      self.R_sub_0 = ROOT.RooFormulaVar("R_sub_0", "Base Rate",\
                                        "365*@4*@5*503.4/(@0*@1)*(@2/0.4)*(@3/230.)", \
                                        ROOT.RooArgList(self.mass_of_target, self.mass_of_wimp,\
                                        self.density_of_dark_matter, self.v_sub_0,\
                                        self.kilograms, self.dQ_over_dE))

      self.R_sub_0.setUnit("pb^{-1} yr^{-1}")

    # The following is dealing with the generation of the dR/dQ
    # NO escape velocity!


    self.r = ROOT.RooFormulaVar("r", "Lewin/Smith r",\
                                "4*@0*@1/((@0+@1)**2)", ROOT.RooArgList(\
                                self.mass_of_wimp, self.mass_of_target))

    self.E_sub_0 = ROOT.RooFormulaVar("E_sub_0", "Lewin/Smith E_sub_0",\
                                      "0.5e6*@0*((@1/@2)**2)", ROOT.RooArgList(\
                                      self.mass_of_wimp, self.v_sub_0, self.speed_of_light))
    # The following is for the total rate from Jungman, including
    # an exponential form factor


    # This if from Lewin, in particular: G.J. Alner et al. / Astroparticle Physics 23 (2005) p. 457, eqn. (13)
    # This is the conversion from sigma to normalized sigma per nucleon
    self.normalization = ROOT.RooFormulaVar("normalization",
                                            "Normalization Constant to WIMP-nucleon xs",
                                            #"(9.1e-3)*((1/@0)**2)/@1",
                                            #ROOT.RooArgList(self.atomic_mass_of_target,
                                            #self.r))
                                            #"((0.932/(@0*@1/(@0+@1)))**2)*(1/@2)**2",
                                            "(((0.932*@1/(0.932+@1))/(@0*@1/(@0+@1)))**2)*(1/@2)**2",
                                            ROOT.RooArgList(self.mass_of_target,
                                            self.mass_of_wimp, self.atomic_mass_of_target))
    self.normalization.setUnit("pb pb^{-1}")


    # Velocity of earth with annual modulation
    self.v_sub_E = pdfs.MGMWimpTimeFunction("v_sub_E", \
                                            "Velocity of the Earth", \
                                            self.v_sub_E_sub_0, self.v_sub_E_sub_1, \
                                            basevars.get_time(), basevars.get_time_offset())
    self.v_sub_E.setUnit( self.v_sub_E_sub_0.getUnit() )

    # Velocity corresponding to minimum energy which can give regoil energy E_R(recoil energy dependent)
    self.v_sub_min = ROOT.RooFormulaVar("v_sub_min", \
                                        "Minimum Velocity of Minimum Energy", \
                                        "sqrt(@0/(@1*@2))*@3", \
                                        ROOT.RooArgList(self.recoil_energy, self.E_sub_0, self.r,\
                                        self.v_sub_0))
    self.v_sub_min.setUnit( self.speed_of_light.getUnit() )
    
    # Woods-Saxon/Helm
    # This is the form-factor we use.
    self.woods_saxon_helm_ff_squared = pdfs.MGMWimpHelmFFSquared("woods_saxon_helm_ff_squared",\
                                                                 "Helm FF^{2} ",\
                                                                 self.q, self.r_sub_n, self.s)

    # Exponential 
    #self.exponential_ff_squared = ROOT.RooGenericPdf("exponential_ff_squared",\
                                                      #"Exponential Form Factor squared",\
                                                      #"exp(-@0/@1)",\
                                                      #ROOT.RooArgList(self.recoil_energy, self.q_sub_0))


    # Defining different WIMP models

    self.final_function = pdfs.MGMWimpDiffRatePdf("WIMPPDF_With_Time", \
                                                  "WIMP Pdf", \
                                                  self.v_sub_0, self.v_sub_min, \
                                                  self.v_sub_E, self.R_sub_0, \
                                                  self.E_sub_0, self.r, self.woods_saxon_helm_ff_squared)

    self.final_function_with_escape = pdfs.MGMWimpDiffRateEscapeVelPdf("WIMPPDF_With_Time_And_Escape_Vel", \
                                                                       "WIMP Pdf (esc velocity)", \
                                                                       self.v_sub_0, self.v_sub_min, \
                                                                       self.v_sub_E, self.R_sub_0, \
                                                                       self.E_sub_0, self.r, \
                                                                       self.v_sub_esc, self.woods_saxon_helm_ff_squared)

    self.final_function_with_escape_no_ff = pdfs.MGMWimpDiffRateEscapeVelPdf("WIMPPDF_With_Time_And_Escape_Vel", \
                                                                             "WIMP Pdf (esc velocity)", \
                                                                             self.v_sub_0, self.v_sub_min, \
                                                                             self.v_sub_E, self.R_sub_0, \
                                                                             self.E_sub_0, self.r, \
                                                                             self.v_sub_esc)


    self.simple_model = pdfs.MGMWimpDiffRateBasicPdf("simple model", 
                                                     "Lewin/Smith simple model",
                                                     self.R_sub_0,
                                                     self.E_sub_0,
                                                     self.recoil_energy,
                                                     self.r)#,
                                                     #self.woods_saxon_helm_ff_squared)


  # Get simple ( only 1-D exponential ) model ( returns RooAbsPdf )
  def get_simple_model(self):
    return self.simple_model

  # Get normalization : multiply cross-section by this to get
  # per - nucleon cross-section (in pb ) ( returns RooRealVar )
  def get_normalization(self): 
    return self.normalization


  # Get different form factors ( returns RooAbsPdf )

  def get_helm_form_factor(self):
    return self.woods_saxon_helm_ff_squared

  #def get_exponential_form_factor(self):
  #  return self.exponential_ff_squared

  # Get different 2-D WIMP models . The following allow
  # changing the value of the WIMP mass by passing
  # in an ( optional ) float parameter    

  # WIMP Model, no escape velocity, Helm form factor
  def get_WIMP_model(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function

  # WIMP Model , escape velocity , Helm form factor
  def get_WIMP_model_with_escape_vel(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_with_escape

  # WIMP Model , escape velocity , form factor = 1
  def get_WIMP_model_with_escape_vel_no_ff(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_with_escape_no_ff

  def get_model(self):
    return self.final_function_with_escape


