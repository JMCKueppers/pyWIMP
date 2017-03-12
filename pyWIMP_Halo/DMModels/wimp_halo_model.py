import ROOT
import pyWIMP_Halo.WIMPPdfs as pdfs

from base_model import BaseModel

class WIMPHaloModel(BaseModel):     #A class wrapping several different WIMP PDFs
  def __init__(self,
               basevars,
               velocity_distribution,    #Arbitrary velocity distribution as histogram
               mass_of_wimp=20):         #Mass of the WIMP (GeV)

    BaseModel.__init__(self, basevars)

    self.velocity_distribution = velocity_distribution  #Has to be a TH1D histogram

    self.speed_of_light = ROOT.RooRealVar("speed_of_light", \
                                          "Speed of Light", 299792.458, "km s^-1")


    self.atomic_mass_of_target = ROOT.RooRealVar("atomic_mass_of_target", \
                                                 "Atomic Mass of Target", 72.6, "amu")  #Atomic mass A
    self.mass_of_target = ROOT.RooFormulaVar("mass_of_target", "Mass of Target", "0.932*@0", \
                                             ROOT.RooArgList(self.atomic_mass_of_target))  #Mass of target M_T
    self.mass_of_target.setUnit("GeV c^-2")

    self.mass_of_detector = ROOT.RooRealVar("mass_of _detector", "Mass of detector", 1, "kg")


    self.mass_of_wimp = ROOT.RooRealVar("mass_of_wimp", "Mass of Wimp", mass_of_wimp, "GeV c^-2")

    self.wimp_nucleon_xs = ROOT.RooRealVar("WIMP-nucleon xs", \
                                           "Cross section of WIMP-nucleon scattering", 1, "pb")
    self.density_of_dark_matter = ROOT.RooRealVar("density_of_dark_matter", \
                                                  "Density of Dark Matter", 0.3, "Gev cm^-3")


    #Reduced masses of WIMP/target nucleus an WIMP/nucleon
    self.mu_wimp_target = ROOT.RooFormulaVar("mu_wimp_target", \
                                           "Reduced mass of WIMP and target nucleus", \
                                           "(@0*@1)/(@0+@1)", \
                                           ROOT.RooArgList(self.mass_of_wimp, self.mass_of_target))
    self.mu_wimp_target.setUnit( self.mass_of_wimp.getUnit() )
    self.mu_wimp_nucleon = ROOT.RooFormulaVar("mu_wimp_nucleon", \
                                           "Reduced mass of WIMP and nucleon", \
                                           "(@0*@1)/(@0*@2+@1)", \
                                           ROOT.RooArgList(self.mass_of_wimp, self.mass_of_target, \
                                                           self.atomic_mass_of_target))
    self.mu_wimp_nucleon.setUnit( self.mass_of_wimp.getUnit() )


    self.constant_factors = ROOT.RooFormulaVar("constant_factors", "Constant factors", \
                                                "((@0*@1*@2)/(2*TMath::Power(@3,2.)*@4*TMath::Power(@5,2.)))*TMath::Power(@6,2.)", \
                                                 ROOT.RooArgList(self.mass_of_detector, self.wimp_nucleon_xs, \
                                                                self.density_of_dark_matter, self.mu_wimp_nucleon, \
                                                                self.mass_of_wimp, self.speed_of_light, self.atomic_mass_of_target))
    self.constant_factors.setUnit("5.61E-16 m c^2 GeV^-1")


    #Base velocity of velocity-distribution
    self.v_sub_0 = ROOT.RooRealVar("v_sub_0", \
                                   "Base Velocity", 220, "km s^-1")
    # Highest velocity possible without leaving galaxy
    self.v_sub_esc = ROOT.RooRealVar("v_sub_esc", \
                                     "Escape Velocity", 544, "km s^-1")
    # Velocity corresponding to minimum energy which can give recoil energy E_R(recoil energy dependent)
    self.v_sub_min = ROOT.RooFormulaVar("v_sub_min", "Minimum Velocity of Minimum Energy", \
                                        "TMath::Sqrt((@0*@1)/(2*TMath::Power(@2,2.)))*TMath::Power(10.,-3.)*@3", \
                                        ROOT.RooArgList(basevars.get_energy(), self.mass_of_target, \
                                                        self.mu_wimp_target, self.speed_of_light))
    self.v_sub_min.setUnit( self.speed_of_light.getUnit() )


    # x-Component of earth velocity
    self.v_sub_E_sub_x = ROOT.RooRealVar("v_sub_E_sub_x", \
                                         "x-Component of earth velocity", 11.1, "km s^-1")

    # y-Component of earth velocity
    self.v_sub_E_sub_y = ROOT.RooRealVar("v_sub_E_sub_y", \
                                         "y-Component of earth velocity", 232.2, "km s^-1")

    # z-Component of earth velocity
    self.v_sub_E_sub_z = ROOT.RooRealVar("v_sub_E_sub_z", \
                                         "z-Component of earth velocity", 7.3, "km s^-1")

    #Modulus of earth velocity
    self.v_sub_E = ROOT.RooFormulaVar("v_sub_E", "Modulus of earth velocity",\
                                      "TMath::Sqrt(@0**2+@1**2+@2**2)",\
                                      ROOT.RooArgList(self.v_sub_E_sub_x, self.v_sub_E_sub_y, \
                                                      self.v_sub_E_sub_z))
    self.v_sub_E.setUnit(self.v_sub_E_sub_x.getUnit())


    # Following is for the Form Factor

    # Momentum transfer as calculated in Lewin/Smith sec. 4
    self.q = ROOT.RooFormulaVar("q", "Momentum Transfer",\
                                "sqrt(2*@0*@1)/197.3", ROOT.RooArgList(\
                                basevars.get_energy(), self.mass_of_target))  #197.3=hbar
    self.q.setUnit("fm^-1")
    # Effective nuclear radius as calculated in Lewin/Smith sec. 4
    self.r_sub_n = ROOT.RooFormulaVar("r_sub_n", "Effective Nuclear Radius",\
                                      "1.14*TMath::Power(@0, 1./3.)",\
                                      ROOT.RooArgList(self.atomic_mass_of_target))
    self.r_sub_n.setUnit("fm")

    self.s = ROOT.RooRealVar("s", "Nuclear Skin Thickness",0.9)
    self.s.setUnit("fm")


    # Woods-Saxon/Helm
    # This is the form-factor we use.
    self.helm_ff_squared = pdfs.MGMWimpHelmFFSquared("helm_ff_squared",\
                                                                 "Helm FF^{2} ",\
                                                                 self.q, self.r_sub_n, self.s)


    # Defining different WIMP models

    self.final_function = pdfs.JKWimpDiffRateEarthPdf("WIMPPDF_With_Escape_Vel", \
                                                              "WIMP Pdf", \
                                                              self.constant_factors, self.v_sub_E_sub_x, \
                                                              self.v_sub_E_sub_y, self.v_sub_E_sub_z, self.v_sub_min, \
                                                              self.v_sub_esc, self.velocity_distribution, \
                                                              self.helm_ff_squared)

    self.final_function_no_ff = pdfs.JKWimpDiffRateEarthPdf("WIMPPDF_With_Escape_Vel_No_ff", \
                                                              "WIMP Pdf no ff", \
                                                              self.constant_factors, self.v_sub_E_sub_x, \
                                                              self.v_sub_E_sub_y, self.v_sub_E_sub_z, self.v_sub_min, \
                                                              self.v_sub_esc, self.velocity_distribution)

    self.final_function_galaxy = pdfs.JKWimpDiffRateGalaxyPdf("WIMPPDF_With_Escape_Vel_Galaxy", \
                                                              "WIMP Pdf galactic system", \
                                                              self.constant_factors, self.v_sub_min, self.v_sub_esc, \
                                                              self.velocity_distribution, self.helm_ff_squared)

    self.final_function_galaxy_no_ff = pdfs.JKWimpDiffRateGalaxyPdf("WIMPPDF_With_Escape_Vel_Galaxy_No_ff", \
                                                              "WIMP Pdf galactic system no ff", \
                                                              self.constant_factors, self.v_sub_min, self.v_sub_esc, \
                                                              self.velocity_distribution)

    self.final_function_ls = pdfs.JKWimpDiffRateLSPdf("WIMPPDF_With_Escape_Vel_LS", \
                                                              "WIMP Pdf Lewin Smith", \
                                                              self.constant_factors, self.v_sub_E, self.v_sub_min, \
                                                              self.v_sub_esc, self.v_sub_0, self.helm_ff_squared)

    self.final_function_ls_no_ff = pdfs.JKWimpDiffRateLSPdf("WIMPPDF_With_Escape_Vel_LS_No_ff", \
                                                              "WIMP Pdf Lewin Smith no ff", \
                                                              self.constant_factors, self.v_sub_E, self.v_sub_min, \
                                                              self.v_sub_esc, self.v_sub_0)

    self.final_function_sfg = pdfs.JKWimpDiffRateSFGPdf("WIMPPDF_With_Escape_Vel_SFG", \
                                                              "WIMP Pdf Savage Freese Gondolo", \
                                                              self.constant_factors, self.v_sub_E, self.v_sub_min, \
                                                              self.v_sub_esc, self.v_sub_0, self.helm_ff_squared)

    self.final_function_sfg_no_ff = pdfs.JKWimpDiffRateSFGPdf("WIMPPDF_With_Escape_Vel_SFG_No_ff", \
                                                              "WIMP Pdf Savage Freese Gondolo no ff", \
                                                              self.constant_factors, self.v_sub_E, self.v_sub_min, \
                                                              self.v_sub_esc, self.v_sub_0)


  def set_v_sub_E(self, v_sub_E):
    scaling = v_sub_E/self.v_sub_E.getVal()
    self.v_sub_E_sub_x.setVal(scaling*self.v_sub_E_sub_x.getVal())
    self.v_sub_E_sub_y.setVal(scaling*self.v_sub_E_sub_y.getVal())
    self.v_sub_E_sub_z.setVal(scaling*self.v_sub_E_sub_z.getVal())


  # Get Woods-Saxon/Helm form factor ( returns RooAbsPdf )
  def get_helm_form_factor(self):
    return self.helm_ff_squared


  # Get different 2-D WIMP models . The following allow
  # changing the value of the WIMP mass by passing
  # in an ( optional ) float parameter

  # WIMP Model, earth rest system, Helm form factor
  def get_WIMP_model(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function

  # WIMP Model, earth rest system, no Helm form factor
  def get_WIMP_model_no_ff(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_no_ff

  # WIMP Model, galactic rest system, Helm form factor
  def get_WIMP_model_galaxy(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_galaxy

  # WIMP Model, galactic rest system, no Helm form factor
  def get_WIMP_model_galaxy_no_ff(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_galaxy_no_ff

  # WIMP Model, Lewin Smith, Helm form factor
  def get_WIMP_model_ls(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_ls

  # WIMP Model, Lewin Smith, no Helm form factor
  def get_WIMP_model_ls_no_ff(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_ls_no_ff

  # WIMP Model, Savage Freese Gondolo, Helm form factor
  def get_WIMP_model_sfg(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_sfg

  # WIMP Model, Savage Freese Gondolo, no Helm form factor
  def get_WIMP_model_sfg_no_ff(self, mass_of_wimp=None):
    if mass_of_wimp: self.mass_of_wimp.setVal(mass_of_wimp)
    return self.final_function_sfg_no_ff
