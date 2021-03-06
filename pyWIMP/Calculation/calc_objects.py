#
# calc_object exports classes which can be 
# used for calculation of sensitivities
#
import ROOT
import ExclusionCalculation as ec
import OscillationSensitivityCalculation as osc
import DataCalculation as dat
from  ..utilities.utilities import unroll_RooAbsPdf
from pyWIMP.DMModels.gaussian_signal import GaussianSignalModel 
from pyWIMP.DMModels.tritium_decay_model import TritiumDecayModel 

class WIMPModel:
    """
    Class handles performing a sensitivity calculation
    for a Ge detector given the input parameters.
    FixME: This class uses AllWIMPModels class, but doesn't
    allow an interface to change those values.
    """
    def get_requested_values(cls):
        """
        Returns requested values plus defaults
        """
        return {'total_time' : ('Total time (year)', 5),
                'threshold' : ('Threshold (keV)', 1),
                'energy_max' : ('Maximum energy (keV)', 50),
                'num_time_bins' : ('Number of time bins, 0 means unbinned', 0),
                'num_energy_bins' : ('Number of energy bins, 0 means unbinned', 0),
                'mass_of_detector' : ('Mass of detector (kg)', 1),
                'background_rate' : ('Background rate (counts/keV/kg/day)', 0.1),
                'tritium_activation_rate' : ('Tritium activation rate (counts/kg/day)', 200),
                'tritium_exposure_time' : ('Tritium exposure time days', 0),
                'wimp_mass' : ('WIMP mass (GeV/c^-2)', 10),
                'confidence_level' : ('Confidence level (0 -> 1)', 0.9),
                'variable_quenching' : ('Set to use variable quenching', False),
                'constant_time' : ('Set time as constant', False),
                'constant_energy' : ('Set energy as constant', False),
                'show_plots' : ('Show plot results of fit', False),
                'do_axioelectric' : ('Do axioelectric fitting', False),
                'axion_mass' : ('Axion mass in keV [default 0]', 0.),
                'print_out_plots' : ("""Print plot results of fit to eps format.  
This will cause the program to continue
instead of stopping once a plot is made.
Batch mode will be set so that no plot
will be displayed during the program.
                                     """, False),
                'debug' : ('Set debug flag, enables verbose output', False) }
    get_requested_values = classmethod(get_requested_values)

    def __init__(self, 
                 output_pipe,
                 exit_manager,
                 num_iterations,
                 input_variables):

        self.input_variables = input_variables
        for akey, val in self.get_requested_values().items():
            aval = val[1]
            if akey in input_variables.keys():
                aval = input_variables[akey]
            setattr(self,akey, aval) 
        self.number_iterations = num_iterations
        self.output_pipe = output_pipe
        self.exit_now = False
        self.is_initialized = False
        self.exit_manager = exit_manager 
        if self.debug:
            self.print_level = 3
            self.verbose = True
        else:
            self.print_level = -1
            self.verbose = False

        self.do_bin_data = False
        if self.num_energy_bins != 0 and self.num_time_bins != 0: self.do_bin_data = True

    # overload this function for derived classes.
    def initialize(self):

        if self.is_initialized: return
        from pyWIMP.DMModels.wimp_model import WIMPModel
        from pyWIMP.DMModels.base_model import BaseVariables
        from pyWIMP.DMModels.flat_model import FlatModel

        self.total_counts = int(self.mass_of_detector*
                                self.background_rate*
                                (self.energy_max-self.threshold)*
                                self.total_time*365)
        self.basevars = BaseVariables(time_beginning=0,
            time_in_years=self.total_time,
            energy_threshold=self.threshold,
            energy_max=self.energy_max)

        if self.num_energy_bins != 0: self.basevars.get_energy().setBins(int(self.num_energy_bins))
        if self.num_time_bins != 0: self.basevars.get_time().setBins(int(self.num_time_bins))
        self.variables = ROOT.RooArgSet()
        if self.constant_time:
            self.basevars.get_time().setVal(0)
            self.basevars.get_time().setConstant(True)
        else:
            self.variables.add(self.basevars.get_time())
        if not self.constant_energy:
            self.variables.add(self.basevars.get_energy())


        self.calculation_class = \
            ec.ExclusionCalculation(self.exit_manager)
 
        # This is where we define our models
        self.backgroundClass = TritiumDecayModel(self.basevars, 
                                                 self.tritium_exposure_time, 
                                                 self.tritium_activation_rate,
                                                 self.mass_of_detector,
                                                 self.background_rate)

        # This model is already an extended model
        self.background_model =  self.backgroundClass.get_model()

        self.background_extend = self.background_model
        if not self.do_axioelectric:
            # The following has not been normalized to per-nucleon yet.
            self.model_normal = ROOT.RooRealVar("model_normal", 
                                                "WIMP-nucleus #sigma", 
                                                0, -1e-15, 0.1, 'pb')
            self.wimpClass = WIMPModel(self.basevars,
                mass_of_wimp=self.wimp_mass,
                kilograms = self.mass_of_detector,
                constant_quenching=(not self.variable_quenching))

            self.model = self.wimpClass.get_model()
            self.norm = self.wimpClass.get_normalization().getVal()
            self.model_extend = ROOT.RooExtendPdf("model_extend", 
                                                   "model_extend", 
                                                   self.model, 
                                                   self.model_normal)
        else:
            # wimpClass is of course a misnomer here, but we use it for now
            self.wimpClass = GaussianSignalModel(self.basevars,
                                                 mean_of_signal=self.axion_mass)
 
            self.model_normal = ROOT.RooRealVar("model_normal", 
                                                "Counts", 
                                                0, -10, 1000)
            # This is where we define our model
            self.model = self.wimpClass.get_model()
            self.norm = self.wimpClass.get_normalization()*self.model.getNorm(
                        ROOT.RooArgSet(self.basevars.get_energy())) 

            # We actually want the inverse of the normaliation, since this is later multiplied 
            # and we want the total counts in a particular gaussian.
            self.norm = 1./self.norm
            self.model_extend = ROOT.RooExtendPdf("model_extend", 
                                                   "model_extend", 
                                                   self.model, 
                                                   self.model_normal)


        self.added_pdf = ROOT.RooAddPdf("b+s", 
                                        "Background + Signal", 
                                        ROOT.RooArgList(
                                        self.background_extend, 
                                        self.model_extend))
 
        self.test_variable = self.model_normal
        self.data_set_model = self.background_extend
        self.fitting_model = self.added_pdf
        self.is_initialized = True
    
    def run(self):
        """
        Do the work.  Perform the fits, and return the results
        This function runs from the base class, so derived classes should
        not need to overload it.
        """
        import pickle
        import signal
        import os

        ROOT.RooRandom.randomGenerator().SetSeed(0)
        self.initialize()
        ROOT.gROOT.SetBatch()
        if self.show_plots or self.print_out_plots:
            self.calculation_class.set_canvas(ROOT.TCanvas())
        if self.show_plots and not self.print_out_plots:
            ROOT.gROOT.SetBatch(0)

        if not self.debug:
            ROOT.RooMsgService.instance().setSilentMode(True)
            ROOT.RooMsgService.instance().setGlobalKillBelow(4)

        # Open the pipe to write back on
        write_pipe = os.fdopen(self.output_pipe, 'w') 

        # Set the numeric integration properties, this is important
        # since some of these integrals are difficult to do
        precision = ROOT.RooNumIntConfig.defaultConfig()
        precision.setEpsRel(1e-8)
        precision.setEpsAbs(1e-8)
        precision.method2D().setLabel("RooIntegrator2D")

        # Perform the integral.  This is a bit of sanity check, it this fails, 
        # then we have a problem 
        norm_integral = self.model.createIntegral(self.variables)

        # FixME, this integral sometimes doesn't converge, set a timeout?
        # This integral is in units of pb^{-1} 
        norm_integral_val = norm_integral.getVal()

        print norm_integral_val
        if norm_integral_val == 0.0:
            print "Integral defined as 0, meaning it is below numerical precision"
            print "Aborting further calculation"
            write_pipe.write(pickle.dumps({}))
            write_pipe.close()
            return
 
        # Set up signal handler
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGUSR1, self.exit_manager.exit_handler)

       

        self.calculation_class.set_debug(self.debug)
        self.calculation_class.set_show_plots(self.show_plots)
        self.calculation_class.set_print_out_plots(self.print_out_plots)
        self.calculation_class.set_input_variables(self.input_variables)
        self.calculation_class.set_plot_base_name(
            "%s WIMP Mass: %g GeV" % (self.__class__.__name__,
                                  self.wimp_mass))
        results_list = \
            self.calculation_class.scan_confidence_value_space_for_model(
                self.fitting_model, 
                self.data_set_model, 
                self.test_variable, 
                self.norm, 
                self.variables, 
                self.do_bin_data, 
                self.number_iterations, 
                self.confidence_level)

        write_pipe.write(pickle.dumps(results_list))
        write_pipe.close()

class DataExclusion(WIMPModel):
    def get_requested_values(cls):
        adict = WIMPModel.get_requested_values()
        del adict['constant_energy']
        del adict['constant_time']
        del adict['background_rate']
        adict['fix_l_line_ratio'] = ('Fix ratio of the Ge and Zn L-lines', False)
        adict['data_file'] = ('Name of data root file', 'temp.root')
        adict['object_name'] = ("""Name of object inside data file. 
This can be a:                                             
                                                                
TH1: the bins will be interprected as energy. 
A RooDataHist will be used to generate a binned fit.  
                                                                
TTree: with branches, 'ee_energy', 'time', 'weight', plus others.
'time' and ee_energy are both optional, it is assumed that if they
don't exist that time and energy are constant.
'weight' is optional, but when present will adjust the weight of each
entry to 'weight' in the RooDataSet.
A RooDataSet will be used to generate an unbinned fit. 
                                """, 'output_data')
        adict['data_set_cuts'] = ("""String of cuts to apply. 
This only applies if the object passed
in is a TTree (see object_name).  These cuts will be used to generate
a subset of the TTree and pass into RooDataSet.
                                  """, '')
        return adict
    get_requested_values = classmethod(get_requested_values)
    def initialize(self):
        # The background model is the same, but now we switch it to the data model
        from pyWIMP.DMModels.low_energy_background import LowEnergyBackgroundModel
        from pyWIMP.DMModels.base_model import BaseVariables
        from pyWIMP.DMModels.wimp_model import WIMPModel
        self.total_counts = None
        open_file = ROOT.TFile(self.data_file)
        self.workspace = open_file.Get(self.object_name)
        # Do some introspection, we can handle TH1s, and RooAbsData 
        self.basevars = BaseVariables(time_beginning=0,
            time_in_years=self.total_time,
            energy_threshold=self.threshold,
            energy_max=self.energy_max,
            use_tag=False)

        self.variables = ROOT.RooArgSet()

        if self.workspace.InheritsFrom(ROOT.TH1.Class()):
            self.variables.add(self.basevars.get_energy())
            self.basevars.get_time().setVal(0)
            self.basevars.get_time().setConstant(True)

            self.data_set_model = ROOT.RooDataHist("data", "data", 
                                    ROOT.RooArgList(self.basevars.get_energy()),
                                    self.workspace)
        elif self.workspace.InheritsFrom(ROOT.TTree.Class()):
            # Default to setting them to constant.
            self.basevars.get_time().setVal(0)
            self.basevars.get_time().setConstant(True)
            self.basevars.get_energy().setVal(0)
            self.basevars.get_energy().setConstant(True)
            if self.num_energy_bins != 0:
                self.basevars.get_energy().setBins(int(self.num_energy_bins))
            if self.num_time_bins != 0:
                self.basevars.get_time().setBins(int(self.num_time_bins))

            branches = self.workspace.GetListOfBranches()
            efficiency = "" 
            branch_arg_list = []

            for i in range(branches.GetEntries()):
                branch_name = branches[i].GetName()
                if branch_name == "ee_energy": 
                    self.basevars.get_energy().setConstant(False)
                    self.variables.add(self.basevars.get_energy())
                    branch_arg_list.append((self.basevars.get_energy(), 
                                            branch_name))
                if branch_name == "time": 
                    self.basevars.get_time().setConstant(False)
                    self.variables.add(self.basevars.get_time())
                    branch_arg_list.append((self.basevars.get_time(), 
                                            branch_name))
                elif branch_name == "weight":
                    efficiency = branch_name 
                    self.variables.add(self.basevars.get_weighting())

            if not self.data_set_cuts:
                # Load the DataSet the easy way
                self.data_set_model = ROOT.RooDataSet("data", "data", 
                                        self.workspace,
                                        self.variables,
                                        efficiency)
            else:
                # Otherwise, we have to get the correct events,
                # which requires stepping through all events
                self.data_set_model = ROOT.RooDataSet("data", "data", 
                                        self.variables,
                                        efficiency)
                # Get an event list with the correct cut events
                ROOT.gROOT.cd()
                el = ROOT.TEventList("el", "el")
                cuts = self.data_set_cuts
                iter = self.variables.createIterator()
                while 1: 
                    obj = iter.Next()
                    if not obj: break
                    if obj.GetTitle() == 'weight': continue
                    cuts += " && ((%s <= %f) && (%s >= %f))" % (obj.GetName(),
                                                                obj.getMax(),
                                                                obj.GetName(),
                                                                obj.getMin())
                self.workspace.Draw(">>%s" % el.GetName(), cuts, "goff")
                event_list = [el.GetEntry(i) for i in range(el.GetN())] 
                #event_list = [i for i in range(self.workspace.GetEntries())] 
                for j in event_list: 
                    self.workspace.GetEntry(j)
                    eff_val = 1.
                    if efficiency: 
                        eff_val = getattr(self.workspace, efficiency)
                    for val in branch_arg_list: 
                        val[0].setVal(getattr(self.workspace, val[1]))
                    if eff_val != 0: self.data_set_model.add(self.variables, eff_val)
                     
        else:
            print "Requested: %s, isn't a TTree or TH1!" % self.data_set_name
            raise TypeError

        if self.num_energy_bins != 0 or self.num_time_bins != 0:
            self.original_data_set = self.data_set_model
            self.data_set_model = self.original_data_set.binnedClone()

        if self.data_set_model.isWeighted(): print "Data set is weighted"
        print "Data set has %i entries." % self.data_set_model.sumEntries()


        if not self.do_axioelectric:
            self.wimpClass = WIMPModel(self.basevars,
                mass_of_wimp=self.wimp_mass,
                kilograms = self.mass_of_detector,
                constant_quenching=(not self.variable_quenching))
 
            # This is where we define our model
            self.model = self.wimpClass.get_model()
            #self.model = self.wimpClass.get_simple_model()
            self.norm = self.wimpClass.get_normalization().getVal()

            # The following has not been normalized to per-nucleon yet.
            self.model_normal = ROOT.RooRealVar("model_normal", 
                                                "WIMP-nucleus #sigma", 
                                                1, -10, 100, 
                                                "pb")
            self.model_extend = ROOT.RooExtendPdf("model_extend", 
                                                   "model_extend", 
                                                   self.model, 
                                                   self.model_normal)
            # Getting the number of events for a model_normal of 1 
            # This gives us number of events per model_normal value
            scaler = self.model_extend.expectedEvents(self.variables)
            self.model_normal.setMax(2*self.data_set_model.sumEntries()/scaler)
        else:
            # wimpClass is of course a misnomer here, but we use it for now
            self.wimpClass = GaussianSignalModel(self.basevars,
                                                 mean_of_signal=self.axion_mass)
 
            self.model_normal = ROOT.RooRealVar("model_normal", 
                                                "Counts", 
                                                0, -10, 1000)
            # This is where we define our model
            self.model = self.wimpClass.get_model()
            self.norm = self.wimpClass.get_normalization()*self.model.getNorm(
                        ROOT.RooArgSet(self.basevars.get_energy())) 

            # We actually want the inverse of the normaliation, since this is later multiplied 
            # and we want the total counts in a particular gaussian.
            self.norm = 1./self.norm
            self.model_extend = ROOT.RooExtendPdf("model_extend", 
                                                   "model_extend", 
                                                   self.model, 
                                                   self.model_normal)

        self.is_initialized = True

        self.calculation_class = \
            dat.DataCalculation(self.exit_manager)
        self.low_energy = LowEnergyBackgroundModel(self.basevars, 
                                                   use_ratio=self.fix_l_line_ratio)

        list_of_models, list_of_coefficients = self.low_energy.get_list_components()
        self.extended_models = []
        i = 0
        while 1: 
            amod = list_of_models.at(i)
            avar = list_of_coefficients.at(i)
            if not amod: break
            i += 1
            extend = ROOT.RooExtendPdf("extend%s" % amod.GetName(),
                                       "extend%s" % amod.GetName(),
                                       amod, avar)
            self.extended_models.append(extend)
        temp_list = ROOT.RooArgList()
        temp_list.add(self.model_extend)
        for amod in self.extended_models:
            temp_list.add(amod)
        self.added_pdf = ROOT.RooAddPdf("b+s", 
                                        "Background + Signal", 
                                        temp_list)
        self.test_variable = self.model_normal
        self.fitting_model = self.added_pdf
        


class OscillationSignalDetection(WIMPModel):
    def get_requested_values(cls):
        adict = WIMPModel.get_requested_values()
        del adict['constant_energy']
        del adict['constant_time']
        del adict['wimp_mass']
        del adict['variable_quenching']
        adict['model_amplitude'] = ('Initial model amplitude', 0.1)
        return adict
    get_requested_values = classmethod(get_requested_values)

    # overload this function for derived classes.
    def initialize(self):

        if self.is_initialized: return
        from pyWIMP.DMModels.base_model import BaseVariables
        from pyWIMP.DMModels.flat_model import FlatModel
        from pyWIMP.DMModels.oscillation_model import OscillationModel
        self.basevars = BaseVariables(time_beginning=0,
            time_in_years=self.total_time,
            energy_threshold=self.threshold,
            energy_max=self.energy_max)

        self.oscClass = OscillationModel(self.basevars)

        self.flatClass = FlatModel(self.basevars)


        self.calculation_class = \
            osc.OscillationSensitivityCalculation(self.exit_manager)

        self.variables = ROOT.RooArgSet()
        self.variables.add(self.basevars.get_time())

        self.basevars.get_energy().setVal(0)
        self.basevars.get_energy().setConstant(True)

        # This is where we define our models
        self.background =  self.flatClass.get_model()
        self.model = self.oscClass.get_model()
        self.norm = 1 
        self.is_initialized = True

        
        if self.model_amplitude > 1: self.model_amplitude = 1
        elif self.model_amplitude < 0: self.model_amplitude = 0
        self.signal_percentage = ROOT.RooRealVar("signal_percentage", 
                                            "signal_percentage", 
                                            self.model_amplitude)
        self.background_model = ROOT.RooAddPdf(
                                "background",
                                "Data Model",
                                self.model,
                                self.background,
                                self.signal_percentage)

        self.model_normal = ROOT.RooRealVar("model_normal", 
                                            "model_normal", 
                                            self.model_amplitude, 
                                            0, 1)

        self.total_fit_counts = ROOT.RooRealVar("total_fit_counts", 
                                            "total_fit_counts", 
                                            self.total_counts, 
                                            0, 3*self.total_counts)
        self.added_pdf = ROOT.RooAddPdf(
                                "added_pdf",
                                "Fit Model",
                                self.model,
                                self.background,
                                self.model_normal)

        self.model_extend = ROOT.RooExtendPdf("model_extend",
                                              "Signal + Background",
                                              self.added_pdf,
                                              self.total_fit_counts)

        self.test_variable = self.model_normal
        self.data_set_model = self.background_model
        self.fitting_model = self.model_extend 


