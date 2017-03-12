import ROOT


def get_efficiency(Energy, E_thresh=1.58, sigma_rec=0.51):
  efficiency = ROOT.TF1("efficiency", "0.5*(1+TMath::Erf((x-[0])/([1]*TMath::Sqrt(2))))", Energy-1., Energy+1.)
  efficiency.SetParameter(0, E_thresh)
  efficiency.SetParameter(1, sigma_rec)
  return efficiency.Eval(Energy)

def get_hist_with_efficiency(hist_in, hist_out_name, E_thresh=1.58, sigma_rec=0.51):
  efficiency = ROOT.TF1("efficiency", "0.5*(1+TMath::Erf((x-[0])/([1]*TMath::Sqrt(2))))", hist_in.GetXaxis().GetXmin(), hist_in.GetXaxis().GetXmax())
  efficiency.SetParameter(0, E_thresh)
  efficiency.SetParameter(1, sigma_rec)
  hist_out = hist_in.Clone(hist_out_name)
  hist_out.Multiply(efficiency)
  return hist_out

def get_efficiency_integral(hist_in, E_thresh=1.58, sigma_rec=0.51):
  hist_temp = get_hist_with_efficiency(hist_in, "temp", E_thresh, sigma_rec)
  return hist_temp.Integral("width")

def get_efficiency_curve(Nbins, x_min, x_max, hist_out_name="efficiency_curve", E_thresh=1.58, sigma_rec=0.51):
  hist_out = ROOT.TH1F(hist_out_name, hist_out_name, Nbins, x_min, x_max)
  for i in range(0,Nbins):
    hist_out.SetBinContent(i, get_efficiency(hist_out.GetBinCenter(i)))
  return hist_out
