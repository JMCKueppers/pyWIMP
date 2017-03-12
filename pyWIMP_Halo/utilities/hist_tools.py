import ROOT
import numpy as np
import os

def hist_from_data(file_path, hist_bins, delimiter=",", correct_hist_gaps = True, normalize = True):
  # Read data from file
  data_array = np.loadtxt(file_path, delimiter=",")
  
  # Define histogram for velocity distribution
  bin_width = (data_array[-1][0]-data_array[0][0])/float(hist_bins)
  hist_min = data_array[0][0]
  hist_max = data_array[-1][0]-bin_width
  hist = ROOT.TH1D(os.path.basename(file_path), "velocity distribution", hist_bins, hist_min, hist_max)

  # Set content of bins
  print "Filling histogram\n"
  for i in range(0, hist_bins):
    bin_min = i*bin_width
    bin_max = (i+1)*bin_width
    bin_content = 0
    bin_data_number = 0
  
    # Fill bin with mean value of all data within bin limits
    for j in range(0, len(data_array)):
      if (data_array[j][0] >= bin_min) and (data_array[j][0] < bin_max):
        bin_content += data_array[j][1]
        bin_data_number += 1
      elif data_array[j][0] >= bin_max:
        break  # Fasten up loop
  
    # For physical reasons, no bin content is below zero
    if bin_content > 0:
      bin_content /= bin_data_number  # Calculate mean value for bin
      hist.SetBinContent(i, bin_content)
    else:
      hist.SetBinContent(i, 0)
  
  # Fill gaps in histogram
  if correct_hist_gaps:
    print "Filling gaps\n"
    upper_limit = hist_bins-1
    lower_limit = 0
    while hist.GetBinContent(upper_limit) <= 0:
      upper_limit -= 1
    while hist.GetBinContent(lower_limit) <= 0:
      lower_limit += 1
  
    # Approximate gap content with linear interpolation
    for i in range(lower_limit, upper_limit):
      if hist.GetBinContent(i) <= 0:
        width = 0
        while hist.GetBinContent(i+width) <= 0:
          width += 1
        bin_content = hist.GetBinContent(i-1)+((hist.GetBinContent(i+width)-hist.GetBinContent(i-1))/width)
        hist.SetBinContent(i, bin_content)

  if normalize:
    print "Normalizing input histogram\n"
    normalization = hist.Integral("width")
    hist.Scale(1./normalization)
  return hist

def get_hist_ratio (h0, h1, hist_out_name, percent_axis=True):
  xmin = np.array([h0.GetXaxis().GetXmin(), h1.GetXaxis().GetXmin()])
  xmax = np.array([h0.GetXaxis().GetXmax(), h1.GetXaxis().GetXmax()])
  nbins = np.array([h0.GetNbinsX(), h1.GetNbinsX()])
  if xmin[0]==xmin[1] and xmax[0]==xmax[1] and nbins[0]==nbins[1]:
    h_temp = [h0.Clone(hist_out_name), h1.Clone("copy1"+hist_out_name)]
  else:
    h_in = [h0.Clone("copy0"+hist_out_name), h1.Clone("copy1"+hist_out_name)]
    h_temp = [ROOT.TH1D("temp0"+hist_out_name, "temp 0", nbins.min(), xmin.min(), xmax.max()), ROOT.TH1D("temp1"+hist_out_name, "temp 1", nbins.min(), xmin.min(), xmax.max())]
    bin_width = (xmax.max()-xmin.min())/nbins.min()
    for m in [0,1]:
      for i in range(0, nbins.min()):
        lower_edge = i*bin_width
        upper_edge = (i+1)*bin_width
        if upper_edge < xmin[m] or lower_edge > xmax[m]:
          h_temp[m].SetBinContent(i,0)
        else:
          bin_content = 0
          bin_amount = 0
          for j in range(0, nbins[m]):
            if h_in[m].GetBinLowEdge(j) >= lower_edge and h_in[m].GetBinLowEdge(j) <= upper_edge:
              bin_content += h_in[m].GetBinContent(j)
              bin_amount += 1
          if bin_amount > 1:
            bin_content /= bin_amount
          h_temp[m].SetBinContent(i, bin_content)
  hist_out = h_temp[0].Clone(hist_out_name)
  hist_out.SetTitle(hist_out_name)
  hist_out.Divide(h_temp[1])
  if percent_axis:
    for i in range(hist_out.GetNbinsX()):
      hist_out.SetBinContent(i, hist_out.GetBinContent(i)*100.)
  return hist_out

def get_hist_deviation (h0, h1, hist_out_name, percent_axis=True):
  hist_out = get_hist_ratio (h0, h1, hist_out_name, percent_axis=False)
  for i in range(hist_out.GetNbinsX()):
    hist_out.SetBinContent(i, hist_out.GetBinContent(i)-1.)
  if percent_axis:
    for i in range(hist_out.GetNbinsX()):
      hist_out.SetBinContent(i, hist_out.GetBinContent(i)*100.)
  return hist_out

def draw_hist(hist_in, line_color = 1, line_style = 1):
  hist_in.SetLineColor(line_color)
  hist_in.SetLineStyle(line_style)
  hist_in.Draw("Same")
  
def export_tf1(function, file_path, bins):
  func_min = function.GetXmin()
  func_max = function.GetXmax()
  bin_width = (func_max-func_min)/bins
  print "Export function\n"
  output_file = open(file_path, "w")
  for i in range(0, bins):
    output_file.write(str(i*bin_width)+", "+str(function.Eval(i*bin_width))+"\n")
  output_file.close()


def draw_spectra(hists1, colors, output_path, hists2=None, x_title="Recoil energy (keV)", y_title=None, x_min=0., x_max=200., y_min=None, y_max=None, legend=None, text=None, logscale=0, percent_axis=False):
  c1 = ROOT.TCanvas()
  c1.SetGrid(0,0)
  c1.SetFrameLineColor(ROOT.kGray+2)
  if logscale == 1:
    c1.SetLogy()
  if logscale == 2:
    c1.SetLogy()
    c1.SetLogx()
  c1.SetLeftMargin(0.15)
  c1.SetRightMargin(0.11)

  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)
  
  for i in range(len(hists1)):
    #hists1[i].GetXaxis().SetLimits(x_min,x_max)
    hists1[i].GetXaxis().SetRangeUser(x_min,x_max)
    if y_min and y_max:
      hists1[i].GetYaxis().SetRangeUser(y_min,y_max)
      #hists1[i].GetYaxis().SetLimits(y_min,y_max)
    hists1[i].GetXaxis().SetTitle(x_title)
    if y_title:
      hists1[i].GetYaxis().SetTitle(y_title)
    else:
      hists1[i].GetYaxis().SetTitle("Events / "+str((hists1[i].GetXaxis().GetXmax()-hists1[i].GetXaxis().GetXmin())/hists1[i].GetNbinsX())+" keV")
    hists1[i].GetXaxis().SetTitleSize(0.03*1./0.7)
    hists1[i].GetYaxis().SetTitleSize(0.03*1./0.7)
    hists1[i].GetXaxis().SetTitleOffset(1.1)
    hists1[i].GetYaxis().SetTitleOffset(1.8)
    hists1[i].GetYaxis().SetLabelSize(0.03*1./0.7)
    hists1[i].GetXaxis().SetLabelSize(0.03*1./0.7)
    hists1[i].GetXaxis().SetAxisColor(ROOT.kGray+2)
    hists1[i].GetYaxis().SetAxisColor(ROOT.kGray+2)
    hists1[i].GetXaxis().SetTickLength(0.01*1.0/0.7)
    hists1[i].GetYaxis().SetTickLength(0.0125)
    hists1[i].SetLineWidth(2)
    draw_hist(hists1[i], colors[i], 1)
  if hists2:
    for i in range(len(hists2)):
      if y_min and y_max:
        hists2[i].GetYaxis().SetRangeUser(y_min,y_max)
      hists2[i].SetLineWidth(2)
      draw_hist(hists2[i], colors[i], 2)

  if text:
    text_field = ROOT.TPaveText(0.15,0.9,0.95,0.95, "NDC")
    text_field.AddText(text)
    text_field.SetFillColorAlpha(0,0)
    text_field.SetTextAlign(12)
    text_field.SetBorderSize(0)
    text_field.SetTextSize(0.03*1./0.7)
    text_field.Draw()

  if legend:
    legend_field = ROOT.TLegend(0.57,0.6,0.95,0.89)
    legend_field.SetFillColorAlpha(0,0)
    legend_field.SetTextAlign(12)
    legend_field.SetBorderSize(0)
    legend_field.SetTextSize(0.03*1./0.7)
    for i in range(len(legend)):
      hists1[i].SetLineColor(colors[i])
      legend_field.AddEntry(hists1[i], legend[i], "L")
    legend_field.Draw()

  if percent_axis:
    rightmax = 105.
    axis = ROOT.TGaxis(3.5, ROOT.gPad.GetUymin(), 3.5, 1.05, 0, rightmax, 510, "+LS")
    axis.SetLineColor(ROOT.kGray+2)
    axis.SetTitle("Deviation (%)")
    axis.SetTitleSize(0.03*1./0.7)
    axis.SetLabelSize(0.03*1./0.7)
    axis.SetTitleFont(42)
    axis.SetLabelFont(42)
    axis.SetTickLength(0.0125)
    axis.SetTitleOffset(1.3)
    axis.SetLabelOffset(0.02)
    axis.Draw()

  c1.Print(os.path.splitext(output_path)[0]+".pdf")
  print "Created file "+os.path.splitext(output_path)[0]+".pdf\n"
  c1.Print(os.path.splitext(output_path)[0]+".root")
  print "Created file "+os.path.splitext(output_path)[0]+".root\n"
  c1.Print(os.path.splitext(output_path)[0]+".xml")
  print "Created file "+os.path.splitext(output_path)[0]+".xml\n"



def draw_spectra_with_ratio(hists1, hists2, ratio_hists, colors, output_path, x_title="Recoil energy (keV)", y_title=None, x_min=0., x_max=200., y_min=None, y_max=None, ry_min=0.8, ry_max=1.2, percent_band=5., legend=None, text=None, logscale=0):
  c1 = ROOT.TCanvas()
  c1.Divide(1,2)
  if logscale == 1:
    c1.GetPad(1).SetLogy()
  if logscale == 2:
    c1.GetPad(1).SetLogy()
    #c1.GetPad(1).SetLogx()  #Logscale x-axis not uniform in plot and ratio-plot
    #c1.GetPad(2).SetLogx()
  c1.GetPad(1).SetPad(0.,0.3,1.,1.)
  c1.GetPad(2).SetPad(0.,0.,1.,0.3)
  c1.GetPad(1).SetBottomMargin(0.02)
  c1.GetPad(2).SetTopMargin(0.05)
  c1.GetPad(2).SetBottomMargin(0.4)
  c1.GetPad(1).SetGrid(0,0)
  c1.GetPad(2).SetGrid(0,1)
  c1.GetPad(1).SetFrameLineColor(ROOT.kGray+2)
  c1.GetPad(2).SetFrameLineColor(ROOT.kGray+2)
  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetOptStat(0)
  
  #Draw 5% band
  c1.cd(2)
  x=[x_min, x_max]
  y=[0., 0.]
  xe=[0., 0.]
  ye=[percent_band, percent_band]
  Area = ROOT.TGraphErrors(2,np.array(x),np.array(y),np.array(xe),np.array(ye))
  #Area.GetXaxis().SetLimits(x_min,x_max)
  Area.GetXaxis().SetRangeUser(x_min,x_max)
  Area.GetYaxis().SetRangeUser(-100.*(1.-ry_min),100.*(ry_max-1.))
  Area.GetYaxis().SetTitle("Deviation (%)")
  Area.GetXaxis().SetTitle(x_title)
  Area.GetXaxis().SetTitleSize(0.03*1./0.3)
  Area.GetYaxis().SetTitleSize(0.03*1./0.3)
  Area.GetXaxis().SetTitleOffset(1.3)
  Area.GetYaxis().SetTitleOffset(0.5)
  Area.GetXaxis().SetLabelOffset(0.03)
  Area.GetXaxis().SetLabelSize(0.03*1./0.3)
  Area.GetYaxis().SetLabelSize(0.03*1./0.3)
  Area.GetXaxis().SetAxisColor(ROOT.kGray+2)
  Area.GetYaxis().SetAxisColor(ROOT.kGray+2)
  Area.GetYaxis().SetNdivisions(505)
  Area.GetXaxis().SetTickLength(0.01*1.0/0.3)
  Area.GetYaxis().SetTickLength(0.02)
  Area.SetFillColor(ROOT.kYellow-8)
  Area.Draw("E3AL")
  Line = ROOT.TLine()
  Line.SetLineColor(ROOT.kGray+1)
  Line.SetLineStyle(2)
  Line.DrawLine(x_min,percent_band, x_max,percent_band)
  Line.DrawLine(x_min,-percent_band, x_max,-percent_band)
  
  for i in range(len(hists1)):
    c1.cd(1)
    hists1[i].GetXaxis().SetRangeUser(x_min,x_max)
    #hists1[i].GetXaxis().SetLimits(x_min,x_max)
    if y_min and y_max:
      hists1[i].GetYaxis().SetRangeUser(y_min,y_max)
    hists1[i].GetXaxis().SetTitle("")
    if y_title:
      hists1[i].GetYaxis().SetTitle(y_title)
    else:
      hists1[i].GetYaxis().SetTitle("Events / "+str((hists1[i].GetXaxis().GetXmax()-hists1[i].GetXaxis().GetXmin())/hists1[i].GetNbinsX())+" keV")
    hists1[i].GetXaxis().SetTitleSize(0.03*1./0.7)
    hists1[i].GetYaxis().SetTitleSize(0.03*1./0.7)
    hists1[i].GetYaxis().SetTitleOffset(1.15)
    hists1[i].GetXaxis().SetLabelSize(0.)
    hists1[i].GetYaxis().SetLabelSize(0.03*1./0.7)
    hists1[i].GetXaxis().SetAxisColor(ROOT.kGray+2)
    hists1[i].GetYaxis().SetAxisColor(ROOT.kGray+2)
    hists1[i].GetXaxis().SetTickLength(0.01*1.0/0.7)
    hists1[i].GetYaxis().SetTickLength(0.0125)
    hists1[i].SetLineWidth(2)
    draw_hist(hists1[i], colors[i], 1)
    c1.cd(2)
    ratio_hists[i].GetXaxis().SetRangeUser(x_min,x_max)
    ratio_hists[i].SetLineWidth(2)
    if ratio_hists[i].GetXaxis().GetXmin() == 0 and ratio_hists[i].GetXaxis().GetXmax() == 0:
      ratio_hists[i].SetLineColor(ROOT.kGray+2)
    else:
      draw_hist(ratio_hists[i], colors[i], 1)
  for i in range(len(hists1)):
    c1.cd(1)
    hists2[i].SetLineWidth(2)
    draw_hist(hists2[i], colors[i], 2)

  if text:
    c1.cd(1)
    text_field = ROOT.TPaveText(0.08,0.9,0.95,0.95, "NDC")
    text_field.AddText(text)
    text_field.SetFillColorAlpha(0,0)
    text_field.SetTextAlign(12)
    text_field.SetBorderSize(0)
    text_field.SetTextSize(0.03*1./0.7)
    text_field.Draw()

  if legend:
    c1.cd(1)
    legend_field = ROOT.TLegend(0.65,0.6,0.95,0.9)
    legend_field.SetFillColorAlpha(0,0)
    legend_field.SetTextAlign(12)
    legend_field.SetBorderSize(0)
    legend_field.SetTextSize(0.03*1./0.7)
    for i in range(len(legend)):
      hists1[i].SetLineColor(colors[i])
      legend_field.AddEntry(hists1[i], legend[i], "L")
    legend_field.Draw()
  
  c1.Print(os.path.splitext(output_path)[0]+".pdf")
  print "Created file "+os.path.splitext(output_path)[0]+".pdf\n"
  c1.Print(os.path.splitext(output_path)[0]+".root")
  print "Created file "+os.path.splitext(output_path)[0]+".root\n"
  c1.Print(os.path.splitext(output_path)[0]+".xml")
  print "Created file "+os.path.splitext(output_path)[0]+".xml\n"
