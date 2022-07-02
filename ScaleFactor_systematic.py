import ROOT
import os
import math
import numpy as np

def SF_systematic(h_nominal, h_shift):

  h_sys = h_nominal.Clone()
  Nbin = h_nominal.GetNbinsX()

  for i in range(Nbin):
    error = abs(h_nominal.GetBinContent(i+1) -  h_shift.GetBinContent(i+1))
    h_sys.SetBinContent(i+1, error)

  return h_sys
  
def non_closure_test(era, indir, OS_hist, SS_hist):

  lumi = 0.
  if(era == '2017'):
    lumi = 41480.
  elif(era == '2018'):
    lumi = 59830.
  elif(era == '2016postapv'):
    lumi = 16810.
  elif(era == '2016apv'):
    lumi = 19520.

#  signal_list = ['DYnlo','TTTo2L','TTTo1L']
  data_list = ['DoubleEG','SingleEG','EGamma']
  pass_list = ['DY']
  back_list = ['TTTo1L']

  h_data_OS = None
  h_data_SS = None
  h_MC_OS = None
  h_MC_SS = None
  h_back_OS = None
  h_back_SS = None

  path = indir + '/era' + era + '/ApplyChargeFlipsf_Nominal/'
  dirs = os.listdir(path)
 
  for f in dirs:
    fname = '_'.join(((f.replace('.root','')).split('_'))[1:])
    fin = ROOT.TFile.Open(path + f)
    isdata = False
    for data in data_list:
      if(data in fname):
        isdata = True
    if(fname in pass_list):
      pass
    elif(fname not in back_list and not isdata):
      if h_MC_OS is None:
        h_MC_OS = fin.Get(OS_hist).Clone()
        h_MC_OS.SetDirectory(ROOT.gROOT)
        h_MC_OS.Scale(lumi)
        h_MC_SS = fin.Get(SS_hist).Clone()
        h_MC_SS.SetDirectory(ROOT.gROOT)
        h_MC_SS.Scale(lumi)
      else:
        h_TT_OS = fin.Get(OS_hist).Clone()
        h_TT_SS = fin.Get(SS_hist).Clone()
        h_MC_OS.Add(h_TT_OS,lumi)
        h_MC_SS.Add(h_TT_SS,lumi)

    elif(isdata):
      if(h_data_OS is None):
        h_data_OS = fin.Get(OS_hist).Clone()
        h_data_OS.SetDirectory(ROOT.gROOT)
      else:
        h = fin.Get(OS_hist).Clone()
        h_data_OS.Add(h,1.)
      if(h_data_SS is None):
        h_data_SS = fin.Get(SS_hist).Clone()
        h_data_SS.SetDirectory(ROOT.gROOT)
      else:
        h = fin.Get(SS_hist).Clone()
        h_data_SS.Add(h,1.)
    else:
      if(h_back_OS is None):
        h_back_OS = (fin.Get(OS_hist)).Clone()
        h_back_OS.SetDirectory(ROOT.gROOT)
        h_back_OS.Scale(lumi)
      else:
        h = fin.Get(OS_hist).Clone()
        h_back_OS.Add(h,lumi)
      if(h_back_SS is None):
        h_back_SS = fin.Get(SS_hist).Clone()
        h_back_SS.SetDirectory(ROOT.gROOT)
        h_back_SS.Scale(lumi)
      else:
        h = fin.Get(SS_hist).Clone()
        h_back_SS.Add(h,lumi)
      fin.Close()

  Nbin = h_data_SS.GetNbinsX()

  uncertainty = 0.
  summation      = 0.

  for i in range(Nbin):
    data_number = h_data_SS.GetBinContent(i+1) 
    back_number = h_back_SS.GetBinContent(i+1)
    MC_number   = h_MC_SS.GetBinContent(i+1)
    lowEdge     = h_data_SS.GetBinLowEdge(i+1)
    highEdge    = h_data_SS.GetBinWidth(i+1)+lowEdge
#    MC_number  += back_number
#    back_number = 0
#    print(data_number)
#    print(MC_number)

    if data_number > 0 and MC_number > 1:
      bin_uncertainty = abs((data_number - back_number - MC_number)/MC_number)
#      print("%.2f: %.1f(%.1f ~ %.1f)"%(bin_uncertainty, data_number,lowEdge,highEdge))
      a = data_number
      b = back_number
      c = MC_number
      weight = 1./(a/(c*c)+b/(c*c)+(c+a-b-c)**2/(c**3))
      #weight = c
      uncertainty +=  bin_uncertainty*weight# unc. / (unc._data)^2 --> unc._data from Poisson dist.
      summation   +=  weight
      print("%.2f: %.4f, %.1f(%.1f ~ %.1f)"%(bin_uncertainty, weight, data_number,lowEdge,highEdge))


  uncertainty /= summation
  return uncertainty


def SF_produce(era):

  fin = ROOT.TFile.Open("data/ChargeFlipSF_%s_MLE.root"%era,"UPDATE")
  channel = ['SS']

  for ch in channel:
    h_sys_list = []

    h_nominal  = fin.Get(ch + "_ChargeFlip_SF")
    h_chi2     = fin.Get(ch + "_ChargeFlip_SF_chi2")
    h_subMC    = fin.Get(ch + "_ChargeFlip_SF_subMC")

    h_sys_subMC =  SF_systematic(h_nominal, h_subMC)    
    
    src_hist = ["SS_ttc_l1_pt", "SS_ttc_l2_pt","SS_ttc_l1_eta","SS_ttc_l2_eta"]
    overall_uncertainty = 0.
    for h in src_hist:
      uncertainty = non_closure_test(era,'validation',h,h)
      print(h+": %f"%uncertainty)
      overall_uncertainty = max(uncertainty,overall_uncertainty)
    print(era + ": %f"%overall_uncertainty)
    print("------------------------------")

    Nbin = h_nominal.GetNbinsX()

    h_final = h_nominal.Clone()
    h_final.SetName("SS_ChargeFlip_SF_AllUnc")

    for i in range(Nbin):
      stat_err = h_nominal.GetBinError(i+1)
      over_err = h_nominal.GetBinContent(i+1)*overall_uncertainty
      #sys_err_subMC = h_sys_subMC.GetBinContent(i+1)
      total_err = (stat_err*stat_err + over_err*over_err)**0.5
      h_final.SetBinError(i+1,total_err)

    fin.cd()
    h_final.Write()
    fin.Close()

      

if __name__ == '__main__':
  Eras = ['2016apv','2016postapv','2017','2018']
#  Eras = ['2018']
  for era in Eras:
    SF_produce(era)



