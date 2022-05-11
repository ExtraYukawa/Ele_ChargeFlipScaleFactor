import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion
import numpy as np
import ChargeFlip_Fit
import json
from ScaleFactor_produce import SF_produce

#histograms name
OS_hists_name = ['OPS_l1_pt','OPS_l1_eta','OPS_l1_phi','OPS_l2_pt','OPS_l2_eta','OPS_l2_phi','OPS_z_pt','OPS_z_eta','OPS_z_phi','OPS_z_mass','OPS_kinematic_region']
SS_hists_name = ['ttc_l1_pt','ttc_l1_eta','ttc_l1_phi','ttc_l2_pt','ttc_l2_eta','ttc_l2_phi','ttc_mll','ttc_kinematic_region']
hists_name = ['OPS_l1_pt','OPS_l1_eta','OPS_l1_phi','OPS_l2_pt','OPS_l2_eta','OPS_l2_phi','OPS_z_pt','OPS_z_eta','OPS_z_phi','OPS_z_mass','OPS_kinematic_region','ttc_l1_pt','ttc_l1_eta','ttc_l1_phi','ttc_l2_pt','ttc_l2_eta','ttc_l2_phi','ttc_mll','ttc_kinematic_region']
#OS_hists_name = ['OPS_z_mass','OPS_kinematic_region']
#SS_hists_name = ['ttc_mll', 'ttc_kinematic_region']
#hists_name = ['OPS_z_mass', 'OPS_kinematic_region', 'ttc_mll', 'ttc_kinematic_region']

def TTC_Analysis(era):

  plotdir = '/eos/user/t/tihsu/plot/Ele_chargeflip_sf/'

  lumi = 0.
  if(era == '2017'):
    lumi = 41480.
  elif(era == '2018'):
    lumi = 59830.
  elif(era == '2016'):
    lumi = 16810.
  elif(era == '2016apv'):
    lumi = 19520.

  signal_list = ['DYnlo','TTTo2L','TTTo1L']
  data_list = ['DoubleEG','SingleEG','EGamma']
  pass_list = ['DY']

  h_data_OS = None
  h_data_SS = None
  h_MC_OS = None
  h_MC_SS = None
  h_back_OS = None
  h_back_SS = None

  path = 'flatten/era' + era + '/NotApplyChargeFlipsf_Nominal/' 
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
    elif(fname in signal_list):
      if h_MC_OS is None:
        h_MC_OS = fin.Get('OS_OPS_kinematic_region').Clone()
        h_MC_OS.SetDirectory(ROOT.gROOT)
        h_MC_OS.Scale(lumi)
        h_MC_SS = fin.Get('SS_ttc_kinematic_region').Clone()
        h_MC_SS.SetDirectory(ROOT.gROOT)
        h_MC_SS.Scale(lumi)
      else:
        h_TT_OS = fin.Get('OS_OPS_kinematic_region').Clone()
        h_TT_SS = fin.Get('SS_ttc_kinematic_region').Clone()
        h_MC_OS.Add(h_TT_OS,lumi)
        h_MC_SS.Add(h_TT_SS,lumi)

    elif(isdata):
      if(h_data_OS is None):
        h_data_OS = fin.Get('OS_OPS_kinematic_region').Clone()
        h_data_OS.SetDirectory(ROOT.gROOT)
      else:
        h = fin.Get('OS_OPS_kinematic_region').Clone()
        h_data_OS.Add(h,1.)
      if(h_data_SS is None):
        h_data_SS = fin.Get('SS_ttc_kinematic_region').Clone()
        h_data_SS.SetDirectory(ROOT.gROOT)
      else:
        h = fin.Get('SS_ttc_kinematic_region').Clone()
        h_data_SS.Add(h,1.)
    else:
      if(h_back_OS is None):
        h_back_OS = (fin.Get('OS_OPS_kinematic_region')).Clone()
        h_back_OS.SetDirectory(ROOT.gROOT)
        h_back_OS.Scale(lumi)
      else:
        h = fin.Get('OS_OPS_kinematic_region').Clone()
        h_back_OS.Add(h,lumi)
      if(h_back_SS is None):
        h_back_SS = fin.Get('SS_ttc_kinematic_region').Clone()
        h_back_SS.SetDirectory(ROOT.gROOT)
        h_back_SS.Scale(lumi)
      else:
        h = fin.Get('SS_ttc_kinematic_region').Clone()
        h_back_SS.Add(h,lumi)
    fin.Close()
  pt_region  = [20., 40., 60., 100., 300.]
  eta_region = [0.,  0.8, 1.479, 2.4]
  h = ChargeFlip_Fit.fit(era, h_MC_OS,h_MC_SS,h_back_OS,h_back_SS,h_data_OS,h_data_SS,pt_region,eta_region,1,1,1)
  SF_nominal = SF_produce(era, h[0], h[1], h[2], h[3], 1, '')
  h2 = ChargeFlip_Fit.fit(era, h_MC_OS, h_MC_SS, h_back_OS, h_back_SS, h_data_OS, h_data_SS, pt_region, eta_region,1,0,0)
  SF_subMC   = SF_produce(era, h2[0], h2[1], h2[2], h2[3], 0, '_subMC')
  h3 = ChargeFlip_Fit.fit(era, h_MC_OS, h_MC_SS, h_back_OS, h_back_SS, h_data_OS, h_data_SS, pt_region, eta_region,0,1,0)
  SF_chi2    = SF_produce(era, h3[0], h3[1], h3[2], h3[3], 0, '_chi2')

  c = ROOT.TCanvas()

  fout = ROOT.TFile.Open("data/ChargeFlipSF_%s_MLE.root"%era,"UPDATE")
  fout.cd()

  for i in range(2):
    legend = ROOT.TLegend(0.1,0.7,0.44448,0.9)
    if(i==1):
      tag="SS"
    else:
      tag="OS"
#    SF_nominal[i].GetYaxis().SetRangeUser(0.95,1.05)
    SF_nominal[i].Draw("pe")

    SF_subMC[i].SetLineColor(2)
    SF_subMC[i].Draw("SAME HIST")

    SF_chi2[i].Draw("SAME HIST")

    legend.AddEntry(SF_nominal[i],"Nominal(MLE)","p")
    legend.AddEntry(SF_subMC[i],"w/o subMC","l")
    legend.AddEntry(SF_chi2[i],"#Chi^2 Fitting","l")
    legend.Draw("SAME")

    c.SaveAs(plotdir + tag + "SF_sys_" + era + ".png" )

    SF_nominal[i].Write()
    SF_subMC[i].Write()
    SF_chi2[i].Write()
    
  fout.Close()
  #ChargeFlip_Fit.fit(h_MC_OS,h_MC_SS,h_back_OS,h_back_SS,h_data_OS,h_data_SS,pt_region,eta_region,0)

if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
  Eras = ['2017','2018']
  for era in Eras:
    TTC_Analysis(era)
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
