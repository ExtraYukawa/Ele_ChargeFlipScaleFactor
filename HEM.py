import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion_pure
import numpy as np
import json
from collections import OrderedDict
#histograms name
#hists_name = ['OS_OPS_l1_pt','OS_OPS_l1_eta','OS_OPS_l1_phi','OS_OPS_l2_pt','OS_OPS_l2_eta','OS_OPS_l2_phi','OS_OPS_z_pt','OS_OPS_z_eta','OS_OPS_z_phi','OS_OPS_z_mass','OS_OPS_kinematic_region','SS_ttc_l1_pt','SS_ttc_l1_eta','SS_ttc_l1_phi','SS_ttc_l2_pt','SS_ttc_l2_eta','SS_ttc_l2_phi','SS_ttc_mll','SS_ttc_kinematic_region']
#OS_hists_name = ['OS_OS_OS_OPS_z_mass','OPS_kinematic_region']
#SS_hists_name = ['ttc_mll', 'ttc_kinematic_region']
#hists_name = ['OS_OS_OS_OPS_z_mass', 'OPS_kinematic_region', 'ttc_mll', 'ttc_kinematic_region']

def Add(h1,h2):
  if h1 is None:
    h1 = h2.Clone()
  else:
    h1.Add(h2.Clone())
  h1.SetDirectory(ROOT.gROOT)
  return h1

def TTC_Analysis(era):

  lumi = 59830.

  path = 'flatten/era' + era + '/NotApplyChargeFlipsf_Nominal/'
  dirs = os.listdir(path)

  Data_list = OrderedDict()
  Data_name = ["EGammaA","EGammaB","EGammaC","EGammaD"]
  for name in Data_name:
    Data_list[name] = [0,0,0,0]
  for data in Data_list:
    for f in dirs:
      if data in f:
        fin = ROOT.TFile.Open("flatten/era" + era + '/NotApplyChargeFlipsf_Nominal/' + f)
        h_OS = fin.Get('OS_OPS_isHEM')
        h_SS = fin.Get('SS_ttc_isHEM')
        for i in range(len(Data_list[data])):
          Data_list[data][i] += h_OS.GetBinContent(i+1)
          Data_list[data][i] += h_SS.GetBinContent(i+1)
  print(Data_list)
  
  Total_number = 0
  Total_isHEM  = 0
  Total_veto   = 0
  lumi_BCD      = 0
  print("--------------- Era ------------------")

  for data in Data_list:
    print("* " + data)
    total = (Data_list[data][0] + Data_list[data][1] + Data_list[data][2] + Data_list[data][3])
    inBCD = (Data_list[data][1] + Data_list[data][3])
    isHEM = (Data_list[data][1] + Data_list[data][2])
    veto  = (Data_list[data][1])
    print("  - Total = %7.0f  %3.2f percent"%(total,total/total*100))
    print("  - inBCD = %7.0f  %3.2f percent"%(inBCD,inBCD/total*100))
    print("  - inHEM = %7.0f  %3.2f +/- %1.3f percent"%(isHEM,isHEM/total*100,isHEM/total*sqrt(1./total + 1./isHEM)*100))
    print("  - veto  = %7.0f  %3.2f percent"%(veto,veto/total*100))

    Total_number += total
    Total_isHEM  += isHEM
    Total_veto   += veto
    lumi_BCD     += inBCD
  print("--------------- summary ------------")
  print("* Total: %8.0f"%(Total_number))
  print("* Total in HEM: %6.0f"%(Total_isHEM))
  print("* Veto rate: %f"%(Total_veto/Total_isHEM))
  print("* (subB+C+D)/Total: %f"%(lumi_BCD/Total_number))

if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
  Eras = ['2018']
  for era in Eras:
    TTC_Analysis(era)
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
