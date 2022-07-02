import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion
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

  indir = 'flatten/'
  
  lumi = 0.
  if(era == '2017'):
    lumi = 41480.
  elif(era == '2018'):
    lumi = 59830.
  elif(era == '2016postapv'):
    lumi = 16810.
  elif(era == '2016apv'):
    lumi = 19520.

  signal_list = ['DYnlo','TTTo2L','TTTo1L']
  data_list = ['DoubleEG','SingleEG','EGamma']
  pass_list = ['DY']

  jsonfile = open(os.path.join('data/sample_' + era + 'UL.json'))
  samples = json.load(jsonfile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
  jsonfile.close()

  path = indir + 'era' + era + '/NotApplyChargeFlipsf_Nominal/'
  dirs = os.listdir(path)
  fin_demo = ROOT.TFile.Open(path + 'h0_DYnlo.root')
  hists_name = []
  for tkey in fin_demo.GetListOfKeys():
    try:
      key = tkey.GetName()
      obj = fin_demo.Get(key)
      if obj.InheritsFrom('TH1'):
        hists_name.append(key)
    except:
      pass
  print(hists_name)
  fin_demo.Close()
  for hist in hists_name:
    print("plotting " + hist + " ...")
    histos = []
    h_MC_woSF = None
    h_MC_p1   = None
    h_MC_m1   = None
    for process,desc in samples:
      h_SF = None
      for f in dirs:
        fname = '_'.join(((f.replace('.root','')).split('_'))[1:])
        if (not desc[1] and process in f) or (desc[1] and process == fname):
          f_woSF = ROOT.TFile.Open(indir + "era" + era + '/NotApplyChargeFlipsf_Nominal/' + f)
          f_SF   = ROOT.TFile.Open(indir + "era" + era + '/ApplyChargeFlipsf_Nominal/' + f)
          f_p1   = ROOT.TFile.Open(indir + "era" + era + '/ApplyChargeFlipsf_UP/' + f)
          f_m1   = ROOT.TFile.Open(indir + "era" + era + '/ApplyChargeFlipsf_DOWN/' +f)
          h_SF = Add(h_SF, f_SF.Get(hist).Clone())
          if desc[1]:
            h_MC_woSF = Add(h_MC_woSF, f_woSF.Get(hist).Clone())
            h_MC_p1 = Add(h_MC_p1, f_p1.Get(hist).Clone())
            h_MC_m1 = Add(h_MC_m1, f_m1.Get(hist).Clone())
      histos.append(h_SF)      
#    for i in range(len(process_list)-3):
#      h_woSF_MC.Add(histos_woSF[i+1],1.0)
#      h_p1_MC.Add(histos_p1[i+1],1.0)
#      h_m1_MC.Add(histos_m1[i+1],1.0)
    h_MC_woSF.Scale(lumi)
    h_MC_p1.Scale(lumi)
    h_MC_m1.Scale(lumi)
    c1 = plot_DYregion.draw_plots(histos, 1, hist, 0,h_MC_woSF,h_MC_p1,h_MC_m1,era)

    del histos[:]


if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
#  Eras = ['2016apv','2016postapv','2017','2018']
  Eras = ['2017','2018']
  for era in Eras:
    TTC_Analysis(era)
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
