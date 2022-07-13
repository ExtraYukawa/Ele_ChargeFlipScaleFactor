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

def weighted(h,eta_region,pt_region,eta_index,pt_index,number):
  endcap = 0.
  barrel = 0.
  n_endcap = 0.
  n_barrel = 0.

  CFrate = h.GetBinContent(pt_index+1,eta_index+1)
  if (eta_region[eta_index] >= 1.479):
    endcap = CFrate*number
    n_endcap = number
  else:
    barrel = CFrate*number
    n_barrel = number
  return barrel,endcap, n_barrel, n_endcap

    

def Fill_histogram(h_data_OS,h_data_SS,h_MC_OS,h_MC_SS,h_back_OS,h_back_SS,era,data_list,signal_list,pass_list):
  lumi = 0.
  if(era == '2017'):
    lumi = 41480.
  elif(era == '2018'):
    lumi = 59830.
  elif(era == '2016postapv'):
    lumi = 16810.
  elif(era == '2016apv'):
    lumi = 19520.
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
  return h_data_OS,h_data_SS,h_MC_OS,h_MC_SS,h_back_OS,h_back_SS

def TTC_Analysis(eras):

  era = eras[0]#[:4]
  plotdir = '/eos/user/t/tihsu/plot/Ele_chargeflip_sf_eratrig/'

  signal_list = ['DYnlo','TTTo2L']
  data_list = ['DoubleEG','SingleEG','EGamma']
  pass_list = ['DY']

  h_data_OS = None
  h_data_SS = None
  h_MC_OS = None
  h_MC_SS = None
  h_back_OS = None
  h_back_SS = None

  for Era in eras:
    h_data_OS,h_data_SS,h_MC_OS,h_MC_SS,h_back_OS,h_back_SS = Fill_histogram(h_data_OS,h_data_SS,h_MC_OS,h_MC_SS,h_back_OS,h_back_SS,Era,data_list,signal_list,pass_list)

  pt_region  = [20., 40., 60., 100., 300.]
  eta_region = [0.,  0.8, 1.479, 2.5]
  h_data, h_MC, h_data_cov, h_MC_cov = ChargeFlip_Fit.fit(era, h_MC_OS,h_MC_SS,h_back_OS,h_back_SS,h_data_OS,h_data_SS,pt_region,eta_region,1,0,0)
  #useLikelihood, subMC, draw

  pt_bins = len(pt_region)-1
  eta_bins = len(eta_region)-1


  sig_OS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  sig_SS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  data_OS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  data_SS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))

  barrel_data_number = 0.
  barrel_MC_number   = 0.
  barrel_data_CFrate = 0.
  barrel_MC_CFrate   = 0.
  endcap_data_number = 0.
  endcap_MC_number   = 0.
  endcap_data_CFrate = 0.
  endcap_MC_CFrate   = 0.

  for i in range(pt_bins):
   for j in range(eta_bins):
     for ii in range(pt_bins):
       for jj in range(eta_bins):
         index = jj+ii*eta_bins+j*eta_bins*pt_bins+i*eta_bins*pt_bins*eta_bins
         MC_number = h_MC_OS.GetBinContent(h_MC_OS.FindBin(index)) + h_MC_SS.GetBinContent(h_MC_SS.FindBin(index))
         data_number = h_data_OS.GetBinContent(h_data_OS.FindBin(index)) + h_data_SS.GetBinContent(h_data_SS.FindBin(index))
         barrel, endcap, n_barrel, n_endcap = weighted(h_MC,eta_region,pt_region,j,i,MC_number)
         barrel_MC_CFrate += barrel
         endcap_MC_CFrate += endcap
         barrel_MC_number += n_barrel
         endcap_MC_number += n_endcap
         barrel, endcap, n_barrel, n_endcap = weighted(h_MC,eta_region,pt_region,jj,ii,MC_number)
         barrel_MC_CFrate += barrel
         endcap_MC_CFrate += endcap
         barrel_MC_number += n_barrel
         endcap_MC_number += n_endcap

         barrel, endcap, n_barrel, n_endcap = weighted(h_data,eta_region,pt_region,j,i,data_number)
         barrel_data_CFrate += barrel
         endcap_data_CFrate += endcap
         barrel_data_number += n_barrel
         endcap_data_number += n_endcap
         barrel, endcap, n_barrel, n_endcap = weighted(h_data,eta_region,pt_region,jj,ii, data_number)
         barrel_data_CFrate += barrel
         endcap_data_CFrate += endcap
         barrel_data_number += n_barrel
         endcap_data_number += n_endcap

        
  barrel_MC_CFrate /= barrel_MC_number
  endcap_MC_CFrate /= endcap_MC_number
  barrel_data_CFrate /= barrel_data_number
  endcap_data_CFrate /= endcap_data_number
  return barrel_MC_CFrate, endcap_MC_CFrate, barrel_data_CFrate, endcap_data_CFrate


if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
  Eras = [['2016postapv'],['2016apv'],['2017'],['2018']]
  Barrel_MC_CFrate = dict()
  Endcap_MC_CFrate = dict()
  Barrel_data_CFrate = dict()
  Endcap_data_CFrate = dict()

  for eras in Eras:
    barrel_MC_CFrate, endcap_MC_CFrate, barrel_data_CFrate, endcap_data_CFrate = TTC_Analysis(eras)
    Barrel_MC_CFrate[eras[0]] = barrel_MC_CFrate
    Endcap_MC_CFrate[eras[0]] = endcap_MC_CFrate
    Barrel_data_CFrate[eras[0]] = barrel_data_CFrate
    Endcap_data_CFrate[eras[0]] = endcap_data_CFrate


  for eras in Eras:
    print("*** "+eras[0])
    print("----- Barrel MC: %.2e"%Barrel_MC_CFrate[eras[0]])
    print("----- Endcap MC: %.2e"%Endcap_MC_CFrate[eras[0]])
    print("----- Barrel data: %.2e"%Barrel_data_CFrate[eras[0]])
    print("----- Endcap data: %.2e"%Endcap_data_CFrate[eras[0]])
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
