import ROOT
import math
from math import sqrt
import numpy as np
import re

def overunder_flowbin(h1):
  h1.SetBinContent(1,h1.GetBinContent(0)+h1.GetBinContent(1))
  h1.SetBinError(1,sqrt(h1.GetBinError(0)*h1.GetBinError(0)+h1.GetBinError(1)*h1.GetBinError(1)))
  h1.SetBinContent(h1.GetNbinsX(),h1.GetBinContent(h1.GetNbinsX())+h1.GetBinContent(h1.GetNbinsX()+1))
  h1.SetBinError(h1.GetNbinsX(),sqrt(h1.GetBinError(h1.GetNbinsX())*h1.GetBinError(h1.GetNbinsX())+h1.GetBinError(h1.GetNbinsX()+1)*h1.GetBinError(h1.GetNbinsX()+1)))
  return h1

def get_mcEventnumber(filename):
  print 'opening file ', filename
  nevent_temp=0
  for i in range(0,len(filename)):
    ftemp=ROOT.TFile.Open(filename[i])
    htemp=ftemp.Get('nEventsGenWeighted')
    nevent_temp=nevent_temp+htemp.GetBinContent(1)
  return nevent_temp

def all_trigger(df,era):
  all_trigger = None
  if("2016" not in era):
    all_trigger = df.Filter("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_IsoMu27 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf")
  else:
    all_trigger = df.Filter("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_IsoMu27 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_passEle32WPTight || HLT_Ele35_WPLoose_Gsf")
  return all_trigger

def for_egamma_trigger_eechannel(df,era):
  ditri_ele_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || (!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf))")
  return ditri_ele_trigger

def for_diele_trigger(df,era):
  ditri_ele_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
  return ditri_ele_trigger

def for_singleele_trigger_eechannel(df,era):
  sin_ele_trigger = None
  if("2016" not in era):
    sin_ele_trigger = df.Filter("!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) && !(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf)")
  else:
    sin_ele_trigger = df.Filter("!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) && !(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPLoose_Gsf)")
  return sin_ele_trigger

def for_singleele_trigger_emuchannel(df):
  sin_ele_trigger = df.Filter("!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ) && !(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ) && !(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf)")
  return sin_ele_trigger

def for_dimuon_trigger(df):
  ditri_mu_trigger = df.Filter("(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)")
  return ditri_mu_trigger

def for_singlemuon_trigger_mumuchannel(df):
  single_mu_trigger = df.Filter("!(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8) && HLT_IsoMu27")
  return single_mu_trigger

def for_singlemuon_trigger_emuchannel(df):
  single_mu_trigger = df.Filter("!(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && HLT_IsoMu27")
  return single_mu_trigger

def for_cross_trigger(df):
  x_trigger = df.Filter("(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ)")
  return x_trigger

def data_trigger(df, channel, subera, dataset, run_dict, com_dict):
   DiLepton_slc_run = dict()
   Run_List = run_dict["Data"][channel][subera]
   for Name in Run_List.keys():
     DiLepton_slc_run[Name] = ROOT.std.vector('int')()
     for i in Run_List[Name]:
       DiLepton_slc_run[Name].push_back(i)
   if dataset == "EGamma":
     Trigger = "(" + str(com_dict[channel]["Data"][subera]["DoubleEG"]) + ")||(" + str(com_dict[channel]["Data"][subera]["SingleEG"]) + ")"
   else:
     Trigger = com_dict[channel]["Data"][subera][dataset]
   p1 = re.compile(r'[{](.*?)[}]', re.S)  
   variables = re.findall(p1,Trigger)
   var_list = []
   for var in variables:
      Trigger = Trigger.replace(var,"")
      runs = eval(var)
      #runs = [1,2,3]
      runs = [str(run) for run in runs]
      runs = ','.join(runs)
      run_command = "{" + runs + "}"
      var_list.append(run_command)
      
   Trigger = Trigger.format(*var_list)
#   Trigger2 = "Triggers(run , HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,{1,2,3})"
#   print(Trigger)
#   print(Trigger == Trigger2)
   
   data_trigger = df.Filter(str(Trigger))
   return data_trigger
def MC_trigger(df, channel, com_dict):
  command = str(com_dict[channel]["MC"])
  print(command)
  return df.Filter(command)
def kinematic(l1_pt, l2_pt, l1_eta, l2_eta):
  pt_array = [20., 50., 100., 10000000.]
  eta_array = [0., 0.8, 1.479, 2.5]
  pt_input = [l1_pt,l2_pt]
  eta_input = [abs(l1_eta), abs(l2_eta)]
  pt_bin = np.digitize(pt_input, pt_array)
  eta_bin = np.digitize(eta_input, eta_array)
  for i in range(2):
    if(pt_bin[i]<1 or pt_bin[i]>=len(pt_array)):
      return -1
    if(eta_bin[i]<1 or eta_bin[i]>=len(eta_array)):
      return -1
  xbin = len(pt_array)-1
  ybin = len(eta_array)-1
  return ((eta_bin[1]-1) + ybin*(pt_bin[1]-1) + ybin*xbin*(eta_bin[0]-1) + ybin*xbin*ybin*(pt_bin[0]-1))

def select_MC_event(df_MC_tree, add_trigger_SF, add_chargeflip_SF, shift, filters, isOS, era, process, channel, com_dict):
  if era == '2018':
    df_MC_tree = df_MC_tree.Define("PrefireWeight", "1.");

  if isOS:

    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_MC_tree = df_MC_tree.Define("eeID_SF", "eleID_sf_ee(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_MC_tree = df_MC_tree.Define("OPS_kinematic_region","kinematic(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    if(add_chargeflip_SF):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(OPS_kinematic_region,1,OPS_genOS,"+str(shift)+")")
    else:
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "1.")

  else:

    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_MC_tree = df_MC_tree.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_MC_tree = df_MC_tree.Define("eeID_SF", "eleID_sf_ee(ttc_l1_pt, ttc_l2_pt, ttc_l2_eta, ttc_l2_eta)")
    if(add_chargeflip_SF):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(ttc_kinematic_region,0,ttc_genOS,"+str(shift)+")")
    else:
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "1.")

  if not add_trigger_SF:

    if isOS:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[OPS_l1_id]*Electron_RECO_SF[OPS_l2_id]*eeID_SF*charge_flip_SF*genWeight/abs(genWeight)")
    else:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[ttc_l1_id]*Electron_RECO_SF[ttc_l2_id]*eeID_SF*charge_flip_SF*genWeight/abs(genWeight)")

  else:

    if isOS:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[OPS_l1_id]*Electron_RECO_SF[OPS_l2_id]*eeID_SF*trigger_SF*charge_flip_SF*genWeight/abs(genWeight)")
    else:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[ttc_l1_id]*Electron_RECO_SF[ttc_l2_id]*eeID_SF*trigger_SF*charge_flip_SF*genWeight/abs(genWeight)")
  df_MC = df_MC_tree.Filter(filters)

  return MC_trigger(df_MC,channel,com_dict)

