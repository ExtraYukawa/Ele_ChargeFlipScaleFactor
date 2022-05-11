import ROOT
import math
from math import sqrt
import numpy as np

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

def all_trigger(df):
  all_trigger = df.Filter("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_IsoMu27 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf")
  return all_trigger

def for_egamma_trigger_eechannel(df):
  ditri_ele_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ || (!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf))")
  return ditri_ele_trigger

def for_diele_trigger(df):
  ditri_ele_trigger = df.Filter("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
  return ditri_ele_trigger

def for_singleele_trigger_eechannel(df):
  sin_ele_trigger = df.Filter("!(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL) && !(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) && (HLT_passEle32WPTight || HLT_Ele35_WPTight_Gsf)")
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

def kinematic(l1_pt, l2_pt, l1_eta, l2_eta):
  pt_array = [20., 50., 100., 10000000.]
  eta_array = [0., 0.8, 1.479, 2.4]
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

def select_MC_event(df_MC_tree, add_trigger_SF, add_chargeflip_SF, shift, filters, isOS, era, process):
  if era == '2018':
    df_MC_tree = df_MC_tree.Define("PrefireWeight", "1.");

  if isOS:

    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_MC_tree = df_MC_tree.Define("eeID_SF", "eleID_sf_ee(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_MC_tree = df_MC_tree.Define("OPS_kinematic_region","kinematic(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    if(add_chargeflip_SF and (('DY' in process) or ('TTTo' in process))):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(OPS_kinematic_region,1,-1,"+str(shift)+")")
    else:
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "1.")

  else:

    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_MC_tree = df_MC_tree.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_MC_tree = df_MC_tree.Define("eeID_SF", "eleID_sf_ee(ttc_l1_pt, ttc_l2_pt, ttc_l2_eta, ttc_l2_eta)")
    if(add_chargeflip_SF and (('DY' in process) or ('TTTo' in process))):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(ttc_kinematic_region,1,1,"+str(shift)+")")
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

  return all_trigger(df_MC)

