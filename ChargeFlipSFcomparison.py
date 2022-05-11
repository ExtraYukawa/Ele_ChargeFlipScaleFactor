import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion
import numpy as np

TTC_header_path = os.path.join("TTC.h")
ROOT.gInterpreter.Declare('#include "{}"'.format(TTC_header_path))


# the EnableImplicitMT option should only use in cluster, at lxplus, it will make the code slower(my experience)
#ROOT.ROOT.EnableImplicitMT()

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

def select_MC_event(df_MC_tree, add_trigger_SF, filters, isOS, sigma, process):

  if isOS:
    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(DY_l1_pt, DY_l2_pt, DY_l1_eta, DY_l2_eta)")
    df_MC_tree = df_MC_tree.Define("DY_kinematic_region","kinematic(DY_l1_pt, DY_l2_pt, DY_l1_eta, DY_l2_eta)")
    if(('DY' in process) or ('TTTo' in process)):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(DY_kinematic_region,1,-1,"+str(sigma)+")")
    else:
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "1.")
  else:
    df_MC_tree = df_MC_tree.Define("trigger_SF", "trigger_sf_ee(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_MC_tree = df_MC_tree.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    if(('DY' in process) or ('TTTo' in process)):
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "chargeflip_sf(ttc_kinematic_region, 1,1, "+str(sigma)+")")
    else:
      df_MC_tree = df_MC_tree.Define("charge_flip_SF", "1.")
  if not add_trigger_SF:
    if isOS:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[DY_l1_id]*Electron_RECO_SF[DY_l2_id]*Electron_CutBased_TightID_SF[DY_l1_id]*Electron_CutBased_TightID_SF[DY_l2_id]*genWeight/abs(genWeight)")
    else:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[ttc_l1_id]*Electron_RECO_SF[ttc_l2_id]*Electron_CutBased_TightID_SF[ttc_l1_id]*Electron_CutBased_TightID_SF[ttc_l2_id]*genWeight/abs(genWeight)")
  else:
    if isOS:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[DY_l1_id]*Electron_RECO_SF[DY_l2_id]*Electron_CutBased_TightID_SF[DY_l1_id]*Electron_CutBased_TightID_SF[DY_l2_id]*trigger_SF*charge_flip_SF*genWeight/abs(genWeight)")
    else:
      df_MC_tree = df_MC_tree.Define("genweight","puWeight*PrefireWeight*Electron_RECO_SF[ttc_l1_id]*Electron_RECO_SF[ttc_l2_id]*Electron_CutBased_TightID_SF[ttc_l1_id]*Electron_CutBased_TightID_SF[ttc_l2_id]*trigger_SF*charge_flip_SF*genWeight/abs(genWeight)")

  df_MC = df_MC_tree.Filter(filters)
  return all_trigger(df_MC)


path='/eos/user/m/melu/TTC_Nanov8_new/'
add_trigger_SF=True

doubleMu_names = ROOT.std.vector('string')()
for f in ["DoubleMuonB.root","DoubleMuonC.root","DoubleMuonD.root","DoubleMuonE.root","DoubleMuonF.root"]:
  doubleMu_names.push_back(path+f)

singleMu_names = ROOT.std.vector('string')()
for f in ["SingleMuonB.root","SingleMuonC.root","SingleMuonD.root","SingleMuonE.root","SingleMuonF.root"]:
  singleMu_names.push_back(path+f)

doubleEle_names = ROOT.std.vector('string')()
for f in ["DoubleEGB.root","DoubleEGC.root","DoubleEGD.root","DoubleEGE.root","DoubleEGF.root"]:
  doubleEle_names.push_back(path+f)

singleEle_names = ROOT.std.vector('string')()
for f in ["SingleEGB.root","SingleEGC.root","SingleEGD.root","SingleEGE.root","SingleEGF.root"]:
  singleEle_names.push_back(path+f)

muonEle_names = ROOT.std.vector('string')()
for f in ["MuonEGB.root","MuonEGC.root","MuonEGD.root","MuonEGE.root","MuonEGF.root"]:
  muonEle_names.push_back(path+f)

DY_list = ROOT.std.vector('string')()
for f in ['DY.root']:
  DY_list.push_back(path+f)

WJet_list = ROOT.std.vector('string')()
for f in ['WJets.root']:
  WJet_list.push_back(path+f)

WW_list = ROOT.std.vector('string')()
for f in ['WW.root']:
  WW_list.push_back(path+f)

WZ_list = ROOT.std.vector('string')()
for f in ['WZ.root']:
  WZ_list.push_back(path+f)

ZZ_list = ROOT.std.vector('string')()
for f in ['ZZ.root']:
  ZZ_list.push_back(path+f)

WWW_list = ROOT.std.vector('string')()
for f in ['WWW.root']:
  WWW_list.push_back(path+f)

WWZ_list = ROOT.std.vector('string')()
for f in ['WWZ.root']:
  WWZ_list.push_back(path+f)

WZZ_list = ROOT.std.vector('string')()
for f in ['WZZ.root']:
  WZZ_list.push_back(path+f)

ZZZ_list = ROOT.std.vector('string')()
for f in ['ZZZ.root']:
  ZZZ_list.push_back(path+f)

tsch_list = ROOT.std.vector('string')()
for f in ['tsch.root']:
  tsch_list.push_back(path+f)

t_tch_list = ROOT.std.vector('string')()
for f in ['t_tch.root']:
  t_tch_list.push_back(path+f)

tbar_tch_list = ROOT.std.vector('string')()
for f in ['tbar_tch.root']:
  tbar_tch_list.push_back(path+f)

tW_list = ROOT.std.vector('string')()
for f in ['tW.root']:
  tW_list.push_back(path+f)

tbarW_list = ROOT.std.vector('string')()
for f in ['tbarW.root']:
  tbarW_list.push_back(path+f)

ttWtoLNu_list = ROOT.std.vector('string')()
for f in ['ttWtoLNu.root']:
  ttWtoLNu_list.push_back(path+f)

ttWtoQQ_list = ROOT.std.vector('string')()
for f in ['ttWtoQQ.root']:
  ttWtoQQ_list.push_back(path+f)

ttZ_list = ROOT.std.vector('string')()
for f in ['ttZ.root']:
  ttZ_list.push_back(path+f)

ttZtoQQ_list = ROOT.std.vector('string')()
for f in ['ttZtoQQ.root']:
  ttZtoQQ_list.push_back(path+f)

ttH_list = ROOT.std.vector('string')()
for f in ['ttH.root']:
  ttH_list.push_back(path+f)

ttWW_list = ROOT.std.vector('string')()
for f in ['ttWW.root']:
  ttWW_list.push_back(path+f)

ttWZ_list = ROOT.std.vector('string')()
for f in ['ttWZ.root']:
  ttWZ_list.push_back(path+f)

ttZZ_list = ROOT.std.vector('string')()
for f in ['ttZZ.root']:
  ttZZ_list.push_back(path+f)

tzq_list = ROOT.std.vector('string')()
for f in ['tzq.root']:
  tzq_list.push_back(path+f)

TTTo2L_list = ROOT.std.vector('string')()
for f in ['TTTo2L.root']:
  TTTo2L_list.push_back(path+f)

TTTo1L_list = ROOT.std.vector('string')()
for f in ['TTTo1L.root']:
  TTTo1L_list.push_back(path+f)

#QCD50to80_list = ROOT.std.vector('string')()
#for f in ['QCD50to80.root']:
#  QCD50to80_list.push_back(path+f)
#
#QCD80to120_list = ROOT.std.vector('string')()
#for f in ['QCD80to120.root']:
#  QCD80to120_list.push_back(path+f)
#
#QCD120to170_list = ROOT.std.vector('string')()
#for f in ['QCD120to170.root']:
#  QCD120to170_list.push_back(path+f)
#
#QCD170to300_list = ROOT.std.vector('string')()
#for f in ['QCD170to300.root']:
#  QCD170to300_list.push_back(path+f)
#
#QCD300toinf_list = ROOT.std.vector('string')()
#for f in ['QCD300toinf.root']:
#  QCD300toinf_list.push_back(path+f)

#histograms name
OS_hists_name = ['DY_l1_pt','DY_l1_eta','DY_l1_phi','DY_l2_pt','DY_l2_eta','DY_l2_phi','DY_z_pt','DY_z_eta','DY_z_phi','DY_z_mass','DY_kinematic_region']
SS_hists_name = ['ttc_l1_pt','ttc_l1_eta','ttc_l1_phi','ttc_l2_pt','ttc_l2_eta','ttc_l2_phi','ttc_mll','ttc_kinematic_region']
hists_name = ['DY_l1_pt','DY_l1_eta','DY_l1_phi','DY_l2_pt','DY_l2_eta','DY_l2_phi','DY_z_pt','DY_z_eta','DY_z_phi','DY_z_mass','DY_kinematic_region','ttc_l1_pt','ttc_l1_eta','ttc_l1_phi','ttc_l2_pt','ttc_l2_eta','ttc_l2_phi','ttc_mll','ttc_kinematic_region']
#OS_hists_name = ['DY_z_mass','DY_kinematic_region']
#SS_hists_name = ['ttc_mll', 'ttc_kinematic_region']
#hists_name = ['DY_z_mass', 'DY_kinematic_region', 'ttc_mll', 'ttc_kinematic_region']

#histograms bins [nbins, low edge, high edge]
#histos_bins = {
#OS_hists_name[0]:[40,70,110],
#OS_hists_name[1]:[150,-1 ,149],
#SS_hists_name[0]:[40,70,110],
#SS_hists_name[1]:[150,-1, 149]
#}

histos_bins = {
OS_hists_name[0]:[20,0,200],
OS_hists_name[1]:[16,-2.4,2.4],
OS_hists_name[2]:[20,-4,4],
OS_hists_name[3]:[20,0,100],
OS_hists_name[4]:[16,-2.4,2.4],
OS_hists_name[5]:[20,-4,4],
OS_hists_name[6]:[50,0,200],
OS_hists_name[7]:[20,-3,3],
OS_hists_name[8]:[20,-4,4],
OS_hists_name[9]:[60,60,120],
OS_hists_name[10]:[154,0,154],
SS_hists_name[0]:[20,0,200],
SS_hists_name[1]:[16,-2.4,2.4],
SS_hists_name[2]:[20,-4,4],
SS_hists_name[3]:[20,0,100],
SS_hists_name[4]:[16,-2.4,2.4],
SS_hists_name[5]:[20,-4,4],
SS_hists_name[6]:[60,60,120],
SS_hists_name[7]:[154,0,154]
}

def TTC_Analysis(sigma):

  histos = []

  lumi = 41480.

  process_list = []

  DY_xs = 6077.22
  DY_ev = get_mcEventnumber(DY_list)
  process_list.append({'name':'DY', 'list':DY_list, 'xs':DY_xs, 'ev':DY_ev, 'isMC':True, 'hist':[]})

  WJet_xs = 61526.7
  WJet_ev = get_mcEventnumber(WJet_list)
  process_list.append({'name':'WJet', 'list':WJet_list, 'xs':WJet_xs, 'ev':WJet_ev, 'isMC':True, 'hist':[]})

  WW_xs = 118.7
  WW_ev = get_mcEventnumber(WW_list)
  process_list.append({'name':'WW', 'list':WW_list, 'xs':WW_xs, 'ev':WW_ev, 'isMC':True, 'hist':[]})

  WZ_xs = 65.5443
  WZ_ev = get_mcEventnumber(WZ_list)
  process_list.append({'name':'WZ', 'list':WZ_list, 'xs':WZ_xs, 'ev':WZ_ev, 'isMC':True, 'hist':[]})

  ZZ_xs = 15.8274
  ZZ_ev = get_mcEventnumber(ZZ_list)
  process_list.append({'name':'ZZ', 'list':ZZ_list, 'xs':ZZ_xs, 'ev':ZZ_ev, 'isMC':True, 'hist':[]})

  WWW_xs = 0.2086
  WWW_ev = get_mcEventnumber(WWW_list)
  process_list.append({'name':'WWW', 'list':WWW_list, 'xs':WWW_xs, 'ev':WWW_ev, 'isMC':True, 'hist':[]})

  WWZ_xs = 0.1707
  WWZ_ev = get_mcEventnumber(WWZ_list)
  process_list.append({'name':'WWZ', 'list':WWZ_list, 'xs':WWZ_xs, 'ev':WWZ_ev, 'isMC':True, 'hist':[]})

  WZZ_xs = 0.05709
  WZZ_ev = get_mcEventnumber(WZZ_list)
  process_list.append({'name':'WZZ', 'list':WZZ_list, 'xs':WZZ_xs, 'ev':WZZ_ev, 'isMC':True, 'hist':[]})

  ZZZ_xs = 0.01476
  ZZZ_ev = get_mcEventnumber(ZZZ_list)
  process_list.append({'name':'ZZZ', 'list':ZZZ_list, 'xs':ZZZ_xs, 'ev':ZZZ_ev, 'isMC':True, 'hist':[]})

  TTTo2L_xs = 88.3419
  TTTo2L_ev = get_mcEventnumber(TTTo2L_list)
  process_list.append({'name':'TTTo2L', 'list':TTTo2L_list, 'xs':TTTo2L_xs, 'ev':TTTo2L_ev, 'isMC':True, 'hist':[]})

  TTTo1L_xs = 365.4574
  TTTo1L_ev = get_mcEventnumber(TTTo1L_list)
  process_list.append({'name':'TTTo1L', 'list':TTTo1L_list, 'xs':TTTo1L_xs, 'ev':TTTo1L_ev, 'isMC':True, 'hist':[]})

  TTH_xs = 0.5269
  TTH_ev = get_mcEventnumber(ttH_list)
  process_list.append({'name':'TTH', 'list':ttH_list, 'xs':TTH_xs, 'ev':TTH_ev, 'isMC':True, 'hist':[]})

  TTWtoLNu_xs = 0.1792
  TTWtoLNu_ev = get_mcEventnumber(ttWtoLNu_list)
  process_list.append({'name':'TTWtoLNu', 'list':ttWtoLNu_list, 'xs':TTWtoLNu_xs, 'ev':TTWtoLNu_ev, 'isMC':True, 'hist':[]})

  TTWtoQQ_xs = 0.3708
  TTWtoQQ_ev = get_mcEventnumber(ttWtoQQ_list)
  process_list.append({'name':'TTWtoQQ', 'list':ttWtoQQ_list, 'xs':TTWtoQQ_xs, 'ev':TTWtoQQ_ev, 'isMC':True, 'hist':[]})

  TTZ_xs = 0.2589
  TTZ_ev = get_mcEventnumber(ttZ_list)
  process_list.append({'name':'TTZ', 'list':ttZ_list, 'xs':TTZ_xs, 'ev':TTZ_ev, 'isMC':True, 'hist':[]})

  TTZtoQQ_xs = 0.6012
  TTZtoQQ_ev = get_mcEventnumber(ttZtoQQ_list)
  process_list.append({'name':'TTZtoQQ', 'list':ttZtoQQ_list, 'xs':TTZtoQQ_xs, 'ev':TTZtoQQ_ev, 'isMC':True, 'hist':[]})

  TTWW_xs = 0.007003
  TTWW_ev = get_mcEventnumber(ttWW_list)
  process_list.append({'name':'TTWW', 'list':ttWW_list, 'xs':TTWW_xs, 'ev':TTWW_ev, 'isMC':True, 'hist':[]})

  TTWZ_xs = 0.002453
  TTWZ_ev = get_mcEventnumber(ttWZ_list)
  process_list.append({'name':'TTWZ', 'list':ttWZ_list, 'xs':TTWZ_xs, 'ev':TTWZ_ev, 'isMC':True, 'hist':[]})

  TTZZ_xs = 0.001386
  TTZZ_ev = get_mcEventnumber(ttZZ_list)
  process_list.append({'name':'TTZZ', 'list':ttZZ_list, 'xs':TTZZ_xs, 'ev':TTZZ_ev, 'isMC':True, 'hist':[]})

  tZq_xs = 0.07561
  tZq_ev = get_mcEventnumber(tzq_list)
  process_list.append({'name':'tZq', 'list':tzq_list, 'xs':tZq_xs, 'ev':tZq_ev, 'isMC':True, 'hist':[]})

  tW_xs = 35.85
  tW_ev = get_mcEventnumber(tW_list)
  process_list.append({'name':'tW', 'list':tW_list, 'xs':tW_xs, 'ev':tW_ev, 'isMC':True, 'hist':[]})

  tbarW_xs = 35.85
  tbarW_ev = get_mcEventnumber(tbarW_list)
  process_list.append({'name':'tbarW', 'list':tbarW_list, 'xs':tbarW_xs, 'ev':tbarW_ev, 'isMC':True, 'hist':[]})

  t_sch_xs = 3.36
  t_sch_ev = get_mcEventnumber(tsch_list)
  process_list.append({'name':'t_sch', 'list':tsch_list, 'xs':t_sch_xs, 'ev':t_sch_ev, 'isMC':True, 'hist':[]})

  t_tch_xs = 136.02
  t_tch_ev = get_mcEventnumber(t_tch_list)
  process_list.append({'name':'t_tch', 'list':t_tch_list, 'xs':t_tch_xs, 'ev':t_tch_ev, 'isMC':True, 'hist':[]})

  tbar_tch_xs = 80.95
  tbar_tch_ev = get_mcEventnumber(tbar_tch_list)
  process_list.append({'name':'tbar_tch', 'list':tbar_tch_list, 'xs':tbar_tch_xs, 'ev':tbar_tch_ev, 'isMC':True, 'hist':[]})


#  QCD50to80_xs = 1984000.0
#  QCD50to80_ev = get_mcEventnumber(QCD50to80_list)
#
#  QCD80to120_xs = 366500.0
#  QCD80to120_ev = get_mcEventnumber(QCD80to120_list)
#
#  QCD120to170_xs = 66490.0
#  QCD120to170_ev = get_mcEventnumber(QCD120to170_list)
#
#  QCD170to300_xs = 16480.0
#  QCD170to300_ev = get_mcEventnumber(QCD170to300_list)
#
#  QCD300toinf_xs = 1097.0
#  QCD300toinf_ev = get_mcEventnumber(QCD300toinf_list)

####################
#  Z MASS FITTING ##
####################

  OS_Zmass = 90.
  SS_Zmass = 90.
  OS_width = 30.
  SS_width = 30.


  # define the filters here, 1:2mu, 2:1e1m, 3:2ele
#  OS_filters="DY_region==3 && DY_z_mass>60 && DY_z_mass<120 && (DY_l1_pt>30 || DY_l2_pt>30) && DY_drll>0.3"
#  SS_filters="ttc_region==3 && ttc_mll>60 && ttc_mll<120 && (ttc_l1_pt>30 || ttc_l2_pt>30) && ttc_drll>0.3"
  OS_filters="DY_region==3 && DY_z_mass>%f && DY_z_mass<%f && MET_T1Smear_pt < 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter"%((OS_Zmass-OS_width),(OS_Zmass+OS_width)) 
  SS_filters="ttc_region==3 && ttc_mll>%f && ttc_mll<%f && MET_T1Smear_pt < 50 && n_tight_jet <3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
#  OS_filters="DY_region==3 && DY_z_mass>60 && MET_T1Smear_pt < 50 && n_tight_jet < 3"
#  SS_filters="ttc_region==3 && ttc_mll>60 && MET_T1Smear_pt < 50 && n_tight_jet < 3"
 
  for process in process_list:
    print("Processing %s..."%process['name'])
    #if not ('t_tch' in process['name']):
    #  continue
    df_OS = ROOT.RDataFrame("Events", process['list'])
    df_SS = ROOT.RDataFrame("Events", process['list'])
    df_OS_tree = select_MC_event(df_OS, add_trigger_SF, OS_filters, 1, sigma, process['name'])
    df_SS_tree = select_MC_event(df_SS, add_trigger_SF, SS_filters, 0, sigma, process['name'])
    MC_histos_tmp = []
    for i in OS_hists_name:
      df_histo = df_OS_tree.Histo1D((process['name']+'_OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      MC_histos_tmp.append(df_histo.GetValue().Clone())
    for i in SS_hists_name:
      df_histo = df_SS_tree.Histo1D((process['name']+'_SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      MC_histos_tmp.append(df_histo.GetValue().Clone())
    process['hist'] = MC_histos_tmp

  OS_filters="DY_region==3 && DY_z_mass>%f && DY_z_mass<%f && MET_T1_pt < 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
  SS_filters="ttc_region==3 && ttc_mll>%f && ttc_mll<%f && MET_T1_pt < 50 && n_tight_jet <3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
#  OS_filters="DY_region==3 && DY_z_mass>60 && MET_T1_pt < 50 && n_tight_jet < 3"
#  SS_filters="ttc_region==3 && ttc_mll>60 && MET_T1_pt < 50 && n_tight_jet < 3"


  df_DoubleEle_tree = ROOT.RDataFrame("Events", doubleEle_names)
  df_OS_DoubleEle = df_DoubleEle_tree.Filter(OS_filters)
  df_SS_DoubleEle = df_DoubleEle_tree.Filter(SS_filters)
  df_OS_DoubleEle = df_OS_DoubleEle.Define("DY_kinematic_region","kinematic(DY_l1_pt, DY_l2_pt, DY_l1_eta, DY_l2_eta)")
  df_SS_DoubleEle = df_SS_DoubleEle.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")

  df_OS_DoubleEle_trigger = for_diele_trigger(df_OS_DoubleEle)
  df_SS_DoubleEle_trigger = for_diele_trigger(df_SS_DoubleEle)
  df_DoubleEle_histos=[]
  for i in OS_hists_name:
    df_DoubleEle_histo = df_OS_DoubleEle_trigger.Histo1D(('DoubleEle_OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
    df_DoubleEle_histos.append(df_DoubleEle_histo)
  for i in SS_hists_name:
    df_DoubleEle_histo = df_SS_DoubleEle_trigger.Histo1D(('DoubleEle_SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
    df_DoubleEle_histos.append(df_DoubleEle_histo)

  df_SingleEle_tree = ROOT.RDataFrame("Events", singleEle_names)
  df_OS_SingleEle = df_SingleEle_tree.Filter(OS_filters)
  df_SS_SingleEle = df_SingleEle_tree.Filter(SS_filters)
  df_OS_SingleEle = df_OS_SingleEle.Define("DY_kinematic_region","kinematic(DY_l1_pt, DY_l2_pt, DY_l1_eta, DY_l2_eta)")
  df_SS_SingleEle = df_SS_SingleEle.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
  df_OS_SingleEle_trigger = for_singleele_trigger_eechannel(df_OS_SingleEle)
  df_SS_SingleEle_trigger = for_singleele_trigger_eechannel(df_SS_SingleEle)
  df_SingleEle_histos=[]
  for i in OS_hists_name:
    df_SingleEle_histo = df_OS_SingleEle_trigger.Histo1D(('SingleEle_OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
    df_SingleEle_histos.append(df_SingleEle_histo)
  for i in SS_hists_name:
    df_SingleEle_histo = df_SS_SingleEle_trigger.Histo1D(('SingleEle_SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
    df_SingleEle_histos.append(df_SingleEle_histo)


  fout = ROOT.TFile.Open("ChargeFlip_comparison_sigma_"+str(sigma)+".root","RECREATE")
  fout.cd()
  for ij in range(0,len(histos_bins)):
# ROOT version 6.14 don;t have function "ROOT.RDF.RunGraphs"
#  ROOT.RDF.RunGraphs({df_ZZG_histo, df_ZZ_histo, df_ggZZ_4e_histo,df_ggZZ_4mu_histo, df_ggZZ_4tau_histo, df_ggZZ_2e2mu_histo,df_ggZZ_2e2tau_histo, df_ggZZ_2mu2tau_histo, df_TTZ_histo,df_TTG_histo, df_WWZ_histo, df_WZG_histo,df_WZZ_histo, df_ZZZ_histo, df_WZTo3L_histo,df_WZTo2L_histo, df_ZG_histo})
    h_DoubleEle = df_DoubleEle_histos[ij].GetValue().Clone()
    h_SingleEle = df_SingleEle_histos[ij].GetValue().Clone()

    for process in process_list:
      h = process['hist'][ij].Clone()
      h.Scale(process['xs']/process['ev'])
      histos.append(h.Clone()) 

    histos.append(h_DoubleEle.Clone()) 
    histos.append(h_SingleEle.Clone())

    for i in range(0,27):
      histos[i]=overunder_flowbin(histos[i])
      histos[i].Write()
   
#    c1 = plot_DYregion.draw_plots(histos, 1, hists_name[ij], 0)
    del histos[:]

if __name__ == "__main__":
  start = time.time()
  start1 = time.clock() 
  TTC_Analysis(0)
  end = time.time()
  end1 = time.clock()
  print "wall time:", end-start
  print "process time:", end1-start1
