import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion 
from common import overunder_flowbin, get_mcEventnumber, all_trigger, for_diele_trigger, for_singleele_trigger_eechannel, kinematic, select_MC_event, for_egamma_trigger_eechannel
import numpy as np
import optparse


# the EnableImplicitMT option should only use in cluster, at lxplus, it will make the code slower(my experience)
#ROOT.ROOT.EnableImplicitMT()

def analysis(region, era, isMC, add_trigger_SF, fin, xs, cfsf, shift, start, end, label):

#############
##  BASIC  ##
#############
 
  TTC_header_path = os.path.join("scripts/TTC_" + era + ".h")
  ROOT.gInterpreter.Declare('#include "{}"'.format(TTC_header_path))

  if(era == "2017"):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/TTC_version9/'
  elif(era == '2018'):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/2018/'
  else:
    print("You must select 2017/2018. Now set default to 2017.")
    path='/eos/cms/store/group/phys_top/ExtraYukawa/TTC_version9/'
  
  outdir = 'validation/' if (region == 'validation') else 'flatten/'

############
##  READ  ##
############

  f_list = ROOT.std.vector('string')()
  f_list.push_back(path+fin)

  if isMC:
    ev = get_mcEventnumber(f_list)

#################
##  HISTOGRAM  ##
#################

  #histograms name
  OS_hists_name = ['OPS_l1_pt','OPS_l1_eta','OPS_l1_phi','OPS_l2_pt','OPS_l2_eta','OPS_l2_phi','OPS_z_pt','OPS_z_eta','OPS_z_phi','OPS_z_mass','OPS_kinematic_region']
  SS_hists_name = ['ttc_l1_pt','ttc_l1_eta','ttc_l1_phi','ttc_l2_pt','ttc_l2_eta','ttc_l2_phi','ttc_mll','ttc_kinematic_region']

  #histograms bins [nbins, low edge, high edge]

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
    OS_hists_name[10]:[150,-1,149],
    SS_hists_name[0]:[20,0,200],
    SS_hists_name[1]:[16,-2.4,2.4],
    SS_hists_name[2]:[20,-4,4],
    SS_hists_name[3]:[20,0,100],
    SS_hists_name[4]:[16,-2.4,2.4],
    SS_hists_name[5]:[20,-4,4],
    SS_hists_name[6]:[60,60,120],
    SS_hists_name[7]:[150,-1,149]
  }

  histos = []

#####################
##  Z MASS REGION  ##
#####################

  OS_Zmass = 90.
  SS_Zmass = 90.
  OS_width = 30.
  SS_width = 30.


#########################
##  *   Selection  *   ##
#########################
  if (region == 'validation'):

    if isMC:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1Smear_pt > 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1Smear_pt > 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.4 && abs(ttc_l2_eta)<2.4&& nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
    else:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1_pt > 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1_pt > 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.4 && abs(ttc_l2_eta)<2.4 && nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))

  else:
    if isMC:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1Smear_pt < 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width)) 
    #OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>60 && OPS_z_mass<120 && OPS_l1_pt>30 && OPS_l2_pt>20 && OPS_drll>0.3 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && nHad_tau==0&& abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4"
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1Smear_pt < 50 && n_tight_jet <3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.4 && abs(ttc_l2_eta)<2.4&& nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
    else:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1_pt < 50 && n_tight_jet < 3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
#    OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>60 && OPS_z_mass<120 && OPS_l1_pt>30 && OPS_l2_pt>20 && OPS_drll>0.3 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && nHad_tau==0 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4"
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1_pt < 50 && n_tight_jet <3 && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_BadPFMuonFilter && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.4 && abs(ttc_l2_eta)<2.4 && nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))

###############
##  Flatten  ##
###############

  df_OS_a = ROOT.RDataFrame("Events", f_list)
  df_SS_a = ROOT.RDataFrame("Events", f_list)
  df_OS = df_OS_a.Range(start,end)
  df_SS = df_SS_a.Range(start,end)

  # HEM veto [ 2018 issue ]

  if era == "2018":
    if isMC:
      df_OS = df_OS.Define("OPS_isHEM", "isHEM(run, Jet_pt, Jet_eta, Jet_phi, 1)")
      df_SS = df_SS.Define("ttc_isHEM", "isHEM(run, Jet_pt, Jet_eta, Jet_phi, 1)")
    else:
      df_OS = df_OS.Define("OPS_isHEM", "isHEM(run, Jet_pt, Jet_eta, Jet_phi, 0)")
      df_SS = df_SS.Define("ttc_isHEM", "isHEM(run, Jet_pt, Jet_eta, Jet_phi, 0)")

    OS_hists_name.append("OPS_isHEM")
    SS_hists_name.append("ttc_isHEM")
    histos_bins["OPS_isHEM"] = [4,0,4]
    histos_bins["ttc_isHEM"] = [4,0,4]
    #OS_filters += " && !(OPS_isHEM==1)"
    #SS_filters += " && !(ttc_isHEM==1)"

  # Basic flatten process

  if isMC:
    df_OS_tree = select_MC_event(df_OS, add_trigger_SF, cfsf, shift, OS_filters, 1, era, fin)
    df_SS_tree = select_MC_event(df_SS, add_trigger_SF, cfsf, shift, SS_filters, 0, era, fin)

    for i in OS_hists_name:
      df_histo = df_OS_tree.Histo1D(('OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      histos.append(df_histo.GetValue().Clone())
    for i in SS_hists_name:
      df_histo = df_SS_tree.Histo1D(('SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      histos.append(df_histo.GetValue().Clone())

  else:
    df_OS_tree = df_OS.Filter(OS_filters)
    df_SS_tree = df_SS.Filter(SS_filters)
    df_OS_tree = df_OS_tree.Define("OPS_kinematic_region","kinematic(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_SS_tree = df_SS_tree.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    if "DoubleEG" in fin:
      df_OS_tree = for_diele_trigger(df_OS_tree)
      df_SS_tree = for_diele_trigger(df_SS_tree)
    elif "SingleEG" in fin:
      df_OS_tree = for_singleele_trigger_eechannel(df_OS_tree)
      df_SS_tree = for_singleele_trigger_eechannel(df_SS_tree)
    elif "EGamma" in fin:
      df_OS_tree = for_egamma_trigger_eechannel(df_OS_tree)
      df_SS_tree = for_egamma_trigger_eechannel(df_SS_tree)

    for i in OS_hists_name:
      df_data_histo = df_OS_tree.Histo1D(('OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
      histos.append(df_data_histo)
    for i in SS_hists_name:
      df_data_histo = df_SS_tree.Histo1D(('SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
      histos.append(df_data_histo)
  if len(histos)==0:
    return 0

##############
##  OUTPUT  ##
##############
  
  cfsf_name = "ApplyChargeFlipsf" if cfsf else "NotApplyChargeFlipsf"
  shift_name = "_Nominal/"
  if shift==1:
    shift_name = "_UP/"
  elif shift == -1:
    shift_name = "_DOWN/"
  fout = ROOT.TFile.Open(outdir + "era" + era + "/" + cfsf_name + shift_name + "h" + str(label) + "_" + fin, "RECREATE")
  fout.cd()
  for ij in range(0,len(histos_bins)):
# ROOT version 6.14 don;t have function "ROOT.RDF.RunGraphs"
#  ROOT.RDF.RunGraphs({df_ZZG_histo, df_ZZ_histo, df_ggZZ_4e_histo,df_ggZZ_4mu_histo, df_ggZZ_4tau_histo, df_ggZZ_2e2mu_histo,df_ggZZ_2e2tau_histo, df_ggZZ_2mu2tau_histo, df_TTZ_histo,df_TTG_histo, df_WWZ_histo, df_WZG_histo,df_WZZ_histo, df_ZZZ_histo, df_WZTo3L_histo,df_WZTo2L_histo, df_ZG_histo})
    h = histos[ij].Clone()
    if isMC:
      h.Scale(xs/ev)
    h = overunder_flowbin(h)
    h.Write()
   
#    c1 = plot_DYregion.draw_plots(histos, 1, hists_name[ij], 0)
  del histos[:]

if __name__ == "__main__":

  # Configuration
  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-e','--era', dest='era', help='era: [2017/2018]', default='2017', type='string')
  parser.add_option('-m','--isMC',dest='isMC',help='isMC = [0/1]',  default = 1, type=int)
  parser.add_option('-i','--fin', dest='fin', help='input file', default = 'tttt.root', type='string') 
  parser.add_option('-x','--xsec',dest='xsec',help='x section',  default = 1,    type=float)
  parser.add_option('-t','--trig',dest='trig',help='Add trigger SF or not', default = 0, type=int)
  parser.add_option('-c','--cfsf',dest='cfsf',help='Add chargeflip SF or not', default =1, type=int)
  parser.add_option('-s','--shift',dest='shift',help='shift chargeflip SF', default =0, type=int)
  parser.add_option('-f','--from',dest='start',help='start index', default = 0, type = int)
  parser.add_option('-n','--to'  ,dest='end',help='end index', default = 100, type=int)
  parser.add_option('-l','--label',dest='label',help='label',  default = 0, type=int)
  parser.add_option('-r','--region',dest='region',help='plot region',default = 'nominal', type='string')
  (args,opt) = parser.parse_args()
  analysis(args.region, args.era, args.isMC, args.trig, args.fin, args.xsec, args.cfsf, args.shift, args.start, args.end, args.label)
