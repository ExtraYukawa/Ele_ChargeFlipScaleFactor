import ROOT
import time
import os
import math
from math import sqrt
import plot_DYregion 
from common import overunder_flowbin, get_mcEventnumber, all_trigger, for_diele_trigger, for_singleele_trigger_eechannel, kinematic, select_MC_event, for_egamma_trigger_eechannel, data_trigger
import numpy as np
import optparse
import json
from collections import OrderedDict

# the EnableImplicitMT option should only use in cluster, at lxplus, it will make the code slower(my experience)
#ROOT.ROOT.EnableImplicitMT()

def analysis(region, era, isMC, add_trigger_SF, fin, xs, cfsf, shift, start, end, label, process, subprocess):

#############
##  BASIC  ##
#############

  if(era == '2016'):
    era = '2016postapv' # convention
  else:
    pass
 
  TTC_header_path = os.path.join("scripts/TTC_" + era + ".h")
  ROOT.gInterpreter.Declare('#include "{}"'.format(TTC_header_path))

  form = "MC" if isMC else "Data"
  channel = "DoubleElectron" 
  subera = fin.replace('.root','')
  if ('2016' in era):
    subera = subera.replace('_','')
  dataset  = None
  if 'EGamma' in subera:
    subera = subera.replace('EGamma','')
    dataset = 'EGamma'
  elif "DoubleEG" in subera:
    subera = subera.replace("DoubleEG","")
    dataset = 'DoubleEG'
  elif "SingleEG" in subera:
    subera = subera.replace("SingleEG","")
    dataset = 'SingleEG'

  # Load Trigger json file
  jsonfile = open(os.path.join('data/DiLeptonTriggers_%s.json'%era))
  trig_list = json.load(jsonfile, encoding='utf-8')
  jsonfile.close()
  # Load Trigger command json file
  jsonfile = open(os.path.join('data/Trigger_command_%s.json'%era))
  trig_command_list = json.load(jsonfile, encoding='utf-8',object_pairs_hook=OrderedDict)
  jsonfile.close()
  # Load MET filters json file
  jsonfile = open(os.path.join('data/%s_MET_Filters.json'%era))
  MET_list = json.load(jsonfile, encoding='utf-8',object_pairs_hook=OrderedDict)
  jsonfile.close()
  if isMC:
    MET_filters = MET_list["MC"][process][subprocess]  
  else:
    MET_filters = MET_list["Data"]

  if(era == "2017"):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/TTC_version9/'
  elif(era == '2018'):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/2018/'
  elif(era == '2016postapv'):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/2016postapvMerged/'
  elif(era == '2016apv'):
    path = '/eos/cms/store/group/phys_top/ExtraYukawa/2016apvMerged/'
  else:
    print("You must select 2016postapv/2016apv/2017/2018. Now set default to 2017.")
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
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1Smear_pt > 50  && n_tight_jet < 3 && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.5 && abs(OPS_l2_eta)<2.5 && nHad_tau==0 && OPS_drll > 0.3"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1Smear_pt > 50 && n_tight_jet < 3 && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.5 && abs(ttc_l2_eta)<2.5 && nHad_tau==0 && ttc_drll > 0.3"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
    else:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1_pt > 50 && n_tight_jet < 3 && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.5 && abs(OPS_l2_eta)<2.5 && nHad_tau==0 && OPS_drll > 0.3"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1_pt > 50 && n_tight_jet < 3 && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.5 && abs(ttc_l2_eta) < 2.5 && nHad_tau==0 && ttc_drll > 0.3"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))

  else:
    if isMC:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1Smear_pt < 35 && n_tight_jet < 3 && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.5 && abs(OPS_l2_eta)<2.5 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width)) 
    #OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>60 && OPS_z_mass<120 && OPS_l1_pt>30 && OPS_l2_pt>20 && OPS_drll>0.3 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && nHad_tau==0&& abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4"
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1Smear_pt < 35 && n_tight_jet <3 && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.5 && abs(ttc_l2_eta)<2.5&& nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))
    else:
      OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>%f && OPS_z_mass<%f && MET_T1_pt < 35 && n_tight_jet < 3 && OPS_l1_pt>20 && OPS_l2_pt>20 && abs(OPS_l1_eta)<2.5 && abs(OPS_l2_eta)<2.5 && nHad_tau==0"%((OS_Zmass-OS_width),(OS_Zmass+OS_width))
#    OS_filters="OPS_region==3 && OPS_2P0F && OPS_z_mass>60 && OPS_z_mass<120 && OPS_l1_pt>30 && OPS_l2_pt>20 && OPS_drll>0.3 && Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && Flag_eeBadScFilter && Flag_ecalBadCalibFilter && nHad_tau==0 && abs(OPS_l1_eta)<2.4 && abs(OPS_l2_eta)<2.4"
      SS_filters="ttc_region==3 && ttc_2P0F && ttc_mll>%f && ttc_mll<%f && MET_T1_pt < 35 && n_tight_jet <3  && ttc_l1_pt>20 && ttc_l2_pt>20 && abs(ttc_l1_eta)<2.5 && abs(ttc_l2_eta)<2.5 && nHad_tau==0"%((SS_Zmass-SS_width),(SS_Zmass+SS_width))

  for MET_filter in MET_filters:
    OS_filters += " && %s"%MET_filter
    SS_filters += " && %s"%MET_filter

  OS_filters = str(OS_filters)
  SS_filters = str(SS_filters)
 
  print(OS_filters)
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
    OS_filters += " && !(OPS_isHEM==1)"
    SS_filters += " && !(ttc_isHEM==1)"


#  if isMC:
#    df_OS.Define("OPS_MET", "MET_T1Smear_pt")

#  else:
#    df_OS.Define("OPS_MET", "MET_T1_pt")
#    df_SS.Define("ttc_MET", "MET_T1_pt")

#  OS_hists_name.append("OPS_MET")
#  SS_hists_name.append("ttc_MET")
#  histos_bins["OPS_MET"] = [20,0,200]
#  histos_bins["ttc_MET"] = [20,0,200]


  # Basic flatten process



  if isMC:
    df_OS_tree = select_MC_event(df_OS, add_trigger_SF, cfsf, shift, OS_filters, 1, era, fin,channel,trig_command_list)
    df_SS_tree = select_MC_event(df_SS, add_trigger_SF, cfsf, shift, SS_filters, 0, era, fin,channel,trig_command_list)

    for i in OS_hists_name:
      print(i)
      df_histo = df_OS_tree.Histo1D(('OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      histos.append(df_histo.GetValue().Clone())
    for i in SS_hists_name:
      df_histo = df_SS_tree.Histo1D(('SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i,'genweight')
      histos.append(df_histo.GetValue().Clone())

    df_histo = df_OS_tree.Histo1D(('OS_MET','',20,0,200),'MET_T1Smear_pt','genweight')
    histos.append(df_histo.GetValue().Clone())
    df_histo = df_SS_tree.Histo1D(('SS_MET','',20,0,200),'MET_T1Smear_pt','genweight')
    histos.append(df_histo.GetValue().Clone())


  else:

  ## Data ##
    df_OS_tree = df_OS.Filter(OS_filters)
    df_SS_tree = df_SS.Filter(SS_filters)
    df_OS_tree = df_OS_tree.Define("OPS_kinematic_region","kinematic(OPS_l1_pt, OPS_l2_pt, OPS_l1_eta, OPS_l2_eta)")
    df_SS_tree = df_SS_tree.Define("ttc_kinematic_region","kinematic(ttc_l1_pt, ttc_l2_pt, ttc_l1_eta, ttc_l2_eta)")
    df_OS_tree = data_trigger(df_OS_tree, channel, subera, dataset, trig_list, trig_command_list)
    df_SS_tree = data_trigger(df_SS_tree, channel, subera, dataset, trig_list, trig_command_list)

    for i in OS_hists_name:
      df_data_histo = df_OS_tree.Histo1D(('OS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
      histos.append(df_data_histo)
    for i in SS_hists_name:
      df_data_histo = df_SS_tree.Histo1D(('SS_'+i,'',histos_bins[i][0],histos_bins[i][1],histos_bins[i][2]), i)
      histos.append(df_data_histo)

    df_histo = df_OS_tree.Histo1D(('OS_MET','',20,0,200),'MET_T1_pt')
    histos.append(df_histo.GetValue().Clone())
    df_histo = df_SS_tree.Histo1D(('SS_MET','',20,0,200),'MET_T1_pt')
    histos.append(df_histo.GetValue().Clone())


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
  print(fout)
  fout.cd()
  for ij in range(0,len(histos)):
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
  parser.add_option('-i','--fin', dest='fin', help='input file', default = 'ttH.root', type='string') 
  parser.add_option('-x','--xsec',dest='xsec',help='x section',  default = 1,    type=float)
  parser.add_option('-t','--trig',dest='trig',help='Add trigger SF or not', default = 0, type=int)
  parser.add_option('-c','--cfsf',dest='cfsf',help='Add chargeflip SF or not', default =1, type=int)
  parser.add_option('-s','--shift',dest='shift',help='shift chargeflip SF', default =0, type=int)
  parser.add_option('-f','--from',dest='start',help='start index', default = 0, type = int)
  parser.add_option('-n','--to'  ,dest='end',help='end index', default = 100, type=int)
  parser.add_option('-l','--label',dest='label',help='label',  default = 0, type=int)
  parser.add_option('-r','--region',dest='region',help='plot region',default = 'nominal', type='string')
  parser.add_option('-p','--process',dest='process',help='main process',default = 'ttXorXX')
  parser.add_option('-u','--subprocess',dest='subprocess',help='sub process',default = 'ttH')
  (args,opt) = parser.parse_args()
  analysis(args.region, args.era, args.isMC, args.trig, args.fin, args.xsec, args.cfsf, args.shift, args.start, args.end, args.label, args.process,args.subprocess)
