import ROOT
import os
import math
import numpy as np


def SF_produce(era, h_data, h_MC, h_data_cov, h_MC_cov, draw, tag):
  # Obtain 2D histogram from root file.
  plotdir = '/eos/user/t/tihsu/plot/Ele_chargeflip_sf/'
  pt_bins = h_data.GetXaxis().GetNbins()
  eta_bins = h_data.GetYaxis().GetNbins()
  print(pt_bins)
  # Calculate SF in respect to SS and OS event.
  h_SF_OS = ROOT.TH1D("OS_ChargeFlip_SF" + tag, ";;", (pt_bins*eta_bins)**2+10,0,(pt_bins*eta_bins)**2+10)
  h_SF_SS = ROOT.TH1D("SS_ChargeFlip_SF" + tag, ";;", (pt_bins*eta_bins)**2+10,0,(pt_bins*eta_bins)**2+10)

  for i in range(pt_bins):
    for j in range(eta_bins):
      for ii in range(pt_bins):
        for jj in range(eta_bins):
          index = jj + ii*eta_bins + j*eta_bins*pt_bins + i*eta_bins*pt_bins*eta_bins
          index_1 = j + i*eta_bins
          index_2 = jj+ ii*eta_bins

          P1_data = h_data.GetBinContent(i+1,j+1)
          P2_data = h_data.GetBinContent(ii+1,jj+1)
          data_error2_p1 = h_data_cov.GetBinContent(index_1,index_1)
          data_error2_p2 = h_data_cov.GetBinContent(index_2,index_2)
          data_error2_p1p2 = h_data_cov.GetBinContent(index_1, index_2)

          P1_MC   = h_MC.GetBinContent(i+1,j+1)
          P2_MC   = h_MC.GetBinContent(ii+1,jj+1)
          MC_error2_p1 = h_MC_cov.GetBinContent(index_1,index_1)
          MC_error2_p2 = h_MC_cov.GetBinContent(index_2,index_2)
          MC_error2_p1p2 = h_MC_cov.GetBinContent(index_1, index_2)

          P_data  = P1_data + P2_data - 2.*P1_data*P2_data
          P_data_err2 = ((1.-2.*P2_data)**2)*data_error2_p1 + ((1.-2.*P1_data)**2)*data_error2_p2 + 2.*(1.-2.*P1_data)*(1.-2.*P2_data)*data_error2_p1p2

          P_MC  = P1_MC + P2_MC - 2.*P1_MC*P2_MC
          P_MC_err2 = ((1.-2.*P2_MC)**2)*MC_error2_p1 + ((1.-2.*P1_MC)**2)*MC_error2_p2 + 2.*(1.-2.*P1_MC)*(1.-2.*P2_MC)*MC_error2_p1p2

          SS_SF = P_data/P_MC
          SS_error = ((P_data_err2/(P_data**2)+P_MC_err2/(P_MC**2))**0.5)*SS_SF

          #print(P_data)
          #print(P_MC)
          OS_SF = (1.-P_data)/(1.-P_MC)
          print(OS_SF)
          OS_error = ((P_data_err2/((1.-P_data)**2)+P_MC_err2/((1.-P_MC)**2))**0.5)*OS_SF
          if(1):
            h_SF_OS.SetBinContent(index+1,OS_SF);
            h_SF_OS.SetBinError(index+1,OS_error);
 
            h_SF_SS.SetBinContent(index+1,SS_SF);
            h_SF_SS.SetBinError(index+1,SS_error);
          
  if draw:
    c = ROOT.TCanvas()
    h_SF_OS.GetYaxis().SetRangeUser(0.95,1.05);
    h_SF_OS.Draw('pe')
    c.SaveAs(plotdir + '/SF_OS_' + era + '.png')
    c.SaveAs(plotdir + '/SF_OS_' + era + '.pdf')
    h_SF_SS.Draw('pe')
    c.SaveAs(plotdir + '/SF_SS_' + era + '.png')
    c.SaveAs(plotdir + '/SF_SS_' + era + '.pdf')
    fout = ROOT.TFile.Open("data/ChargeFlipSF_%s_MLE.root"%era,"RECREATE")
    fout.cd()
    h_SF_OS.Write()
    h_SF_SS.Write()
    h_SF_SS_sys = h_SF_SS.Clone()
    h_SF_SS_sys.SetName("SS_ChargeFlip_SF_sys"+tag)
    h_SF_SS_sys.Write()
    h_SF_OS.SetDirectory(0)
    h_SF_SS.SetDirectory(0)
    fout.Close()
  return h_SF_OS, h_SF_SS
if __name__ == '__main__':
  Eras = ['2017', '2018']
  for era in Eras:
    fin = ROOT.TFile.Open("data/ChargeFlipProbability_" + era + "_MLE.root")
    h_data = fin.Get("data_CFRate")
    h_MC   = fin.Get("MC_CFRate")
    h_data_cov = fin.Get("data_CovMatrix")
    h_MC_cov   = fin.Get("MC_CovMatrix")
    SF_produce(era, h_data, h_MC, h_data_cov, h_MC_cov, 1, '')




