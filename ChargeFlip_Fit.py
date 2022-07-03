import ROOT
import numpy as np
import math
import ctypes

def fit(era, h_sig_OS, h_sig_SS, h_back_OS, h_back_SS, h_data_OS, h_data_SS, pt_region, eta_region,useLikelihood, subMC, draw):

  print("sig_OS: %f"%h_sig_OS.Integral())
  print("data_OS: %f"%h_data_OS.Integral())
  print("sig_SS: %f"%h_sig_SS.Integral())
  print("data_SS: %f"%h_data_SS.Integral())
  print("back_OS: %f"%h_back_OS.Integral())
  print("back_SS: %f"%h_back_SS.Integral())

  plotdir = '/eos/user/t/tihsu/plot/Ele_chargeflip_sf_eratrig/'
  method = 'chi2'
  if(useLikelihood):
    method = 'MLE'
  pt_region  = np.array(pt_region)
  eta_region = np.array(eta_region)
  pt_bins = len(pt_region)-1
  eta_bins = len(eta_region)-1
  channel = ['MC_os','MC_ss','data_os','data_ss']
  Number = dict()
  Pij = dict()
  Vij = dict()
  P_fit = dict()

  sig_OS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  sig_SS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  back_OS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  back_SS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  data_OS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
  data_SS = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))


  for i in range(pt_bins):
   for j in range(eta_bins):
     for ii in range(pt_bins):
       for jj in range(eta_bins):
         index = jj+ii*eta_bins+j*eta_bins*pt_bins+i*eta_bins*pt_bins*eta_bins
         sig_OS[i][j][ii][jj] = h_sig_OS.GetBinContent(index+2)
         sig_SS[i][j][ii][jj] = h_sig_SS.GetBinContent(index+2)
         back_OS[i][j][ii][jj] = h_back_OS.GetBinContent(index+2)
         back_SS[i][j][ii][jj] = h_back_SS.GetBinContent(index+2)
         data_OS[i][j][ii][jj] = h_data_OS.GetBinContent(index+2)
         data_SS[i][j][ii][jj] = h_data_SS.GetBinContent(index+2)
  
  Number['MC_os'] = sig_OS
  Number['MC_ss'] = sig_SS
  Number['data_os'] = data_OS - (back_OS * float(subMC))
  Number['data_ss'] = data_SS - (back_SS * float(subMC))

  w = 600
  he = 600
  c = ROOT.TCanvas('c','c',10,10,w,he)
  ROOT.gStyle.SetOptStat("kFALSE")
  ROOT.gStyle.SetPaintTextFormat(".2e")
  ROOT.gStyle.SetPalette(69);
  ROOT.gStyle.SetCanvasBorderSize(0)
  ROOT.gStyle.SetCanvasBorderMode(0)  
  ROOT.gStyle.SetFrameBorderMode(0)
  ROOT.gROOT.SetStyle("Plain");
  c.SetRightMargin(0.15);
  c.SetTopMargin(0.15);
  c.Update()


  for h in Number:
    for i in range(pt_bins):
      for j in range(eta_bins):
        for ii in range(pt_bins):
          for jj in range(eta_bins):
            if(ii>i or jj>j):
              Number[h][i][j][ii][jj]+=Number[h][ii][jj][i][j]
              Number[h][ii][jj][i][j] = 0
  if draw:
    Fout = ROOT.TFile.Open('data/ChargeFlipProbability_' + era + '_'+method+'.root',"Recreate")
  Data_MC = ['data','MC']
  for h in Data_MC:
    Pij[h] = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
    Vij[h] = np.zeros((pt_bins,eta_bins,pt_bins,eta_bins))
    P_fit[h] = np.zeros((pt_bins,eta_bins))
    for i in range(pt_bins):
      for j in range(eta_bins):
        for ii in range(pt_bins):
          for jj in range(eta_bins):
            N_ss = Number[h+'_ss'][i][j][ii][jj]
            N_T = N_ss + Number[h+'_os'][i][j][ii][jj]
#            if(N_ss < 0.): N_ss = 0.
            if(N_ss>=0. and N_T>0. and N_T>=N_ss):
              Pij[h][i][j][ii][jj] = N_ss/N_T
              Vij[h][i][j][ii][jj] = N_ss/(N_T**2)+(N_ss**2)/(N_T**3)-2*(N_ss**1.5)/(N_T**2.5)
    gMinuit = ROOT.TMinuit(pt_bins*eta_bins)
    
    def fcn(npar, gin, f,par,iflag):
      chi2 = 0.
      Likelihood = 0.
      v = 0
      for i in range(pt_bins):
        for j in range(eta_bins):
          for ii in range(pt_bins):
            for jj in range(eta_bins):
              P1 = par[i*eta_bins+j]
              P2 = par[ii*eta_bins+jj]
              eP = P1+P2-2.*P1*P2
              N = int(round(Number[h+'_os'][i][j][ii][jj]+Number[h+'_ss'][i][j][ii][jj]))
              Nsc = int(round(Number[h+'_ss'][i][j][ii][jj]))
              if ((not Vij[h][i][j][ii][jj]==0.) and Number[h+'_os'][i][j][ii][jj]>0. and Number[h+'_ss'][i][j][ii][jj]>=1.):
                chi2 += ((Pij[h][i][j][ii][jj]-(P1+P2-2.*P1*P2))**2)/(Vij[h][i][j][ii][jj])
                Likelihood += Nsc*np.math.log((N*(eP)))-N*(eP)-np.math.log(np.math.factorial(Nsc))
                v+=1
      f[0] = chi2/float(v-9)
      if(useLikelihood):
        f[0] = -1.*Likelihood

    gMinuit.SetFCN(fcn)
    for i in range(pt_bins):
      for j in range(eta_bins):
        init_val = 0.0001
        if Number[h+'_os'][i][j][i][j]>5000:
          init_val = ((1.-(1.-2.*Pij[h][i][j][i][j])**0.5)/2.)
        gMinuit.DefineParameter(i*eta_bins+j,"P"+str(i)+str(j),init_val,0.00000001,0.,0.1)
    gMinuit.Command("Minuit2")

    w = 600
    he = 600
    c = ROOT.TCanvas(h+'c','c',10,10,w,he)
    ROOT.gStyle.SetOptStat("kFALSE")
    ROOT.gStyle.SetPaintTextFormat(".2e")
    ROOT.gStyle.SetPalette(69);
    ROOT.gStyle.SetCanvasBorderSize(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gROOT.SetStyle("Plain");
    c.SetRightMargin(0.15);
    c.SetTopMargin(0.15);
    c.Update()


    f_value = ROOT.Double(0.)
    fedm = ROOT.Double(0.)
    errdef = ROOT.Double(0.)
    npari = ctypes.c_int(0)
    nparx = ctypes.c_int(0)
    istat = ctypes.c_int(0)
    gMinuit.mnstat(f_value,fedm,errdef,npari,nparx,istat)
    print("Goodness of fit %f"%f_value)

    npar = gMinuit.GetNumPars()
    matrix = ROOT.TMatrixD(npar,npar);
    gMinuit.mnemat(matrix.GetMatrixArray(),npar)
    matrix.Print();

    h_chargeflip = ROOT.TH2D(h+"_CFRate",";P_{T}[GeV] ; |\eta|};",pt_bins,pt_region,eta_bins,eta_region)
    h_chargeflip_cov = ROOT.TH2D(h+"_CovMatrix",";;",npar,0,9,npar,0,9)
    for i in range(npar):
      for j in range(npar):
        h_chargeflip_cov.SetBinContent(i+1,j+1,matrix[i][j])
    h_chargeflip_cov.Draw("COLZTEXT e")
    if draw:
      c.SaveAs(plotdir + 'CovMatrix_'+h+'_'+method+ '_' + era + '.png')

    result = [[0. for j in range(eta_bins)] for i in range(pt_bins)]
    error = [[0. for j in range(eta_bins)] for i in range(pt_bins)]

    for i in range(pt_bins):
      for j in range(eta_bins):
        result = ROOT.double(0.)
        error = ROOT.double(0.)
        error_plus = ROOT.double(0.)
        error_minus= ROOT.double(0.)
        eparab     = ROOT.double(0.)
        gcc        = ROOT.double(0.)
        gMinuit.GetParameter(i*eta_bins+j,result,error)
        gMinuit.mnerrs(i*eta_bins+j,error_plus,error_minus,eparab,gcc)
        P_fit[h][i][j] = result
        h_chargeflip.SetBinContent(i+1,j+1,result)
        h_chargeflip.SetBinError(i+1,j+1,error)
    if h=='data':
      h_data = h_chargeflip.Clone()
      h_data_cov = h_chargeflip_cov.Clone()
    else:
      h_MC = h_chargeflip.Clone()
      h_MC_cov = h_chargeflip_cov.Clone()
    c.SetLogx()
    c.SetLogz()
    ROOT.gStyle.SetOptStat("kFALSE")
    ROOT.gStyle.SetPaintTextFormat(".2e")
    ROOT.gROOT.SetStyle("Plain");
    h_chargeflip.Draw('COLZTEXT e')
    c.Update()
    if draw:
      c.SaveAs(plotdir+h+'_CFRate_'+method+'_' + era + '.png')
      c.SaveAs(plotdir+h+'_CFRate_'+method+'_' + era + '.pdf')

# Fout = ROOT.TFile.Open(plotdir+'ChargeFlipProbability_2017_'+method+'.root',"Recreate")
  h_data.Draw('COLZTEXT e')
  c.Update()
  if draw:
    c.SaveAs(plotdir+'Data_CFRate_'+method+'_' + era + '.png')
    c.SaveAs(plotdir+'Data_CFRate_'+method+'_' + era + '.pdf')
    h_data.Write()
    h_data_cov.Write()
    h_MC.Write()
    h_MC_cov.Write()

    c.SetLogz(0)
    h_SF = h_data.Clone()
    h_SF.Divide(h_SF,h_MC)
    h_SF.Draw('COLZTEXT e')
    c.Update()
    c.SaveAs(plotdir+'CFRate_'+method+'_MDRatio_' + era + '.png')
    c.SaveAs(plotdir+'CFRate_'+method+'_MDRatio_' + era + '.pdf')

  SF = np.ones((pt_bins,eta_bins))
  for i in range(pt_bins):
    for j in range(eta_bins):
      SF[i][j] = P_fit['data'][i][j]/P_fit['MC'][i][j]

#  print(P_fit['MC'])
#  print(SF)
  h_data.SetDirectory(0)
  h_MC.SetDirectory(0)
  h_data_cov.SetDirectory(0)
  h_MC_cov.SetDirectory(0)
  if draw:
    Fout.Close()
  return h_data, h_MC, h_data_cov, h_MC_cov
