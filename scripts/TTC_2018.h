#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "ROOT/RVec.hxx"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "TH1D.h"

TString era = "2018";
TFile*f=TFile::Open("data/TriggerSF_"+era+"UL.root");
TH2D*h1_ee=(TH2D*)f->Get("h2D_SF_ee_SF_l1l2pt");
TH2D*h1_mm=(TH2D*)f->Get("h2D_SF_mumu_SF_l1l2pt");
TH2D*h1_em=(TH2D*)f->Get("h2D_SF_emu_SF_l1l2pt");

TFile*f_cf=TFile::Open("data/ChargeFlipSF_" + era + "_MLE.root");
TH1D*h_OS=(TH1D*)f_cf->Get("OS_ChargeFlip_SF");
TH1D*h_SS  =(TH1D*)f_cf->Get("SS_ChargeFlip_SF_AllUnc");

TFile*f_cfregion=TFile::Open("data/ChargeFlipProbability_" + era + "_MLE.root");
TH2D*h_data = (TH2D*) f_cfregion->Get("data_CFRate");

TFile*fele=TFile::Open("data/EleIDSF_" + era + ".root");
TH2D*h_eleSF=(TH2D*)fele->Get("EleIDSF");

//srand(12345);

float eleID_sf_ee(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>500) l1_pt=499.;
	if(l2_pt>500) l2_pt=499.;
	float sf_l1=h_eleSF->GetBinContent(h_eleSF->FindBin(l1_pt,fabs(l1_eta)));
	float sf_l2=h_eleSF->GetBinContent(h_eleSF->FindBin(l2_pt,fabs(l2_eta)));
	return sf_l1*sf_l2;
}

float trigger_sf_ee(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_ee->GetBinContent(h1_ee->FindBin(l1_pt,l2_pt));
	return sf_l1;
}

float trigger_sf_mm(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_mm->GetBinContent(h1_mm->FindBin(l1_pt,l2_pt));
	return sf_l1;
}

float trigger_sf_em(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
	if(l1_pt>200) l1_pt=199;
	if(l2_pt>200) l2_pt=199;
	float sf_l1=h1_em->GetBinContent(h1_em->FindBin(l1_pt,l2_pt));
	return sf_l1;
}

int kinematic(float l1_pt, float l2_pt, float l1_eta, float l2_eta){
  std::vector<Float_t> pt_region  = {20., 40., 60., 100., 100000000000.};
  std::vector<Float_t> eta_region = {0.,  0.8, 1.479, 2.5};
  int pt_bins  = pt_region.size() - 1;
  int eta_bins = eta_region.size() - 1;
  int l1_pt_index = -1;
  int l2_pt_index = -1;
  int l1_eta_index = -1;
  int l2_eta_index = -1;
  for(int i = 0; i < pt_bins; i++){
    if(pt_region[i] <= l1_pt && l1_pt < pt_region[i+1]) l1_pt_index = i;
    if(pt_region[i] <= l2_pt && l2_pt < pt_region[i+1]) l2_pt_index = i;
  }
  for(int i = 0; i < eta_bins; i++){
    if(eta_region[i] <= fabs(l1_eta) && fabs(l1_eta) < eta_region[i+1]) l1_eta_index = i;
    if(eta_region[i] <= fabs(l2_eta) && fabs(l2_eta) < eta_region[i+1]) l2_eta_index = i;
  }
  if(l1_pt_index <0 || l2_pt_index <0 || l1_eta_index <0 || l2_eta_index<0) return -1;
  return l2_eta_index + eta_bins*(l2_pt_index) + eta_bins*pt_bins*(l1_eta_index) + eta_bins*pt_bins*eta_bins*(l1_pt_index);
}
     

float chargeflip_sf(int kinematic_region, int charge1, int charge2, float sigma){
    int index = kinematic_region+1;
    float sf = 1.0;
    if((charge1*charge2)<0){
      float OS_SF = h_OS->GetBinContent(index);
      float OS_sigma = h_OS->GetBinError(index);
      sf = OS_SF - sigma*OS_sigma;
    }
    else if((charge1*charge2)>0){
      float SS_SF = h_SS->GetBinContent(index);
      float SS_sigma = h_SS->GetBinError(index);
      sf = SS_SF + sigma*SS_sigma;
    }
    if(sf<0.) sf=0.;
    return sf;
}

int isHEM(int run, ROOT::VecOps::RVec<float> Jet_pt, ROOT::VecOps::RVec<Float_t> Jet_eta, ROOT::VecOps::RVec<Float_t> Jet_phi, int isMC){

  int nJet = Jet_pt.size();
  int HEM = 0; 

  // HEM: 0: not in HEM region 1:in HEM region and era 2:only in HEM region but not era 3: not in HEM region but in era
  for(int i = 0; i < nJet; i++){
    if((Jet_pt[i] > 15) && (Jet_eta[i] < -1.3) && (-3.2 < Jet_eta[i]) && (-1.57 < Jet_phi[i]) && (Jet_phi[i] < -0.87)){
      HEM = 2;
      continue;
    }
  }
  if (not isMC){
    if(HEM == 2 && run >= 319077){
      HEM = 1;
    }
    else if(HEM == 0 && run >= 319077){
      HEM = 3;
    }
  }
  else{
    float probability = (float)rand()/(float)RAND_MAX;
    if(probability < 0.648189){
      if(HEM == 2){
        HEM = 1; 
      }
      else HEM = 3;
    }
  }  
  return HEM;
}
bool Triggers(int run, bool triggers, std::vector<int> vec){
    if(!triggers || vec.at(0) == -1 ){
        return triggers;
    }
    for(auto v : vec){
        if(run == v){
            return triggers;
        }
    }
    return false;
}

