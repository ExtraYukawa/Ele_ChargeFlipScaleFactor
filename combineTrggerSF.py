import ROOT

sourceDir = "/eos/user/z/zhenggan/ExtraYukawa/TriggerSF/"

def get_hist(era,channel,subtitle,fout):
  infile = sourceDir + "year" + era + "/" + channel + "/files/SF_" + subtitle + ".root"
  fin = ROOT.TFile.Open(infile,"READ")
  h = fin.Get(subtitle)
  h_copy = ROOT.TH2D()
  h.Copy(h_copy)
  fout.cd()
  channelname = ''
  subtitlename = ''
  if(channel == "DoubleMuon"):
     channelname = "mumu_"
  elif(channel == "DoubleElectron"):
     channelname = 'ee_'
  elif(channel == "ElectronMuon"):
     channelname = 'emu_'
  if(subtitle == 'l1pteta'):
     subtitlename = 'lep1pteta'
  elif(subtitle == 'l2pteta'):
     subtitlename = 'lep2pteta'
  h.Write("h2D_SF_"+channelname+subtitlename)
  fin.Close()
def produce_combine_file(era):
  channels = ["DoubleMuon", "DoubleElectron", "ElectronMuon"]
  subtitles = ["l1pteta", "l2pteta"]
  fout = ROOT.TFile.Open("data/TriggerSF_"+era+"UL.root","RECREATE")
  for channel in channels:
    for subtitle in subtitles:
       get_hist(era,channel,subtitle,fout)
  fout.Close()
if __name__ == "__main__":
  eras = ["2016apv","2016postapv","2017","2018"]
  for era in eras:
    produce_combine_file(era)

