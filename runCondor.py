import os
import sys
import optparse
import subprocess
import json
import ROOT
from collections import OrderedDict

def prepare_range(path, fin):
  try:
    f_read = ROOT.TFile.Open(path + fin)
    entries = (f_read.Get('Events')).GetEntriesFast()
    step = 500000
    init = 0
    index = []
    while(init < entries):
      index.append(init)
      init += step
    index.append(int(entries))
    f_read.Close()
    return index
  except:
    print("%s%s fail to process."%(path,fin))
    return None
def prepare_shell(region, era, isMC, trig_sf, cfsf, shift, fin, xs, start, end, label, FarmDir,condor, process, subprocess,tag_dir,check):
  cmsswBase = os.environ['CMSSW_BASE']
  shell_name = region + "_" + era + "_" + str(cfsf) + str(shift) + fin + str(label) + ".sh"
  with open('%s/%s'%(FarmDir,shell_name),'w') as shell:
    shell.write('#!/bin/bash\n')
    shell.write('WORKDIR=%s/src/Ele_ChargeFlipScaleFactor\n'%cmsswBase)
    shell.write('cd %s\n'%cmsswBase)
    shell.write('eval `scram r -sh`\n')
    shell.write('cd ${WORKDIR}\n')
    shell.write('python %s/src/Ele_ChargeFlipScaleFactor/flatten.py '%cmsswBase)
    shell.write('--region %s '%region)
    shell.write('--era %s '%era)
    shell.write('--isMC %d '%isMC)
    shell.write('--trig %d '%trig_sf)
    shell.write('--fin %s '%fin)
    shell.write('--xsec %f '%xs)
    shell.write('--cfsf %d '%cfsf)
    shell.write('--shift %d '%shift)
    shell.write('--from %d '%start)
    shell.write('--to %d '%end)
    shell.write('--label %d '%label)
    shell.write('--process %s '%process)
    shell.write('--subprocess %s '%subprocess)
  condor.write('cfgFile=%s\n'%shell_name)
  condor.write('queue 1\n')

  cfsf_name = "ApplyChargeFlipsf" if cfsf else "NotApplyChargeFlipsf"
  shift_name = "_Nominal/"
  if shift==1:
    shift_name = "_UP/"
  elif shift == -1:
    shift_name = "_DOWN/"
  fout_name = tag_dir + "era" + era + "/" + cfsf_name + shift_name + "h" + str(label) + "_" + fin
  if check == 1:
    try:
      fout = ROOT.TFile.Open(fout_name, "READ")
    except:
      print("Can not produce %s"%fout_name)
 

if __name__=='__main__':

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('-m', '--method', dest='method', help='[I] flatten w/o SF [II] flatten w SF [III] plot validation region', default='I', type='string')
  parser.add_option('-c', '--check', dest ='check', help='check if everything goes well', default='0', type='int')
  (args,opt) = parser.parse_args()

  FarmDir = os.environ['CMSSW_BASE'] + "/Farm/"
  os.system('mkdir -p %s'%FarmDir)
#  os.system('rm %s/*'%FarmDir)
  condor = open('%s/condor.sub'%FarmDir,'w')
  condor.write('output = %s/job_common.out\n'%FarmDir)
  condor.write('error  = %s/job_common.err\n'%FarmDir)
  condor.write('log    = %s/job_common.log\n'%FarmDir)
  condor.write('executable = %s/$(cfgFile)\n'%FarmDir)
  condor.write('requirements = (OpSysAndVer =?= "CentOS7")\n')
  condor.write('+JobFlavour = "tomorrow"\n')
  condor.write('+MaxRuntime = 7200\n')

#  Eras = ['2017','2018']
  Eras = ['2016postapv','2016apv','2017','2018']
  trigger_sf = 1

  chargeflip_sf = []
  directory = []
  tag_dir = ''
  region = ''

  if args.method == 'I':
    tag_dir = 'flatten/'
    region  = 'nominal'
    directory = ["NotApplyChargeFlipsf_Nominal"]
    chargeflip_sf = [(0,0)]

  elif (args.method == 'II' or args.method == 'V'):
    tag_dir = 'flatten/'
    region  = 'nominal'
    chargeflip_sf = [(1,1),(1,0),(1,-1)]
    directory = ["ApplyChargeFlipsf_Nominal", "ApplyChargeFlipsf_UP", "ApplyChargeFlipsf_DOWN"]

  elif (args.method == 'III'):
    tag_dir = 'validation/'
    region  = 'validation'
    chargeflip_sf = [(1,0)]
    directory = ["ApplyChargeFlipsf_Nominal"]

  elif (args.method == 'IV' or args.method == 'V'):
    tag_dir = 'validation/'
    region  = 'validation'
    chargeflip_sf = [(0,0),(1,1),(1,-1)]
    directory = ["NotApplyChargeFlipsf_Nominal", "ApplyChargeFlipsf_UP", "ApplyChargeFlipsf_DOWN"]

  for era in Eras:
    for di in directory:
      os.system('mkdir -p ' + tag_dir + 'era%s/%s'%(era,di))
      if not args.check:
        os.system('rm ' + tag_dir + 'era%s/%s/*'%(era,di))
    if(era == "2017"):
      path = '/eos/cms/store/group/phys_top/ExtraYukawa/TTC_version9/'
    elif(era == '2018'):
      path = '/eos/cms/store/group/phys_top/ExtraYukawa/2018/'
    elif(era == '2016postapv'):
      path = '/eos/cms/store/group/phys_top/ExtraYukawa/2016postapvMerged/'
    elif(era == '2016apv'):
      path = '/eos/cms/store/group/phys_top/ExtraYukawa/2016apvMerged/'

    jsonfile = open(os.path.join('data/sample_' + era + 'UL.json'))
    samples = json.load(jsonfile, encoding='utf-8', object_pairs_hook=OrderedDict).items()
    jsonfile.close()


    for process, desc in samples:
      #if not desc[2] == "WJet":
      #  continue
      dirs = os.listdir(path)
      for f in dirs:
        if((desc[1] and process == f.replace('.root','')) or (not desc[1] and process in f)):
          range_list = prepare_range(path,f)
          if(range_list is None):
            continue
          for i in range(len(range_list)-1):
            start = range_list[i]
            end = range_list[i+1]
            for cfsf,shift in chargeflip_sf:
              prepare_shell(region, era,desc[1],trigger_sf,cfsf,shift,f,desc[0],start,end,i, FarmDir,condor, desc[2],process,tag_dir,args.check)
  condor.close()
  if not args.check:
    os.system('condor_submit %s/condor.sub'%FarmDir)

