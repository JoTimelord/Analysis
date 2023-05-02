import os
import json
import glob
import uproot
from scale1fb import sumGenWeights, nevents

# output jobs file name
output_name='/home/users/joytzphysics/Analysis/jobs/.jobs_fullyMerge_v2.txt'

# input MC sample directories
bkg_mc_dir='/ceph/cms/store/user/jguiang/VBSVHSkim/bkg_2oslep_2ak4_1ak8_v2'
sig_mc_dir='/ceph/cms/store/user/jguiang/VBSVHSkim/sig_2oslep_2ak4_1ak8_v2'

# bkg/sig year prefixes & luminosity
bkg_years={"20UL16NanoAODv9": 16.81,
           "20UL16NanoAODAPVv9": 19.52, 
           "20UL17NanoAODv9": 41.48,
           "20UL18NanoAODv9": 59.83}

sig_years={"20UL16APV-106X":19.52,
           "20UL16-106X": 16.81,
           "20UL17-106X": 41.48,
           "20UL18-106X": 59.83}

# save all desired background file name keys
with open('/home/users/joytzphysics/Analysis/data/xsections.json',"r") as file:
    bkg_keys=json.load(file)

# signal file name keys: xsections
sig_keys={"VBSWZH_Inclusive_4f":1.67403,
          "VBSWWH_Inclusive_4f":1.50253,
          "VBSZZH_Inclusive_4f":1.06028,
          "VBSOSWWH_Inclusive_4f":2.52444}

with open(output_name,'w') as f:
    for key, xsec in bkg_keys.items():
        print(f"dealing with {key}")
        for yr, lumi in bkg_years.items():
            summed_wgt=0
            raw_events=0
            bkg_file_lists=glob.glob(f"{bkg_mc_dir}/{key}*{yr}*")
            for file in bkg_file_lists:
                print(f"dealing with {file}")
                raw_events+=nevents(file)
                summed_wgt+=sumGenWeights(file)
            if len(bkg_file_lists)==0: 
                print(f"Process with year {yr} and key {key} cannot be found.")
                continue
            scale1fb=xsec*lumi*1000/summed_wgt
            i=0
            for filepath in glob.iglob(f"{bkg_mc_dir}/{key}*{yr}*"+"**/*.root", recursive=True):
                f.write(f"./fullyMerge -t Events -d outputs/output_fullyMerge_2oslep_2ak4_1ak8_v2 -s {scale1fb} -n {key}_{yr}_{i} -T tree {filepath}\n")
                i+=1
    for key,xsec in sig_keys.items():
        print(f"dealing with {key}")
        for yr, lumi in sig_years.items():
            summed_wgt=0
            raw_events=0
            sig_file_lists=glob.glob(f"{sig_mc_dir}/{key}*{yr}*")
            for file in sig_file_lists:
                print(f"dealing with {file}")
                raw_events+=nevents(file)
                summed_wgt+=sumGenWeights(file)
            if len(sig_file_lists)==0: 
                print(f"Process with year {yr} and key {key} cannot be found.")
                continue
            scale1fb=xsec*lumi*1000/summed_wgt
            i=0
            for filepath in glob.iglob(f"{sig_mc_dir}/{key}*{yr}*"+"**/*.root", recursive=True):
                f.write(f"./fullyMerge -t Events -d outputs/output_fullyMerge_2oslep_2ak4_1ak8_v2 -s {scale1fb} -n {key}_{yr}_{i} -T tree {filepath}\n")
                i+=1
