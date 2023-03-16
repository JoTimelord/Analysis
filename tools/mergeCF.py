#!/usr/bin/env python3

import os
import glob
import argparse
from tqdm import tqdm
from multiprocessing import Pool
from subprocess import Popen, PIPE 
from cutflow import Cutflow, CutflowCollection

BKG_SAMPLE_MAP = {
    "DYJets": {
        "20UL16NanoAODv9*":["DY*"],
        "20UL16NanoAODAPVv9*":["DY*"],
        "20UL17NanoAODv9*":["DY*"],
        "20UL18NanoAODv9*":["DY*"],
    },
    "Others": {
        "20UL16NanoAODv9*":["WW*", "WZ*", "ZZ*", "EWKW*WToQQ*", "EWKZ*ZToNuNu*", "EWKZ*ZToLL*", "EWKZ*ZToQQ*", "EWKW*WToLNu*", 
                            "TTToSemiLep*", "ST*", "VHToNonbb*", "WminusH*", "WplusH*", "ZH_HToBB*", "ggZH_HToBB*", "WJets*"],
        "20UL16NanoAODAPVv9*":["WW*", "WZ*", "ZZ*", "EWKW*WToQQ*", "EWKZ*ZToNuNu*", "EWKZ*ZToLL*", "EWKZ*ZToQQ*", "EWKW*WToLNu*",
                               "TTToSemiLep*", "ST*", "VHToNonbb*", "WminusH*", "WplusH*", "ZH_HToBB*", "ggZH_HToBB*", "WJets*"],
        "20UL17NanoAODv9*":["WW*", "WZ*", "ZZ*", "EWKW*WToQQ*", "EWKZ*ZToNuNu*", "EWKZ*ZToLL*", "EWKZ*ZToQQ*", "EWKW*WToLNu*", 
                            "TTToSemiLep*", "ST*", "VHToNonbb*", "WminusH*", "WplusH*", "ZH_HToBB*", "ggZH_HToBB*", "WJets*"],
        "20UL18NanoAODv9*":["WW*", "WZ*", "ZZ*", "EWKW*WToQQ*", "EWKZ*ZToNuNu*", "EWKZ*ZToLL*", "EWKZ*ZToQQ*", "EWKW*WToLNu*", 
                            "TTToSemiLep*", "ST*", "VHToNonbb*", "WminusH*", "WplusH*", "ZH_HToBB*", "ggZH_HToBB*", "WJets*"],
    },
    "TTX": {
        "20UL16NanoAODv9*":["ttH*", "TTW*", "TTZ*", "TTbb*", "TTToHadronic*", "TTTo2L*"],
        "20UL16NanoAODAPVv9*":["ttH*", "TTW*", "TTZ*", "TTbb*", "TTToHadronic*", "TTTo2L*"],
        "20UL17NanoAODv9*":["ttH*", "TTW*", "TTZ*", "TTbb*", "TTToHadronic*", "TTTo2L*"],
        "20UL18NanoAODv9*":["ttH*", "TTW*", "TTZ*", "TTbb*", "TTToHadronic*", "TTTo2L*"],
    },
}

SIG_SAMPLE_MAP_2 = {
    "WZH": {
        "20UL16APV*":["WZH*"],
        "20UL16*":["WZH*"],
        "20UL17*":["WZH*"],
        "20UL18*":["WZH*"],
    },
    "OSWWH": {
        "20UL16APV*":["OSWWH*"],
        "20UL16*":["OSWWH*"],
        "20UL17*":["OSWWH*"],
        "20UL18*":["OSWWH*"],
    },
    "ZZH": {
        "20UL16APV*":["ZZH*"],
        "20UL16*":["ZZH*"],
        "20UL17*":["ZZH*"],
        "20UL18*":["ZZH*"],
    },
    "WWH": {
        "20UL16APV*":["WWH*"],
        "20UL16*":["WWH*"],
        "20UL17*":["WWH*"],
        "20UL18*":["WWH*"],
    },
}

SIG_SAMPLE_MAP_1 = {
    "OSWWH_C2V_3": {
        "20UL16APV*": ["VBSOSWWH_incl_C2V_3*"],
        "20UL16*": ["VBSOSWWH_incl_C2V_3*"],
        "20UL17*": ["VBSOSWWH_incl_C2V_3*"],
        "20UL18*": ["VBSOSWWH_incl_C2V_3*"],
    },
    "WZH_C2V_3": {
        "20UL16APV*": ["VBSWZH_incl_C2V_3*"],
        "20UL16*": ["VBSWZH_incl_C2V_3*"],
        "20UL17*": ["VBSWZH_incl_C2V_3*"],
        "20UL18*": ["VBSWZH_incl_C2V_3*"],
    },
    "ZZH_C2V_3": {
        "20UL16APV*": ["VBSZZH_incl_C2V_3*"],
        "20UL16*": ["VBSZZH_incl_C2V_3*"],
        "20UL17*": ["VBSZZH_incl_C2V_3*"],
        "20UL18*": ["VBSZZH_incl_C2V_3*"],
    }
}

def hadd_job(split_cmd):
    hadd = Popen(split_cmd, stdout=PIPE, stderr=PIPE)
    hadd.wait()

def merge(output_dir, sample_map, n_hadders=8):
    # Collect cutflows and stage merge jobs
    hadd_cmds = []
    merged_cutflows = {}
    for group_name, group_map in sample_map.items():
        groups = {}
        root_files_to_merge = []
        for year, sample_list in group_map.items():
            groups[year] = []
            for sample_name in sample_list:
                cflow_files = glob.glob(f"{output_dir}/{sample_name}{year}_Cutflow.cflow")
                for cflow_file in cflow_files:
                    if group_name in merged_cutflows.keys():
                        merged_cutflows[group_name] += Cutflow.from_file(cflow_file)
                    else:
                        merged_cutflows[group_name] = Cutflow.from_file(cflow_file)
                root_files = glob.glob(f"{output_dir}/{sample_name}{year}.root")
                for root_file in root_files:
                    root_files_to_merge.append(root_file)
                    groups[year].append(root_file.split("/")[-1].replace(".root", ""))
        # Stage merge (hadd) jobs
        output_file = f"{output_dir}/{group_name}.root"
        hadd_cmds.append(["hadd", output_file] + root_files_to_merge)
        if os.path.exists(output_file):
            os.remove(output_file)
        # Write list of files that were hadded
        with open(output_file.replace(".root", ".txt"), "w") as f_out:
            for year, group in groups.items():
                f_out.write("{0}:\n{1}\n\n".format(year, '\n'.join(sorted(group))))

    # Run hadd jobs
    with Pool(processes=n_hadders) as pool:
        list(tqdm(pool.imap(hadd_job, hadd_cmds), total=len(hadd_cmds), desc="Executing hadds"))

    if merged_cutflows:
        return CutflowCollection(merged_cutflows)
    else:
        return None



if __name__ == "__main__":
    # create hadded output directory
    # output_dir="/home/users/joytzphysics/Analysis/outputs/output_semiMerge_2oslep_2ak4_1ak8"
    # output_dir="/home/users/joytzphysics/Analysis/outputs/output_semiMerge_2oslep_2ak4_1ak8_v2"
    # output_dir="/home/users/joytzphysics/Analysis/outputs/output_fullyMerge_1lep_1ak8_2ak4_v1"
    output_dir="/home/users/joytzphysics/Analysis/outputs/output_fullyMerge_2oslep_2ak4_1ak8_v2"
    os.makedirs(output_dir, exist_ok=True)

    # Get Cutflow objects for background samples
    cutflows = merge(output_dir, BKG_SAMPLE_MAP)
    cutflows["TotalBkg"] = cutflows.sum()
    # Get Cutflow objects for signal samples
    # cutflows = merge(output_dir, SIG_SAMPLE_MAP_1)
    cutflows_sig= merge(output_dir, SIG_SAMPLE_MAP_2)
    cutflows += cutflows_sig
    cutflows["TotalSig"] = cutflows_sig.sum()
    cutflows.reorder(["WWH","WZH","OSWWH", "ZZH", "TotalSig", "TotalBkg", "DYJets", "TTX", "Others"])
    # cutflows.reorder(["WWDilep","WZHDilep","OSWWHDilep", "ZZHDilep","OSWWH_C2V_3", "WZH_C2V_3", "ZZH_C2V_3"])

    # Write .cflow files
    for group_name, cutflow in cutflows.items():
        cutflow.write_cflow(f"{output_dir}/{group_name}_cutflow.cflow")
    
    # Write to CSV files
    for terminal_cut_name in cutflows.terminal_cut_names:
        cutflows.write_csv(f"{output_dir}/cutflow_{terminal_cut_name}.csv", terminal_cut_name)
        cutflows.write_txt(f"{output_dir}/cutflow_{terminal_cut_name}.txt", terminal_cut_name)
        cutflows.write_tex(f"{output_dir}/cutflow_{terminal_cut_name}.tex", terminal_cut_name)
