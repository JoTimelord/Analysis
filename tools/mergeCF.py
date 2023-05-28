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

SIG_SAMPLE_MAP = {
    "WZH": {
        "20UL16APV*":["VBSWZH*"],
        "20UL16-*":["VBSWZH*"],
        "20UL17*":["VBSWZH*"],
        "20UL18*":["VBSWZH*"],
    },
    "OSWWH": {
        "20UL16APV*":["VBSOSWWH*"],
        "20UL16-*":["VBSOSWWH*"],
        "20UL17*":["VBSOSWWH*"],
        "20UL18*":["VBSOSWWH*"],
    },
    "ZZH": {
        "20UL16APV*":["VBSZZH*"],
        "20UL16-*":["VBSZZH*"],
        "20UL17*":["VBSZZH*"],
        "20UL18*":["VBSZZH*"],
    },
    "WWH": {
        "20UL16APV*":["VBSWWH*"],
        "20UL16-*":["VBSWWH*"],
        "20UL17*":["VBSWWH*"],
        "20UL18*":["VBSWWH*"],
    },
}

def hadd_job(split_cmd):
    hadd = Popen(split_cmd, stdout=PIPE, stderr=PIPE)
    hadd.wait()

def merge(output_dir, merge_dir, sample_map, n_hadders=8):
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
        merge_file = f"{merge_dir}/{group_name}.root"
        hadd_cmds.append(["hadd", merge_file] + root_files_to_merge)
        if os.path.exists(merge_file):
            os.remove(merge_file)
        # Write list of files that were hadded
        with open(merge_file.replace(".root", ".txt"), "w") as f_out:
            for year, group in groups.items():
                f_out.write("{0}:\n{1}\n\n".format(year, '\n'.join(sorted(group))))

    # Run hadd jobs
    with Pool(processes=n_hadders) as pool:
        list(tqdm(pool.imap(hadd_job, hadd_cmds), total=len(hadd_cmds), desc="Executing hadds"))

    if merged_cutflows:
        return CutflowCollection(merged_cutflows)
    else:
        return None

def mergeAll(merge_dir, training_dir, sample_map, is_signal=False):
    # stage merge jobs
    hadd_cmds = []
    if (is_signal==True): merge_file=f"{training_dir}/signal.root"
    else: merge_file=f"{training_dir}/bkg.root"
    root_files_to_merge = []
    for group_name, group_map in sample_map.items():
        root_files = glob.glob(f"{merge_dir}/{group_name}.root")
        root_files_to_merge.extend(root_files)
    hadd_cmds.append(["hadd", merge_file] + root_files_to_merge)
    if os.path.exists(merge_file):
        os.remove(merge_file)
    print("Merging files: ", root_files_to_merge)
    with Pool(processes=1) as pool:
        list(tqdm(pool.imap(hadd_job, hadd_cmds), total=len(hadd_cmds), desc="Executing hadds"))

if __name__ == "__main__":
    # create hadded output directory
    bkg_output_dir="/home/users/joytzphysics/Analysis/outputs/raw_outputs"
    sig_output_dir="/home/users/joytzphysics/Analysis/outputs/raw_outputs"
    merge_dir="/home/users/joytzphysics/Analysis/outputs/hadded"
    os.makedirs(merge_dir, exist_ok=True)

    # Get Cutflow objects for background samples
    cutflows = merge(bkg_output_dir, merge_dir, BKG_SAMPLE_MAP)
    cutflows["TotalBkg"] = cutflows.sum()

    # Get Cutflow objects for signal samples
    cutflows_sig= merge(sig_output_dir, merge_dir, SIG_SAMPLE_MAP)
    cutflows += cutflows_sig
    cutflows["TotalSig"] = cutflows_sig.sum()
    cutflows.reorder(["WWH","WZH","OSWWH", "ZZH", "TotalSig", "TotalBkg", "DYJets", "TTX", "Others"])

    # hadd files for training input (signal vs. bkg)
    train_input_dir="/home/users/joytzphysics/Analysis/outputs/trainingInputs"
    os.makedirs(train_input_dir, exist_ok=True)
    mergeAll(merge_dir,train_input_dir,SIG_SAMPLE_MAP,True)

    # Write .cflow files
    for group_name, cutflow in cutflows.items():
        cutflow.write_cflow(f"{merge_dir}/{group_name}_cutflow.cflow")
    
    # Write to CSV files
    for terminal_cut_name in cutflows.terminal_cut_names:
        cutflows.write_csv(f"{merge_dir}/cutflow_{terminal_cut_name}.csv", terminal_cut_name)
        cutflows.write_txt(f"{merge_dir}/cutflow_{terminal_cut_name}.txt", terminal_cut_name)
        cutflows.write_tex(f"{merge_dir}/cutflow_{terminal_cut_name}.tex", terminal_cut_name)
