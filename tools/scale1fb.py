#!/bin/env python3
import os, sys
import argparse
import uproot


def sumGenWeights(dirname):
    nevents = 0
    if os.path.isfile(dirname):
        with uproot.open(dirname) as f:
            t = f.get("Runs")
            nevents = t["genEventSumw"].array(library="np").sum()
    elif os.path.isdir(dirname):
        for root, dirs, files in os.walk(dirname):
            for file in files:
                name = os.path.join(root, file)
                with uproot.open(name) as f:
                    t = f.get("Runs")
                    n_events = t["genEventSumw"].array(library="np").sum()
                    #print ("Event " + name + " produced ", n_events, " events.")
                    nevents += n_events
    else:
        raise SystemError
    return nevents

def nevents(dirname, tree_name):
    nevents = 0
    if os.path.isfile(dirname):
        if dirname.endswith(".root"):
            f = uproot.open(dirname)
            t = f[tree_name]
            nevents = t.num_entries
    elif os.path.isdir(dirname):
        for root, dirs, files in os.walk(dirname):
            for file in files:
                if file.endswith(".root"):
                    name = os.path.join(root, file)
                    # print (name)
                    fi = uproot.open(name)
                    tree = fi[tree_name]
                    # print (tree.GetEntries())
                    nevents += tree.num_entries
    # print (nevents)
    return nevents


if __name__ == "__main__":
# Initialize parser
    parser = argparse.ArgumentParser()

# Adding optional argument
    parser.add_argument("-f", "--file", help = "Enter directory/path name")
    parser.add_argument("-c", "--choice", help = "Raw event (R), Weighted event (W)")
    parser.add_argument("-t", "--tree", help="Tree name")

# Read arguments from command line
    args = parser.parse_args()
    option = args.choice
    path = args.file
    tree_name = args.tree
    number = -999
    if option == 'R':
        number = nevents(path, tree_name)
    elif option == 'W':
        number = sumGenWeights(path)
    print (number)


