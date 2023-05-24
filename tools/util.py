import uproot
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

# check if there are duplicate events in input directory
# pass in directory which stores the root files and tree name
# return (dataframe, if_repeated truth series) and prints out if there's duplicate events
def checkDuplicateinRoot(dir_name, tree_name):
    list_of_df=[]
    list_of_files=glob.glob(dir_name+"/*.root")
    for file in list_of_files:
        path=file+f":{tree_name}"
        tree=uproot.open(path)
        df_append=tree.arrays(library='pd')
        list_of_df.append(df_append)
    df=pd.concat(list_of_df)
    repeated=df.duplicated()
    if repeated.any()==True: print("There are duplicate events.")
    else: print("There are no duplicate events.")
    return (df,repeated)


def unique(data):
    unique_cols, indexes, counts = np.unique(data.numpy(), axis=1, return_index=True, return_counts=True)
    is_duplicate = np.size(indexes) != np.size(data.numpy(), axis=1)
    num_duplicate_cols = np.max(counts)

    # Print whether or not there are duplicate columns
    if is_duplicate:
        print("There are duplicate columns.")
        print("There are at most", num_duplicate_cols, "duplicate columns.")
    else:
        print("There are no duplicate columns.")


