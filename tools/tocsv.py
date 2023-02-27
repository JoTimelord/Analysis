# importing panda library
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename", help="Enter file name")

args=parser.parse_args()
file=args.filename


dataframe=pd.read_csv(file+".txt", header=None, index_col=False)

# setting column names
dataframe.columns = ['Path', 'Key', 'Year', 'Luminosity', 'Xsection', 'WeightedEventSum', 'Raw event count', 'Per event weight']
  
# storing this dataframe in a csv file
dataframe.to_csv(file+".csv", index=None)
