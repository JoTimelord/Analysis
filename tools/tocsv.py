# importing panda library
import pandas as pd
  
# readinag given csv file
# and creating dataframe
dataframe1 = pd.read_csv("weightsbkg.txt", header=None, index_col=False)

# setting column names
dataframe1.columns = ['Path', 'Key', 'Year', 'Luminosity', 'Xsection', 'WeightedEventSum', 'Raw event count', 'Per event weight']
  
# storing this dataframe in a csv file
dataframe1.to_csv('weightsbkg.csv',
                  index=None)


dataframe2 = pd.read_csv("weightssig.txt", header=None, index_col=False)

# setting column names
dataframe2.columns = ['Path', 'Key', 'Year', 'Luminosity', 'Xsection', 'WeightedEventSum', 'Raw event count', 'Per event weight']
  
# storing this dataframe in a csv file
dataframe2.to_csv('weightssig.csv',
                  index = None)