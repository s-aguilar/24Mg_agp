import numpy as np
import pandas as pd


# Read in the data into dataframe
df1 = pd.read_csv('peakAreasP1.csv')
df2 = pd.read_csv('peakAreasP2.csv')
df = pd.read_excel('24MgYields.xlsx',sheet_name='Yields_24Mg_ap')

pRun = df['Run'].values
pEalpha = df['Ea (keV)'].values
runToEalpha = {pRun[i]:pEalpha[i] for i in range(len(pRun))}

Ealpha = []

for ind, val in df2['Run'].iteritems():
    val = np.int64(val[4:])
    Ealpha.append(runToEalpha[val])


# Append new columns (E alpha) to dataframe, preserving the index
df1 = df1.assign(Ea=pd.Series(Ealpha,index=df1.index).values)
df2 = df2.assign(Ea=pd.Series(Ealpha,index=df2.index).values)

df1.to_csv('P1Yields.csv', sep=',')
df2.to_csv('P2Yields.csv', sep=',')
