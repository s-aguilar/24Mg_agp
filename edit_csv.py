'''
This script reads in the yield files produced from `xxyields_pa.C` as well as
the experimental online logbook file. A pandas DataFrame is used, and modified
to have the yields, the efficiency corrected yields, run number, and beam
energy. The DataFrame is then `cleaned` for bad data and when finalized it is
then saved to an excel file, for subsequent analysis and plotting in
crossSectionV2.py
'''
import time
import numpy as np
import pandas as pd


# Read in the data into dataframe
df1 = pd.read_csv('Yields/P1/_P1.csv')
df2 = pd.read_csv('Yields/P2/_P2.csv')
df3 = pd.read_csv('Yields/A1/_A1.csv')
df = pd.read_excel('24MgYields.xlsx',sheet_name='Yields_24Mg_ap')

pRun = df['Run'].values
pEalpha = df['Ea (keV)'].values
pEalpha = pEalpha.round(3) # round to 2 decimal point
runToEalpha = {pRun[i]:pEalpha[i] for i in range(len(pRun))}

Ealpha = []

for ind, val in df1['Run'].iteritems():
    val = np.int64(val[4:])
    Ealpha.append(runToEalpha[val])


# Append new columns (E alpha) to dataframe, preserving the index
df1 = df1.assign(Ea=pd.Series(Ealpha,index=df1.index).values)
df2 = df2.assign(Ea=pd.Series(Ealpha,index=df2.index).values)
df3 = df3.assign(Ea=pd.Series(Ealpha,index=df3.index).values)

df1.to_csv('Yields/P1/p1Yields.csv', sep=',')
df2.to_csv('Yields/P2/p2Yields.csv', sep=',')
df3.to_csv('Yields/A1/a1Yields.csv', sep=',')
print('DONE!')
