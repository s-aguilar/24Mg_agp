'''
This script reads in the yield files produced from `xxyields.C` as well as
the experimental online logbook file. A pandas DataFrame is used, and modified
to have the yields, the efficiency corrected yields, run number, and beam
energy. The DataFrame is then `cleaned` for bad data and when finalized it is
then saved to an excel file, for subsequent analysis and plotting in
crossSectionExp.py
'''
import time
import numpy as np
import pandas as pd


# VERBOSITY
#   0 = off
#   1 = some
#   2 = all
verbose = 0

# Read in the data into dataframe
df1 = pd.read_csv('Yields/P1/_P1.csv')
df2 = pd.read_csv('Yields/P2/_P2.csv')
df3 = pd.read_csv('Yields/A1/_A1.csv')

# Read online analysis logbook
df = pd.read_excel('24MgYields.xlsx',sheet_name='Yields_24Mg_ap')


# Read in efficiencies
dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')

pRun = df['Run'].values
pEalpha = df['Ea (keV)'].values
pEalpha = pEalpha.round(1) # round to 1 decimal point (100s of eV precision)


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
