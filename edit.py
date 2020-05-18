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


# Read online analysis logbook
df_log = pd.read_excel('24MgYields.xlsx',sheet_name='Yields_24Mg_ap')

channels = ['P1','P2','P3','A1','O17_1','O17_ng']
for chan in channels:

    # Read in the data into dataframe
    df = pd.read_csv('Yields/%s/_%s.csv'%(chan.upper(),chan))


    # Read in efficiencies
    # dfeff = pd.read_csv('calibration/csv/detectorEfficiencies.csv')

    pRun = df_log['Run'].values
    pEalpha = df_log['Ea (keV)'].values
    pEalpha = pEalpha.round(1) # round to 1 decimal point (100s of eV precision)


    runToEalpha = {pRun[i]:pEalpha[i] for i in range(len(pRun))}

    Ealpha = []
    runlist = []

    for ind, val in df['Run'].iteritems():
        val = np.int64(val[4:])
        Ealpha.append(runToEalpha[val])
        runlist.append(val)


    # Append new columns (E alpha) to dataframe, preserving the index
    df = df.assign(Ea=pd.Series(Ealpha,index=df.index).values)
    df = df.assign(runnum=pd.Series(runlist,index=df.index).values)###

    # Shift energy by +11.5 keV
    df = df.apply(lambda x: x + 11.5 if x.name == 'Ea' else x)

    df.to_csv('Yields/%s/%sYields.csv'%(chan,chan.lower()), sep=',')
    print('Finished with: %s'%chan)

print('DONE!')
