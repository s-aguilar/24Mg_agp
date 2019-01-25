import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

print('Attempting to create required directories: ')
try:
    os.mkdir('yieldPlots')
    os.mkdir('yieldPlots/P1')
    os.mkdir('yieldPlots/P2')
    os.mkdir('yieldPlots/A1')
except OSError:
    print ("Directories already exist!")
else:
    print('DONE!')



# Read in the data into dataframe
df1 = pd.read_csv('Yields/P1/p1Yields.csv')
df1 = df1.sort_values(by=['Ea'])
# df2 = pd.read_csv('Yields/P2/p2Yields.csv')
# df2 = df2.sort_values(by=['Ea'])
# """
# Extract the columns of the DataFrame as numpy arrays
p1Run = df1['Run'].values
p1Det = df1['Detector'].values
p1Yield = df1['Yield'].values
p1Yield_err = df1['Yield err'].values
p1Fit = df1['Fit Status'].values
p1Ealpha = df1['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
mask1Fit = (df1['Fit Status'] == 0)
# """

detectors = ['det_h0-0','det_h0-1','det_h0-2','det_h0-3','det_h0-4','det_h0-5',
            'det_h0-6','det_h0-7','det_h1-0','det_h1-1','det_h1-2','det_h1-3']

# """
for det in detectors:

    maskDet = ( (df1['Detector']==det) & mask1Fit )

    plt.clf()
    plt.plot(p1Ealpha[maskDet],p1Yield[maskDet],c='b',marker='.')
    plt.scatter(p1Ealpha[maskDet],p1Yield[maskDet],c='b',marker='.')
    # plt.errorbar(p1Ealpha[maskDet],p1Yield[maskDet],yerr=p1Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.ylim(0,.5)
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('p1 %s'%det)
    plt.savefig('yieldPlots/P1/p1_%s.png'%det,dpi=900)
    # plt.show()
# """

"""
# Extract the columns of the DataFrame as numpy arrays
p2Run = df2['Run'].values
p2Det = df2['Detector'].values
p2Yield = df2['Yield'].values
p2Yield_err = df2['Yield err'].values
p2Fit = df2['Fit Status'].values
p2Ealpha = df2['Ea'].values/1000    # Convert keV to MeV

# Fit Status == 0 -> Good Fit
# Fit Status == 1 -> Bad Fit
#
# Mask for which the fit was bad
mask2Fit = (df2['Fit Status'] == 0)

detectors = ['det_h0-0','det_h0-1','det_h0-2','det_h0-3','det_h0-4','det_h0-5',
            'det_h0-6','det_h0-7','det_h1-0','det_h1-1','det_h1-2','det_h1-3']

for det in detectors:

    maskDet = ( (df2['Detector']==det) & mask2Fit )

    plt.clf()
    plt.plot(p2Ealpha[maskDet],p2Yield[maskDet],c='b',marker='.')
    plt.errorbar(p2Ealpha[maskDet],p2Yield[maskDet],yerr=p2Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('p2 %s'%det)
    plt.savefig('yieldPlots/P2/p2_%s.png'%det,dpi=900)
    # plt.show()
# """
