import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Read in the data into dataframe
df1 = pd.read_csv('P1Yields.csv')
df2 = pd.read_csv('P2Yields.csv')
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
    plt.scatter(p1Ealpha[maskDet],p1Yield[maskDet],c='b',marker='.')
    # plt.errorbar(p1Ealpha[maskDet],p1Yield[maskDet],yerr=p1Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('p1 %s'%det)
    plt.savefig('yieldPlots/p1/p1_%s.png'%det,dpi=1200)
    plt.clf()
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
    plt.scatter(p2Ealpha[maskDet],p2Yield[maskDet],c='b',marker='.')
    plt.errorbar(p2Ealpha[maskDet],p2Yield[maskDet],yerr=p2Yield_err[maskDet],fmt='b.')
    # plt.yscale('log')
    plt.xlabel('$E_{\\alpha}$ (MeV)')
    plt.ylabel('Yield')
    plt.title('p2 %s'%det)
    plt.savefig('yieldPlots/p2/p2_%s.png'%det,dpi=1200)
    plt.clf()
    # plt.show()
# """
