'''
This script reads in the yield files produced from `thicknessAnalysis.C` as well as
the experimental online logbook file. A pandas DataFrame is used, and modified
to have the yields, the efficiency corrected yields, run number, and beam
energy. The DataFrame is then `cleaned` for bad data and when finalized it is
then saved to an excel file, for subsequent analysis and plotting in
crossSectionV2.py
'''
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pyFitFcns import gaussFitPlusBack,lorenFitPlusBack


# VERBOSITY
#   0 = off
#   1 = some
#   2 = all
verbose = 0


# Read online analysis logbook (only relevant columns)
df = pd.read_excel('27Al_p_a.xlsx',sheet_name='Sheet4',usecols=['Run','Ep (keV)'])

# Keep only relevant data
df.query('Run < 97 and Run > 76',inplace=True)


pRun = df['Run'].values
pEproton = df['Ep (keV)'].values
pEproton = pEproton.round(1) # round to 2 decimal points

# Thin target yield
dtypeDict = {'Yield':np.float,'Yield err':np.float,'Area':np.float,'Area err':np.float,
            'X2NDF':np.float,'isValid':np.int,'Status':np.int,'Q int':np.int,'Time':np.int}
ttyYield = pd.read_csv('Yields/Thick/Thick.csv',dtype=dtypeDict)

# Check to see if any fits are bad
print('\n\nCheck for bad fits (Empty is Good!):')
print(ttyYield.query('IsValid !=1 or Status != 0'))
ttyYield.query('IsValid == 1 or Status == 0',inplace=True)
ttyYield = ttyYield.drop(columns = ['IsValid','Status'])
print('\n\n')

print('Check for low stats (Empty is Good!):')
print(ttyYield.query('Area < 100'))
ttyYield.query('Area > 100',inplace=True)
print('\n\n')

print('Check for large uncertainties (Empty is Good!):')
print(ttyYield.query('`Area err` > Area'))
ttyYield.query('`Area err` < Area',inplace=True)
print('\n\n')

# Create dictionary, associating Run number to Proton energy
runToEproton = {pRun[i]:pEproton[i] for i in range(len(pRun))}

Eproton = []
RunNum = []
for ind, val in ttyYield['Run'].iteritems():
    val = np.int64(val[4:])
    Eproton.append(runToEproton[val])
    RunNum.append(val)

DetNum = []
for ind, val in ttyYield['Detector'].iteritems():
    val = np.int64(val[3:])
    DetNum.append(val)

# Columns are not useful as type string, get rid of them and reindex them
ttyYield = ttyYield.drop(columns = ['Run','Detector'])
ttyYield = ttyYield.assign(Ep=pd.Series(Eproton,index=ttyYield.index).values)
ttyYield = ttyYield.assign(Run=pd.Series(RunNum,index=ttyYield.index).values)
ttyYield = ttyYield.assign(Detector=pd.Series(DetNum,index=ttyYield.index).values)

# Convert from pulses to Coulombs to make it a true yield
q_e = 1.602e-19
scale = 1e-8    # 10^-8 C/pulse
q_corr = scale/(q_e)
yield_cor = ttyYield['Yield'].values/ q_corr
yield_err_cor = ttyYield['Yield err'].values/ q_corr


# Update the DataFrame with efficiency corrected yields
ttyYield = ttyYield.assign(Yield_cor=pd.Series(yield_cor,index=ttyYield.index).values)
ttyYield = ttyYield.assign(Yield_err_cor=pd.Series(yield_err_cor,index=ttyYield.index).values)
ttyYield = ttyYield.rename(columns={'Yield_cor': 'Yield cor', 'Yield_err_cor': 'Yield err cor'})
# print(ttyYield.head())

# ttyYield.query('Run!=77',inplace=True)  # Don't use run 77 (low stats)


# Detector names to loop through
for det in range(13):
    # print(det)
    det_df = ttyYield.query('Detector == %s'%det)
    # print(det_df.head(12))

    ePoints = np.linspace(min(det_df['Ep']),max(det_df['Ep']),1000)

    # pG = np.array([992.7,5e-13,5e-1,1e-18,1e-16]) # initial guess of parameters
    # poptG, pcovG = curve_fit(f=gaussFitPlusBack, xdata=det_df['Ep'],
    #     ydata=det_df['Yield cor'],sigma=det_df['Yield err cor'],
    #     bounds=([991,1e-15,1e-1,-np.inf,-np.inf],[993,1e-11,1.3,np.inf,np.inf]),
    #     method='trf',absolute_sigma=True)
    # print(poptG)
    # yfitG = gaussFitPlusBack(ePoints, *poptG)    # Fit line

    pL = np.array([992.7,5,7e-1,0,0]) # initial guess of parameters
    poptL, pcovL = curve_fit(f=lorenFitPlusBack, xdata=det_df['Ep'],
        ydata=det_df['Area']/det_df['Time'],p0=pL,sigma=det_df['Area err']/det_df['Time'],
        absolute_sigma=True,bounds=([991,1,4e-1,-np.inf,-np.inf],[993,15,2,np.inf,np.inf]),
        method='trf')
    yfitL = lorenFitPlusBack(ePoints, *poptL)    # Fit line
    # print(poptL)
    # print(max(yfitL))
    plt.errorbar(det_df['Ep'],det_df['Area']/det_df['Time'],yerr=det_df['Area err']/det_df['Time'],
                    fmt='.',c='k',label='Data')
    # plt.plot(ePoints,yfitG,c='r',label='G+B fit')
    plt.plot(ePoints,yfitL,c='r',label='L+B fit')
    plt.xlabel('E$_{p} (keV)$')
    plt.ylabel('Normlalized counts')
    plt.title('Width: %4.3f keV'%poptL[2])
    plt.legend()
    # plt.show()
    plt.savefig('plots/thinTargetYield_det%s.png'%det)
    plt.clf()
    # exit()

    print('%4.2f'%(poptL[2]/(47.26*.0037011)),'ug/cm2')





exit()
