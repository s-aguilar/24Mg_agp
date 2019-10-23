import numpy as np
import pandas as pd
from scipy import special
import matplotlib.pyplot as plt
from AtomicMassTable import GetElement

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# Use to switch directory paths to use (for the main or TESTING)
import os
currentDir = os.path.realpath('')
parentDir = os.path.realpath('..')



# CHECK WITH JAMES ABOUT THAT FACTOR OF 2 BEING SEENNN
# CHECK THE TO MAKE SURE


"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir


"""
Calculate the reaction rates given the AZURE2 extrapolated angle integrated
cross section as a function of energy
"""


# GetElement returns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]
_barnsTocm2 = 1E-24         # cm2
_c = 2.99792458E10          # cm / s
_u = 931.494/_c**2          # MeV / c^2
_NA = 6.02214E23            # particles / mol


rxnRateCONST = _barnsTocm2 * _NA * (8/(np.pi*_u))**.5 * (8.617E-2)**-1.5


# Read in the angle integrated cross-sections calculated via azure calculate w/o data
colNames = ['E_Cm','E_Ex','Angle_Cm',
                'Fit_Cm_Cross','Fit_Cm_S_factor']
p0Path = os.path.join(desiredDir,'rMatrix/AZUREOut_aa=2_R=1_angInt.extrap')
dfp0 = pd.read_table(p0Path,sep='\s+',names=colNames)

p1Path = os.path.join(desiredDir,'rMatrix/AZUREOut_aa=2_R=3_angInt.extrap')
dfp1 = pd.read_table(p1Path,sep='\s+',names=colNames)

p2Path = os.path.join(desiredDir,'rMatrix/AZUREOut_aa=2_R=4_angInt.extrap')
dfp2 = pd.read_table(p2Path,sep='\s+',names=colNames)

fname = 'rate_Temps.dat'
temperature = np.loadtxt( fname )

# Reduced mass
mu = _m24Mg*_m4He/(_m24Mg+_m4He)

dict_channels = {'p0':dfp0,'p1':dfp1,'p2':dfp2}
channels = ['p0','p1','p2']
for ch in channels:

    plt.clf()

    cross = dict_channels[ch]['Fit_Cm_Cross'].values    # angle integrated cross
    E = dict_channels[ch]['E_Cm'].values        # E_lab
    # E = E*(_m4He+_m24Mg)/_m24Mg                 # E_cm

    rxnRate = []
    dE = .001       # 1 keV steps
    for T in temperature:
        integrand = np.trapz(cross*np.exp(-11.604*E/T)*E,E)
        rate =  rxnRateCONST * mu**(-.5) * T**(-1.5) * integrand
        rxnRate.append(rate)

    rxnRate = np.array(rxnRate)

    col = ["T9","Rate"]

    myp0rates = pd.read_excel('legendre_out/DATA/p0/a0/p0_rates.xlsx',sheet_name='Sheet1',header=0)
    myp1rates = pd.read_excel('legendre_out/DATA/p1/a0/p1_rates.xlsx',sheet_name='Sheet1',header=0)
    myp2rates = pd.read_excel('legendre_out/DATA/p2/a0/p2_rates.xlsx',sheet_name='Sheet1',header=0)

    azp0rates = pd.read_table('p0rates.out',sep='\s+')
    azp1rates = pd.read_table('p1rates.out',sep='\s+')
    azp2rates = pd.read_table('p2rates.out',sep='\s+')

    # print(myp1rates.columns)
    # print(azp1rates.columns)

    temp0 = myp0rates['T9'].values
    temp1 = myp1rates['T9'].values
    temp2 = myp2rates['T9'].values
    aztemp0 = azp0rates['T9'].values
    aztemp1 = azp1rates['T9'].values
    aztemp2 = azp2rates['T9'].values

    rate0 = myp0rates['Rate'].values
    rate1 = myp1rates['Rate'].values
    rate2 = myp2rates['Rate'].values
    azrate0 = azp0rates['Rate'].values
    azrate1 = azp1rates['Rate'].values
    azrate2 = azp2rates['Rate'].values

    # plt.scatter(temp0,rate0,c='silver')
    plt.plot(temp0,rate0,color='silver',label='az extrap p0')

    # plt.scatter(temp1,rate1,c='k')
    plt.plot(temp1,rate1,color='k',label='az extrap p1')

    # plt.scatter(temp2,rate2,c='grey')
    plt.plot(temp2,rate2,color='grey',label='az extrap p2')

    # plt.scatter(aztemp0,azrate0,c='g')
    plt.plot(aztemp0,azrate0,color='g',label='azure rate p0')

    # plt.scatter(aztemp1,azrate1,c='fuchsia')
    plt.plot(aztemp1,azrate1,color='fuchsia',label='azure rate p1')

    # plt.scatter(aztemp2,azrate2,c='b')
    plt.plot(aztemp2,azrate2,color='b',label='azure rate p2')

    plt.grid(b=True, which='both', axis='both')


    # plt.scatter(temperature,rxnRate,c='coral')
    plt.plot(temperature,rxnRate,label='trapz',color='coral')
    plt.legend()
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-30,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.savefig('%s_ReactionRate_azureAngleIntegrated.png'%ch,dpi=900)



    # Plot ratio of azures angle integrated extrapolated cross-section calculated reaction rate (a la Sebastian)
    # vs azures default calculated reaction rate
    if ch == 'p0':
        plt.clf()
        y = rxnRate/azrate0
        plt.plot(aztemp0,y,color='k')
        plt.xlim(0,10)
        plt.ylim(0,5)
        plt.ylabel('Ratio of Reaction Rate',fontsize=14)
        plt.xlabel('Temperature (T9)',fontsize=14)
        plt.title('p$_{0}$ Reaction Rates Ratio',fontsize=20)
        eq = r'$\frac{Rate~p_{{0}_{angInt}}}{Rate~p_{0}}$'
        plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
        plt.grid(b=True, which='both', axis='both')

        savePath = os.path.join(desiredDir,'p0_RatioReactionRates_azureAngleIntegrated.png')
        plt.savefig(savePath,dpi=300)



    df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)
    df.to_excel('legendre_out/DATA/%s/a0/%s_rates_extrap_angInt.xlsx'%(ch,ch))
