import numpy as np
import pandas as pd
from scipy import special
import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation
from AtomicMassTable import GetElement
import imageio

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
Calculate the reaction rates given the AZURE2 extrapolated angle integrated FINE
cross section as a function of energy, for the energy ranges matching data
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

fname = 'iliadis_Temps.dat'
temperature = np.loadtxt( fname )

# Reduced mass
mu = _m24Mg*_m4He/(_m24Mg+_m4He)

dict_channels = {'p0':dfp0,'p1':dfp1,'p2':dfp2}
channels = ['p0','p1','p2']
for ch in channels:

    plt.clf()

    cross = dict_channels[ch]['Fit_Cm_Cross'].values    # angle integrated cross
    E = dict_channels[ch]['E_Cm'].values        # E_Cm


    # e_mask1 = E <= 3.459
    # e_mask2 = E >= 3.66
    # e_mask1 = E >= 3.42454
    e_mask1 = (E <= 4.66271)
    e_mask2 = (E >= 4.66271)
    # not_e_mask1 = (E <= 3.42454)
    # not_e_mask2 = (E >= 4.66271)
    # Plot the angle-integrated AZURE2 extrapolated cross-section
    plt.clf()
    plt.plot(E[e_mask1],cross[e_mask1],c='b',label='UP TO')
    plt.plot(E[e_mask2],cross[e_mask2],c='b',label='ABOVE')
    # plt.scatter(E,cross,c='b')
    plt.xlabel('E_CM (MeV)',fontsize=14)
    plt.ylabel('Angle Integrated Cross-Section (barns)', fontsize=14)
    plt.yscale('log')
    plt.legend()
    plt.title('%s channel - Angle-Integrated Finer Excitation Curve' % ch,fontsize=20)
    plotPath = os.path.join(desiredDir,'%s_angIntExtrapFiner10ExCurveJAMES.png' % ch)
    plt.savefig(plotPath,dpi=300)
    # plt.show()
    plt.clf()


    rxnRateBELOW = []
    rxnRateABOVE = []
    dE = 0
    if ch == 'p0':
        dE = .0006       # 1 keV steps IN LAB energy
        dE = dE * (_m24Mg/(_m24Mg+_m4He))         # In CENTER-OF-MASS energy
    elif ch == 'p1':
        dE = .0006       # 1 keV steps IN LAB energy
        dE = dE * (_m24Mg/(_m24Mg+_m4He))         # In CENTER-OF-MASS energy
    elif ch == 'p2':
        dE = .0001       # 1 keV steps IN LAB energy
        dE = dE * (_m24Mg/(_m24Mg+_m4He))         # In CENTER-OF-MASS energy

    _integralList = []
    filenames = []

    for T in temperature:
        integrand = cross*np.exp(-11.604*E/T)*E
        # integrand_err = (cross*.1) * E * np.exp(-11.604*E/T) # MAYBE? give it 10%
        integralBELOW = np.trapz(integrand[e_mask1],E[e_mask1])
        integralABOVE = np.trapz(integrand[e_mask2],E[e_mask2])

        rateBELOW =  rxnRateCONST * mu**(-.5) * T**(-1.5) * integralBELOW
        rxnRateBELOW.append(rateBELOW)

        rateABOVE =  rxnRateCONST * mu**(-.5) * T**(-1.5) * integralABOVE
        rxnRateABOVE.append(rateABOVE)


    plt.plot(temperature,rxnRateBELOW,label='%s-angInt-BELOW'%ch,color='b')
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-20,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('AZURE2 BELOW %s - Reaction Rate'%ch,fontsize=20)
    plt.legend()
    plt.savefig('%s_ReactionRate_azureAngleIntegratedFiner10JAMESBELOW.png'%ch,dpi=900)

    plt.plot(temperature,rxnRateABOVE,label='%s-angInt-ABOVE'%ch,color='r')
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-20,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('AZURE2 OUT %s - Reaction Rate'%ch,fontsize=20)
    plt.legend()
    plt.savefig('%s_ReactionRate_azureAngleIntegratedFiner10JAMESABOVE.png'%ch,dpi=900)






    # Plot ratio of azures angle integrated extrapolated cross-section calculated reaction rate for inner and outer
    # vs Iliadis reaction rate

    pathToREACLIB = os.path.join(desiredDir,'24MgREACLIB.xlsx')
    df = pd.read_excel(pathToREACLIB,sheet_name='il10',header=0)
    _tempIliadis = df['T9'].values
    _rateIliadis = df['Rate'].values

    plt.clf()
    y = rxnRateBELOW/_rateIliadis
    plt.plot(_tempIliadis,y,color='k')
    plt.xlim(0,10)
    yy = [1,1]
    xx = np.linspace(0,10,2)
    plt.plot(xx,yy)
    plt.yscale('log')
    plt.ylim(1e-3,1)
    plt.ylabel('Ratio of Reaction Rate',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('p$_{%s}$ Reaction Rates Ratio - BELOW'%ch,fontsize=20)
    eq = r'$\frac{Rate~p_{{%s}_{angInt}}}{Rate~p_{0}}$'%ch
    plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
    # plt.grid(b=True, which='both', axis='both')

    savePath = os.path.join(desiredDir,'%s_RatioReactionRates_JAMESBELLOW.png'%ch)
    plt.savefig(savePath,dpi=300)
    plt.clf()


    y = rxnRateABOVE/_rateIliadis
    plt.plot(_tempIliadis,y,color='k')
    plt.xlim(0,10)
    yy = [1,1]
    xx = np.linspace(0,10,2)
    plt.plot(xx,yy)
    plt.yscale('log')
    plt.ylim(1e-3,1)
    plt.ylabel('Ratio of Reaction Rate',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('p$_{%s}$ Reaction Rates Ratio - ABOVE'%ch,fontsize=20)
    eq = r'$\frac{Rate~p_{{%s}_{angInt}}}{Rate~p_{0}}$'%ch
    plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
    # plt.grid(b=True, which='both', axis='both')

    savePath = os.path.join(desiredDir,'%s_RatioReactionRates_JAMESABOVE.png'%ch)
    plt.savefig(savePath,dpi=300)
    plt.clf()


    # df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    # df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)
    # df.to_excel('legendre_out/DATA/%s/a0/%s_rates_extrap_angIntFiner10.xlsx'%(ch,ch))
