import numpy as np
import pandas as pd
from scipy import special
import matplotlib.pyplot as plt
from AtomicMassTable import GetElement

"""
Calculate the reaction rates given the angle integrated cross section as
function of energy
"""


# GetElement returns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]
_u = 931.494    # MeV/c^2


# Read in the angle integrated cross-sections calculated via azure calculate w/o data
colNames = ['E_Cm','E_Ex','Angle_Cm',
                'Fit_Cm_Cross','Fit_Cm_S_factor']

dfp1 = pd.read_table('AZUREOut_aa=2_R=3.extrap',sep='\s+',names=colNames)
dfp2 = pd.read_table('AZUREOut_aa=2_R=4.extrap',sep='\s+',names=colNames)

fname = 'iliadis_Temps.dat' #rate_Temps
temperature = np.loadtxt( fname )

# Reduced mass
mu = _m24Mg*_m4He/(_m27Al+_m1H) * _u

dict_channels = {'p1':dfp1,'p2':dfp2}
channels = ['p1','p2']
for ch in channels:


    cross = dict_channels[ch]['Fit_Cm_Cross'].values
    E = dict_channels[ch]['E_Cm'].values     # E_lab
    E = E*(_m4He+_m24Mg)/_m24Mg                 # E_cm
    rxnRate = []
    dE = .001       # 1 keV steps
    for T in temperature:
        integrand = np.trapz(cross*np.exp(-11.604*E/T)*E,E)#*2*np.pi
        rate =  1E-24*6.022E23*3E10*(8/np.pi)**.5*(8.617**-1.5)* mu**(-.5) * T**(-1.5) * integrand
        rxnRate.append(rate)

    rxnRate = np.array(rxnRate)
    col = ["T9","Rate"]

    plt.scatter(temperature,rxnRate,c='coral')
    plt.plot(temperature,rxnRate,label='low extrap',color='coral')
    plt.yscale('log')
    # plt.xlim(0,10)
    # plt.ylim(1e-30,1e8)
    plt.ylim(1e-30,1e-1)
    plt.xlim(1,2)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
    plt.xlabel('Temperature (T9)')
    plt.title('%s Rates'%ch)

    # df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    # df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)

    df = pd.read_excel('legendre_out/DATA/%s/a0/%s_rates.xlsx'%(ch,ch),sheet_name='Sheet1',header=0)

    mytemps = df['T9'].to_numpy()
    myrates = df['Rate'].to_numpy()


    # My rates from the data
    plt.scatter(temperature,myrates,c='k')
    plt.plot(temperature,myrates,label='me',color='k')

    # Combined rates of my data plus azure extrap
    comb_rates = rxnRate+myrates
    plt.scatter(temperature,comb_rates,c='b')
    plt.plot(temperature,comb_rates,label='low extrap + me',color='b')


    plt.legend()
    # plt.savefig('images/%s_ReactionRate_azureAngleIntegrated.png'%ch,dpi=900)
    plt.show()
