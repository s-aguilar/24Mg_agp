import os
import openpyxl # conda install -c anaconda openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from AtomicMassTable import GetElement
from chi2 import chi2_det,chi2_mat
from scipy.interpolate import CubicSpline
from scipy import special
from scipy import constants


"""
Read in AZURE outputs, plot the R-Matrix fit
"""

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# GetElement resturns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]



# Read in the data into dataframe
colNames = ['E_Cm','E_Ex','Angle_Cm',
                'Fit_Cm_Cross','Fit_Cm_S_factor',
                'Data_Cm_Cross','Data_Cm_Cross_Unc',
                'Data_Cm_S_factor','Data_Cm_S_factor_Unc']
p1 = pd.read_table('rMatrix/AZUREOut_aa=2_R=3.out',sep='\s+',names=colNames)
p2 = pd.read_table('rMatrix/AZUREOut_aa=2_R=4.out',sep='\s+',names=colNames)
# print(p1.head())


# Put angle in lab and Energy in lab E
numOfDataPoints = 228
ang0 = [0]*numOfDataPoints
ang15 = [15]*numOfDataPoints
ang30 = [30]*numOfDataPoints
ang45 = [45]*numOfDataPoints
ang60 = [60]*numOfDataPoints
ang75 = [75]*numOfDataPoints
ang90 = [90]*numOfDataPoints

ang = np.array(ang0+ang15+ang30+ang45+ang60+ang75+ang90,dtype=np.int64)

p1 = p1.assign(Angle=pd.Series(ang,index=p1.index).values)
p2 = p2.assign(Angle=pd.Series(ang,index=p2.index).values)


# Conversion of E_Cm to E_lab
sep_alpha = 9.984   # Separation energy of alpha in MeV
energy = p1['E_Cm'].values
E_Lab = energy*(_m4He+_m24Mg)/_m24Mg
# E_Lab = energy


p1 = p1.assign(Energy=pd.Series(E_Lab,index=p1.index).values)
p2 = p2.assign(Energy=pd.Series(E_Lab,index=p2.index).values)
# print(p1.head())


dict_channels = {'p1':p1,'p2':p2}



# Perform the analysis over all the channels
channels = ['p1','p2']
for ch in channels:

    angle = np.array([0,15,30,45,60,75,90])

    chan = dict_channels[ch]

    energy_chan = chan['Energy'].values     # Lab Energy
    angle_chan = chan['Angle'].values
    cross_chan = chan['Fit_Cm_Cross'].values
    data_chan = chan['Data_Cm_Cross'].values
    dataErr_chan = chan['Data_Cm_Cross_Unc'].values



    for ang in angle:
        mask = angle_chan == ang

        plt.clf()
        plt.errorbar(energy_chan[mask],data_chan[mask],yerr=dataErr_chan[mask],c='b',fmt='.',markersize='2')
        plt.plot(energy_chan[mask],cross_chan[mask],color='r',label='R-Matrix')
        plt.yscale('log')
        plt.ylim(1e-6,1)
        plt.xlim(4,5.6)
        plt.xlabel('$E_{\\alpha}$ (MeV)', fontsize=14)
        plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
        plt.title('$^{24}$Mg($\\alpha$,p$_{%s}\\gamma$) - %d$^{\circ}$'%(ch[-1],ang), fontsize=20)
        plt.legend(prop={'size': 14})
        plt.savefig('crossSection/P%s/RMatrix_p%s_%d.png'%(ch[-1],ch[-1],ang),dpi=600)
        # plt.show()

        plt.clf()
