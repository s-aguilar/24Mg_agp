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
    e_mask1 = ((E >= 3.42454) & (E <= 3.45806))
    e_mask2 = ((E >= 3.67652) & (E <= 4.66271))
    not_e_mask1 = (E <= 3.42454)
    not_e_mask2 = (E >= 4.66271)
    # Plot the angle-integrated AZURE2 extrapolated cross-section
    plt.clf()
    plt.plot(E[e_mask1],cross[e_mask1],c='b')
    plt.plot(E[e_mask2],cross[e_mask2],c='b',label='IN')
    plt.plot(E[not_e_mask1],cross[not_e_mask1],c='r')
    plt.plot(E[not_e_mask2],cross[not_e_mask2],c='r',label='OUT')
    # plt.scatter(E,cross,c='b')
    plt.xlabel('E_CM (MeV)',fontsize=14)
    plt.ylabel('Angle Integrated Cross-Section (barns)', fontsize=14)
    plt.yscale('log')
    plt.legend()
    plt.title('%s channel - Angle-Integrated Finer Excitation Curve' % ch,fontsize=20)
    plotPath = os.path.join(desiredDir,'%s_angIntExtrapFiner10ExCurve.png' % ch)
    plt.savefig(plotPath,dpi=300)
    # plt.show()
    plt.clf()


    rxnRateIN = []
    rxnRateOUT = []
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
        integralIN = np.trapz(integrand[e_mask1],E[e_mask1])
        integralIN += np.trapz(integrand[e_mask2],E[e_mask2])

        # colNames = ['splineX','splineY']
        # p1 = pd.read_excel('legendre_out/DATA/p1/p1IntegrandSplinePoints%f.xlsx'%(T),sheet_name='Sheet1',header=0)
        # p2 = pd.read_excel('legendre_out/DATA/p2/p2IntegrandSplinePoints%f.xlsx'%(T),sheet_name='Sheet1',header=0)
        #
        # p1intx = p1['splineX'].values
        # p1inty = p1['splineY'].values
        # p2intx = p2['splineX'].values
        # p2inty = p2['splineY'].values
        #
        # ee_mask1 = ((p1intx >= 3.42454) & (p1intx <= 3.45806))
        # ee_mask2 = ((p1intx >= 3.67652) & (p1intx <= 4.66271))
        #
        # plt.clf()
        # plt.plot(E[e_mask1],integrand[e_mask1],c='k')
        # plt.plot(E[e_mask2],integrand[e_mask2],c='k',label='azure p0')
        #
        # plt.plot(p1intx[ee_mask1],p1inty[ee_mask1],c='r')
        # plt.plot(p1intx[ee_mask2],p1inty[ee_mask2],c='r',label='data p1')
        #
        # plt.plot(p2intx[ee_mask1],p2inty[ee_mask1],c='g')
        # plt.plot(p2intx[ee_mask2],p2inty[ee_mask2],c='g',label='data p2')
        #
        #
        # plt.xlabel('E_CM (MeV)',fontsize=14)
        # plt.yscale('log')
        # plt.ylabel('Rxn Rate Integrand',fontsize=14)
        # plt.title('T = %fGK'%T)
        # plt.legend()
        # plt.savefig('integrandOverlay%f.png'%T,dpi=300)
        # filenames.append('integrandOverlay%f.png'%T)
        #
        #
        # plt.clf()


        _integralList.append(integralIN)


        rateIN =  rxnRateCONST * mu**(-.5) * T**(-1.5) * integralIN
        rxnRateIN.append(rateIN)

        integralOUT = np.trapz(integrand[not_e_mask1],E[not_e_mask1])
        integralOUT += np.trapz(integrand[not_e_mask2],E[not_e_mask2])

        # print(integralIN,integralOUT)
        # input()
        rateOUT =  rxnRateCONST * mu**(-.5) * T**(-1.5) * integralOUT
        rxnRateOUT.append(rateOUT)


    # images = []
    # for filename in filenames:
    #     images.append(imageio.imread(filename))
    # imageio.mimsave('integrandOverlay.gif', images,duration=1)

    # with imageio.get_writer('/path/to/movie.gif', mode='I') as writer:
    #     for filename in filenames:
    #         image = imageio.imread(filename)
    #         writer.append_data(image)

    print('Done with movie')


    # # Just plot this once
    # if ch == 'p0':
    #     # Plot the integral as fcn of E_cm Curve
    #     plt.clf()
    #     plt.plot(temperature,_integralList,c='b')
    #     plt.xlabel('Temperature (T9)',fontsize=14)
    #     plt.ylabel('Integral Term', fontsize=14)
    #     plt.yscale('log')
    #     plt.title('Integral Term vs T9',fontsize=20)
    #     plotPath = os.path.join(desiredDir,'integralTerm.png'%T)
    #     plt.savefig(plotPath,dpi=300)
    #     plt.clf()



    # rxnRate = np.array(rxnRate)
    #
    # col = ["T9","Rate"]
    #
    # myp0rates = pd.read_excel('legendre_out/DATA/p0/a0/p0_ratesEXTRAP.xlsx',sheet_name='Sheet1',header=0)
    # myp1rates = pd.read_excel('legendre_out/DATA/p1/a0/p1_ratesEXTRAP.xlsx',sheet_name='Sheet1',header=0)
    # myp2rates = pd.read_excel('legendre_out/DATA/p2/a0/p2_ratesEXTRAP.xlsx',sheet_name='Sheet1',header=0)
    #
    # azp0rates = pd.read_table('p0rates.out',sep='\s+')
    # azp1rates = pd.read_table('p1rates.out',sep='\s+')
    # azp2rates = pd.read_table('p2rates.out',sep='\s+')
    #
    # # print(myp1rates.columns)
    # # print(azp1rates.columns)
    #
    # temp0 = myp0rates['T9'].values
    # temp1 = myp1rates['T9'].values
    # temp2 = myp2rates['T9'].values
    # aztemp0 = azp0rates['T9'].values
    # aztemp1 = azp1rates['T9'].values
    # aztemp2 = azp2rates['T9'].values
    #
    # rate0 = myp0rates['Rate'].values
    # rate1 = myp1rates['Rate'].values
    # rate2 = myp2rates['Rate'].values
    # azrate0 = azp0rates['Rate'].values
    # azrate1 = azp1rates['Rate'].values
    # azrate2 = azp2rates['Rate'].values
    #
    # # plt.scatter(temp0,rate0,c='silver')
    # # plt.plot(temp0,rate0,color='silver',label='az extrap p0')
    #
    # if ch == 'p1':
    #     # plt.scatter(temp1,rate1,c='k')
    #     plt.plot(temp1,rate1,color='k',label='az extrap p1')
    #
    # if ch == 'p2':
    #     # plt.scatter(temp2,rate2,c='grey')
    #     plt.plot(temp2,rate2,color='grey',label='az extrap p2')
    #
    #
    # if ch == 'p0':
    #     # plt.scatter(aztemp0,azrate0,c='g')
    #     plt.plot(aztemp0,azrate0,color='g',label='azure rate p0')
    #
    # if ch == 'p1':
    #     # plt.scatter(aztemp1,azrate1,c='fuchsia')
    #     plt.plot(aztemp1,azrate1,color='fuchsia',label='azure rate p1')
    #
    # if ch == 'p2':
    #     # plt.scatter(aztemp2,azrate2,c='b')
    #     plt.plot(aztemp2,azrate2,color='b',label='azure rate p2')


    plt.plot(temperature,rxnRateIN,label='%s-angInt-IN'%ch,color='b')
    # plt.legend()
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-20,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('AZURE2 IN %s - Reaction Rate'%ch,fontsize=20)
    plt.savefig('%s_ReactionRate_azureAngleIntegratedFiner10IN.png'%ch,dpi=900)

    plt.plot(temperature,rxnRateOUT,label='%s-angInt-OUT'%ch,color='r')
    # plt.legend()
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-20,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('AZURE2 OUT %s - Reaction Rate'%ch,fontsize=20)
    plt.legend()
    plt.savefig('%s_ReactionRate_azureAngleIntegratedFiner10OUT.png'%ch,dpi=900)






    # Plot ratio of azures angle integrated extrapolated cross-section calculated reaction rate for inner and outer
    # vs Iliadis reaction rate

    pathToREACLIB = os.path.join(desiredDir,'24MgREACLIB.xlsx')
    df = pd.read_excel(pathToREACLIB,sheet_name='il10',header=0)
    _tempIliadis = df['T9'].values
    _rateIliadis = df['Rate'].values

    plt.clf()
    y = rxnRateIN/_rateIliadis
    plt.plot(_tempIliadis,y,color='k')
    plt.xlim(0,10)
    yy = [1,1]
    xx = np.linspace(0,10,2)
    plt.plot(xx,yy)
    # plt.ylim(.9,1.)
    plt.ylabel('Ratio of Reaction Rate',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('p$_{%s}$ Reaction Rates Ratio - IN'%ch,fontsize=20)
    eq = r'$\frac{Rate~p_{{%s}_{angInt}}}{Rate~p_{0}}$'%ch
    plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
    # plt.grid(b=True, which='both', axis='both')

    savePath = os.path.join(desiredDir,'%s_RatioReactionRates_IN.png'%ch)
    plt.savefig(savePath,dpi=300)
    plt.clf()


    y = rxnRateOUT/_rateIliadis
    plt.plot(_tempIliadis,y,color='k')
    plt.xlim(0,10)
    yy = [1,1]
    xx = np.linspace(0,10,2)
    plt.plot(xx,yy)
    # plt.ylim(.9,1.)
    plt.ylabel('Ratio of Reaction Rate',fontsize=14)
    plt.xlabel('Temperature (T9)',fontsize=14)
    plt.title('p$_{%s}$ Reaction Rates Ratio - OUT'%ch,fontsize=20)
    eq = r'$\frac{Rate~p_{{%s}_{angInt}}}{Rate~p_{0}}$'%ch
    plt.text(9, 1, eq, {'color': 'k', 'fontsize': 18}, va="top", ha="right")
    # plt.grid(b=True, which='both', axis='both')

    savePath = os.path.join(desiredDir,'%s_RatioReactionRates_OUT.png'%ch)
    plt.savefig(savePath,dpi=300)
    plt.clf()


    # df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    # df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)
    # df.to_excel('legendre_out/DATA/%s/a0/%s_rates_extrap_angIntFiner10.xlsx'%(ch,ch))
