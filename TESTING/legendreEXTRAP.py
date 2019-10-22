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
from scipy import integrate


"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.

The 'plot_a' variable enables/disables the plots demonstrating the different
order legendre polynomial fits for 'a' analytic fit
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
# test = 3.7313E10    # constant from literature
# BOTH ARE EQUAL

# print(rxnRateCONST)
# print(test)
# exit()




def convert(old):
    """
    Convert the 'old' legendre coefficient outputs (only the even terms) into a
    version which can then be used by the 'np.polynomial.legendre.legval()'
    function (Both even and odd, so pad the odd terms with 0's)
    """

    # loop through the set of coefficients appending 0's
    new = []
    for _ in old:
        new.append(_)
        new.append(0)

    # Pop the last 0 that is unnecessary
    new.pop()

    return np.array(new)


# Read in the data into dataframe

# labels for azure fit
# colNames = ['E_Cm','E_Ex','Angle_Cm',
#                 'Fit_Cm_Cross','Fit_Cm_S_factor',
#                 'Data_Cm_Cross','Data_Cm_Cross_Unc',
#                 'Data_Cm_S_factor','Data_Cm_S_factor_Unc']

#labels for azure extrapolate segment without data
colNames = ['E_Cm','E_Ex','Angle_Cm',
                'Fit_Cm_Cross','Fit_Cm_S_factor']
p0 = pd.read_table('rMatrix/AZUREOut_aa=2_R=1.extrap',sep='\s+',names=colNames)
p1 = pd.read_table('rMatrix/AZUREOut_aa=2_R=3.extrap',sep='\s+',names=colNames)
p2 = pd.read_table('rMatrix/AZUREOut_aa=2_R=4.extrap',sep='\s+',names=colNames)



# 'seg' is just the number of slices/indices in the azure2 '.extrap' file
# 'unc' is just an artificial uncertainty that is 5% of the measurement, need
# an uncertainty to analytically solve for the parameters


# Put angle in lab and Energy in lab E
seg = 4100
ang0 = [0]*seg
ang15 = [15]*seg
ang30 = [30]*seg
ang45 = [45]*seg
ang60 = [60]*seg
ang75 = [75]*seg
ang90 = [90]*seg
ang = np.array(ang0+ang15+ang30+ang45+ang60+ang75+ang90,dtype=np.int64)
unc = p0['Fit_Cm_Cross'].values*.05
# print(len(p0.index))
p0 = p0.assign(Angle=pd.Series(ang,index=p0.index).values)
p0 = p0.assign(Data_Cm_Cross_Unc=pd.Series(unc,index=p0.index).values)


# Put angle in lab and Energy in lab E
seg = 3100*2+1
ang0 = [0]*seg
ang15 = [15]*seg
ang30 = [30]*seg
ang45 = [45]*seg
ang60 = [60]*seg
ang75 = [75]*seg
ang90 = [90]*seg
ang = np.array(ang0+ang15+ang30+ang45+ang60+ang75+ang90,dtype=np.int64)
unc = p1['Fit_Cm_Cross'].values*.05
# print(len(p1.index))
p1 = p1.assign(Angle=pd.Series(ang,index=p1.index).values)
p1 = p1.assign(Data_Cm_Cross_Unc=pd.Series(unc,index=p1.index).values)


# Put angle in lab and Energy in lab E
seg = 2800*2+1
ang0 = [0]*seg
ang15 = [15]*seg
ang30 = [30]*seg
ang45 = [45]*seg
ang60 = [60]*seg
ang75 = [75]*seg
ang90 = [90]*seg
ang = np.array(ang0+ang15+ang30+ang45+ang60+ang75+ang90,dtype=np.int64)
unc = p2['Fit_Cm_Cross'].values*.05
p2 = p2.assign(Angle=pd.Series(ang,index=p2.index).values)
p2 = p2.assign(Data_Cm_Cross_Unc=pd.Series(unc,index=p2.index).values)


# Conversion of E_Cm to E_lab
sep_alpha = 9.984   # Separation energy of alpha in MeV

energy = p0['E_Cm'].values
# E_Lab = energy*(_m4He+_m24Mg)/_m24Mg
E_Lab = energy
p0 = p0.assign(Energy=pd.Series(E_Lab,index=p0.index).values)

energy = p1['E_Cm'].values
# E_Lab = energy*(_m4He+_m24Mg)/_m24Mg
E_Lab = energy
p1 = p1.assign(Energy=pd.Series(E_Lab,index=p1.index).values)

energy = p2['E_Cm'].values
# E_Lab = energy*(_m4He+_m24Mg)/_m24Mg
E_Lab = energy
p2 = p2.assign(Energy=pd.Series(E_Lab,index=p2.index).values)
# print(p1.head())


dict_channels = {'p0':p0,'p1':p1,'p2':p2}


# If plot is 0 then no plotting, any other value it will generate the plots
# 'a' for analytic solution plots
plot_a = 1


# Perform the analysis over all the channels
channels = ['p0','p1','p2']
for ch in channels:

    print(ch,' channel:\n-Working on Legendre fitting')
    angle = np.array([0,15,30,45,60,75,90])

    chan = dict_channels[ch]

    energy_chan = chan['Energy'].values
    angle_chan = chan['Angle'].values
    cross_chan = chan['Fit_Cm_Cross'].values
    crossErr_chan = chan['Data_Cm_Cross_Unc'].values


    legendre_order = [[0],[0,2],[0,2,4],[0,2,4,6],[0,2,4,6,8],[0,2,4,6,8,10]]
    color_order = ['b','g','r','c','m','y']


    # Lists for python minimization solutions
    ord_0 = []
    ord_2 = []
    ord_4 = []
    ord_6 = []
    ord_8 = []
    ord_10 = []


    # Dict for python minimization solutions
    dict_ord = {'0':ord_0,'2':ord_2,'4':ord_4,'6':ord_6,'8':ord_8,'10':ord_10}


    # Lists for analytic solutions
    ord_00,ord_02,ord_04,ord_06 = [],[],[],[]
    ord_err_00,ord_err_02,ord_err_04,ord_err_06 = [],[],[],[]
    ord_00_x2,ord_02_x2,ord_04_x2,ord_06_x2 = [],[],[],[]
    ord_00_x2ndf,ord_02_x2ndf,ord_04_x2ndf,ord_06_x2ndf = [],[],[],[]
    ord_00_pVal,ord_02_pVal,ord_04_pVal,ord_06_pVal = [],[],[],[]


    # Dicts for analytic solutions
    dict_ord_a = {'0':ord_00,'2':ord_02,'4':ord_04,'6':ord_06}
    dict_ord_err_a = {'0':ord_err_00,'2':ord_err_02,'4':ord_err_04,'6':ord_err_06}
    dict_ord_x2_a = {'0':ord_00_x2,'2':ord_02_x2,'4':ord_04_x2,'6':ord_06_x2}
    dict_ord_x2ndf_a = {'0':ord_00_x2ndf,'2':ord_02_x2ndf,'4':ord_04_x2ndf,'6':ord_06_x2ndf}
    dict_ord_pVal_a = {'0':ord_00_pVal,'2':ord_02_pVal,'4':ord_04_pVal,'6':ord_06_pVal}


    # Passing list to a set does not preserve the original order (ascending) so
    # manually sort it back into a new list
    temp_nrg = set()
    energyList = np.array([x for x in energy_chan if not (x in temp_nrg or temp_nrg.add(x))],dtype=np.float64)


    # perform Legendre fit over the angular distribution at each energy
    for nrg in energyList:
        is_nrg = chan['Energy']==nrg      # create boolean mask
        chan_temp = chan[is_nrg]

        masked_nrg_chan = chan_temp['Energy'].values
        masked_ang_chan = np.radians(chan_temp['Angle'].values)
        masked_cross_chan = chan_temp['Fit_Cm_Cross'].values
        masked_err_chan = chan_temp['Data_Cm_Cross_Unc'].values

        temp_weights = [1/_**2 for _ in masked_err_chan]
        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED


        # Up to Legendre Polynomial "BLAH" in legendre_order[list] index starts from 0, can go up to 5
        leg_ord = 2


        # Loop through the number of terms used in Legendre expansion:
        # -- a0, a0+a2, ... , a0+a2+a4+a6+a8+a10
        for ind in range(0,leg_ord+1):
            results = chi2_mat(masked_cross_chan,masked_err_chan,ind)

            # TURN ON/OFF if you want to debug values
            # stop = input('...Pause...\n\n')


            dict_ord_a[str(legendre_order[ind][-1])].append(results[0])
            dict_ord_err_a[str(legendre_order[ind][-1])].append(results[1])
            dict_ord_x2_a[str(legendre_order[ind][-1])].append(results[2])
            dict_ord_x2ndf_a[str(legendre_order[ind][-1])].append(results[3])
            dict_ord_pVal_a[str(legendre_order[ind][-1])].append(results[4])

            if plot_a:
                # Note: The legendre.legval function requires the coefficients
                # of all terms up to whatever order you want, must pad the even
                # terms with 0's
                temp_coef = convert(results[0])
                temp_coef_err = convert(results[1])

        # # Angular distribution plots
        # if  plot_a:
        #     plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
        #     plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
        #     plt.title('%s - %5.4f MeV - Angular Distribution' % (ch,nrg))
        #     plt.legend()
        #     plt.savefig('legendre_out/fits/%s/%5.4fMeVFit.png'%(ch,nrg),dpi=100)
        #     plt.clf()


    # Convert lists into arrays
    for ind in range(0,leg_ord+1):
        dict_ord_a[str(legendre_order[ind][-1])] = np.array(dict_ord_a[str(legendre_order[ind][-1])])
        dict_ord_err_a[str(legendre_order[ind][-1])] = np.array(dict_ord_err_a[str(legendre_order[ind][-1])])


    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    # colLabels = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','chi2','p-value','chi2ndf']
    colLabels = ['a0','a2','a4','a6','a8','a10']
    colLabelss = ['a0_err','a2_err','a4_err','a6_err','a8_err','a10_err']
    for _ in range(0,leg_ord+1):
        stuff = dict_ord_a[str(legendre_order[_][-1])]
        stuff2 = dict_ord_err_a[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_+1])
        df1 =  pd.DataFrame(data=stuff2,index=energyList,columns=colLabelss[0:_+1])
        df = pd.concat([df, df1], axis=1, sort=False)
        df = df.assign(Chi2=pd.Series(dict_ord_x2_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Pval=pd.Series(dict_ord_pVal_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Chi2NDF=pd.Series(dict_ord_x2ndf_a[str(legendre_order[_][-1])],index=df.index).values)
        df.to_csv('legendre_out/DATA/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        df.to_excel('legendre_out/DATA/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))


    # Now I have all the 'a_i' coefficients for all the many order fits, need
    # to select/filter the necessary 'a0' terms from the fits. This will be
    # done through the p-value.


    # Read in all my a0s and p-values and compare each energies fit side by side
    # starting from lowest order to highest, select order which pvalue is maximal,
    # append the a0 term + continue and repeat for the next energy step

    a0_final = []
    a0_err_final = []
    p_final = []
    x2_final = []
    order = []


    for ind in range(len(ord_00)):
        a0_final.append(dict_ord_a[str(legendre_order[2][-1])][ind][0])
        a0_err_final.append(dict_ord_err_a[str(legendre_order[2][-1])][ind][0])
        p_final.append(dict_ord_pVal_a[str(legendre_order[2][-1])][ind])
        order.append(str(legendre_order[2][-1]))


    # Now spline the 'a0' coefficients and plot as function of energy and overlay
    # the original 'a0' coefficents and make sure the splining is done well.

    solidAngle = 4*np.pi
    # solidAngle = 1
    a0_final = solidAngle*np.array(a0_final)
    a0_err_final = solidAngle*np.array(a0_err_final)


    spline1 = CubicSpline(energyList,a0_final)
    y_val1 = spline1(energyList)


    # # Plot the splines and the points
    # plt.clf()
    # plt.plot(energyList,y_val1,c='b')
    # plt.scatter(energyList,a0_final)
    # plt.yscale('log')
    # plt.title('%s channel - a0 vs E - Cross-section' % ch)
    # plt.savefig('legendre_out/excitationCurve/%s/%s_a0Curve.png' % (ch,ch),dpi=900)
    # plt.clf()


    # Save the points, first merge the two splines
    df = pd.DataFrame(data=a0_final,index=energyList,columns=['a0'])
    df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).values)
    df = df.assign(fitOrder=pd.Series(order,index=df.index).values)
    df = df.assign(splineX=pd.Series(energyList,index=df.index).values)
    df = df.assign(splineY=pd.Series(y_val1,index=df.index).values)
    df.to_excel('legendre_out/DATA/%s/%sFitPoints.xlsx'%(ch,ch))


    # Calculate the reaction rates
    fname = 'rate_Temps.dat'
    temperature = np.loadtxt( fname )

    mu = _m24Mg*_m4He/(_m27Al+_m1H) #* _u

    if ch == 'a1':
        mu = _m24Mg*_m4He/(_m24Mg+_m4He) #* _u   # reduced mass


    # rxnRate = []
    # rxnRate1 = []
    rxnRate2 = []
    rxnRate3 = []
    rxnRate4 = []
    rxnRate5 = []
    rxnRate6 = []
    rxnRate7 = []

    dE = .001
    # dE = .0005 # half keV steps
    for T in temperature:

        # Convert Lab E to E_Cm
        E1 = energyList
        # E1 = energyList*_m24Mg/(_m24Mg+_m4He)
        # print(E1)
        integrand1 = CubicSpline(E1,a0_final*np.exp(-11.604*E1/T)*E1)

        # print(integrand1(E1))
        # plt.plot(E1,integrand1(E1),color='r')
        # plt.scatter(E1,a0_final*np.exp(-11.604*E1/T)*E1,alpha=.01,c='b')
        # plt.yscale('log')
        # plt.show()

        # intg = integrate.romberg(integrand1,min(energyList),max(energyList),divmax=40)
        # intg1 = integrate.quad(integrand1,min(energyList),max(energyList),limit=200)
        intg2 = 0
        for ind,elem in enumerate(E1):
            intg2 +=  a0_final[ind]*np.exp(-11.604*E1[ind]/T)*E1[ind]*dE

        intg3 = integrand1.integrate(min(E1),max(E1))
        intg4 = integrate.cumtrapz(a0_final*np.exp(-11.604*E1/T)*E1,E1)[-1]
        intg5 = integrate.trapz(a0_final*np.exp(-11.604*E1/T)*E1,E1)
        intg6 = integrate.simps(a0_final*np.exp(-11.604*E1/T)*E1,E1)
        intg7 = np.trapz(a0_final*np.exp(-11.604*E1/T)*E1,E1)



        # rate =  rxnRateCONST * T**(-1.5) * intg
        # rate1 = rxnRateCONST * T**(-1.5) * intg1[0]
        rate2 = rxnRateCONST * T**(-1.5) * intg2
        rate3 = rxnRateCONST * T**(-1.5) * intg3
        rate4 = rxnRateCONST * T**(-1.5) * intg4
        rate5 = rxnRateCONST * T**(-1.5) * intg5
        rate6 = rxnRateCONST * T**(-1.5) * intg6
        rate7 = rxnRateCONST * T**(-1.5) * intg7

        # rxnRate.append(rate)
        # rxnRate1.append(rate1)
        rxnRate2.append(rate2)
        rxnRate3.append(rate3)
        rxnRate4.append(rate4)
        rxnRate5.append(rate5)
        rxnRate6.append(rate6)
        rxnRate7.append(rate7)

    # rxnRate = np.array(rxnRate)
    # rxnRate1 = np.array(rxnRate1)
    rxnRate2 = np.array(rxnRate2)
    rxnRate3 = np.array(rxnRate3)
    rxnRate4 = np.array(rxnRate4)
    rxnRate5 = np.array(rxnRate5)
    rxnRate6 = np.array(rxnRate6)
    rxnRate7 = np.array(rxnRate7)

    # plt.scatter(temperature,rxnRate,c='b')
    # plt.plot(temperature,rxnRate,label='romberg',color='b')
    # plt.scatter(temperature,rxnRate1,c='g')
    # plt.plot(temperature,rxnRate1,label='quad',color='g')
    plt.scatter(temperature,rxnRate2,c='r')
    plt.plot(temperature,rxnRate2,label='Sum',color='r')
    plt.scatter(temperature,rxnRate3,c='y')
    plt.plot(temperature,rxnRate3,label='Integrate spl',color='y')
    plt.scatter(temperature,rxnRate4,c='grey')
    plt.plot(temperature,rxnRate4,label='cumtrapz',color='grey')
    plt.scatter(temperature,rxnRate5,c='maroon')
    plt.plot(temperature,rxnRate5,label='trapz',color='maroon')
    plt.scatter(temperature,rxnRate6,c='coral')
    plt.plot(temperature,rxnRate6,label='np trapz',color='coral')
    plt.scatter(temperature,rxnRate7,c='fuchsia')
    plt.plot(temperature,rxnRate7,label='simps',color='fuchsia')
    plt.legend()
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-30,1e8)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
    plt.xlabel('Temperature (T9)')
    plt.savefig('%s_ReactionRate.png'%ch,dpi=900)
    plt.clf()
    df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    df = df.assign(Rate=pd.Series(rxnRate6,index=df.index).values) # np.trapz
    df.to_excel('legendre_out/DATA/%s/a%d/%s_rates.xlsx'%(ch,legendre_order[0][-1],ch))

    print('\n\n')
# """

print('DONE!')
