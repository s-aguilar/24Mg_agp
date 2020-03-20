import os
import openpyxl # conda install -c anaconda openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from AtomicMassTable import GetElement
from chi2 import chi2_det,chi2_mat
from scipy.interpolate import LSQUnivariateSpline
from scipy import special
from scipy import interpolate
from scipy import constants
from scipy import integrate
import pdb

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12

# Use to switch directory paths to use (for the main or TESTING)
import os
currentDir = os.path.realpath('')
parentDir = os.path.realpath('..')



"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir


"""
ALWAYS SET THIS PATH TO THE ONE YOU WANT TO USE!
"""
desiredDir = currentDir

"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.

The 'plot_a' variable enables/disables the plots demonstrating the different
order legendre polynomial fits 'a' analytic fit

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
colNames = ['Energy','Angle','Cross-section','Error']
o17_1 = pd.read_table('rMatrix/24Mg_rMatrix_o17_1.dat',names=colNames)


dict_channels = {'o17_1':o17_1}

# If plot is 0 then no plotting, any other value it will generate the plots
# 'a' for analytic solution plots
plot_a = 0


# Perform the analysis over all the channels
channels = ['o17_1']

for ch in channels:

    print(ch,' channel:\n-Working on Legendre fitting')
    angle = np.array([0,15,30,45,60,75,90])

    chan = dict_channels[ch]

    # is in Lab energy, convert to center-of-mass
    energyCM_chan = chan['Energy'].to_numpy()#LAB    *(_m24Mg/(_m24Mg+_m4He))   # Now its in E_cm
    chan = chan.assign(E_CM=pd.Series(energyCM_chan,index=chan.index).to_numpy())


    # Mask low stats fits
    maskk = ( (chan['Cross-section']>1e-13) & (chan['Error']<chan['Cross-section']) &  (chan['E_CM'] != 4.3172))
    print(len(chan))
    chan = chan[maskk]
    print(len(chan))
    # exit()


    legendre_order = [[0],[0,2],[0,2,4],[0,2,4,6],[0,2,4,6,8],[0,2,4,6,8,10]]
    color_order = ['b','g','r','c','m','y']


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
    energyList = np.array([x for x in chan['E_CM'] if not (x in temp_nrg or temp_nrg.add(x))],dtype=np.float64)

    # Sort by energy, keeping others consistent!
    ind = energyList.argsort()
    energyList = energyList[ind]

    '''
    Degenerate energy point at 5.275 MeV (manually deleted one of them)
    and at 5.3299
    '''


    # perform Legendre fit over the angular distribution at each energy
    for nrg in energyList:
        is_nrg = chan['E_CM']==nrg      # create boolean mask
        chan_temp = chan[is_nrg]

        print(nrg)

        # Energy is in center-of-mass
        masked_nrg_chan = chan_temp['E_CM'].to_numpy()
        masked_ang_chan = np.radians(chan_temp['Angle'].to_numpy())
        masked_cross_chan = chan_temp['Cross-section'].to_numpy()
        masked_err_chan = chan_temp['Error'].to_numpy()

        temp_weights = [1/_**2 for _ in masked_err_chan]
        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED


        # Up to Legendre Polynomial "BLAH" in legendre_order[list] index starts from 0, can go up to 5
        leg_ord = 2


        # Loop through the number of terms used in Legendre expansion:
        # -- a0, a0+a2, ... , a0+a2+a4+a6+a8+a10
        for ind in range(0,leg_ord+1):
            # print(nrg)
            results = chi2_mat(masked_cross_chan,masked_err_chan,masked_ang_chan,ind)

            dict_ord_a[str(legendre_order[ind][-1])].append(results[0])#*37/1e6 * 1/23.985*6.022e23*1/1e-24)
            dict_ord_err_a[str(legendre_order[ind][-1])].append(results[1])#*37/1e6 * 1/23.985*6.022e23*1/1e-24)
            dict_ord_x2_a[str(legendre_order[ind][-1])].append(results[2])
            dict_ord_x2ndf_a[str(legendre_order[ind][-1])].append(results[3])
            dict_ord_pVal_a[str(legendre_order[ind][-1])].append(results[4])

            if plot_a:
                # Note: The legendre.legval function requires the coefficients
                # of all terms up to whatever order you want, must pad the even
                # terms with 0's, this is done through the 'convert' function
                temp_coef = convert(results[0])
                temp_coef_err = convert(results[1])


                # Scatter has bug when the y-values are small, the y-axis does not autoscale
                # plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])
                if ind ==2:
                    plt.plot(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],label='%d$^{th}$ Order'%legendre_order[ind][-1])

        # Angular distribution plots
        if plot_a:
            plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
            plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
            plt.xlabel('Angle',fontsize=14)
            plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
            plt.title('$^{24}$Mg($\\alpha$,p$_{%s}\\gamma$) - %5.4f MeV' % (ch[-1],nrg),fontsize=20)
            plt.legend(prop={'size': 14})
            plt.savefig('legendre_out/fits/%s/%5.4fMeVFit.png'%(ch,nrg),dpi=300)
            plt.clf()


    leg_ord = 2
    # Convert lists into arrays
    for ind in range(0,leg_ord+1):
        dict_ord_a[str(legendre_order[ind][-1])] = np.array(dict_ord_a[str(legendre_order[ind][-1])])
        dict_ord_err_a[str(legendre_order[ind][-1])] = np.array(dict_ord_err_a[str(legendre_order[ind][-1])])


    # Now we have all the coefficients for each order fit, lets plot them
    # plot 'a' coefficients for each order as function of energy as well as
    # the chi2ndf
    ord_step = 0
    print('-Working on plotting the coefficients')
    for ind in range(0,leg_ord+1):
        for ord_step in range(0,leg_ord+1):
            plt.clf()
            if ind > ord_step: continue
            a = [ (_[ind],) for _ in dict_ord_a[str(legendre_order[ord_step][-1])]]     # analytic solution
            plt.scatter(energyList,a,s=8)
            plt.title('%s channel - Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1],legendre_order[ord_step][-1]))
            plt.savefig('legendre_out/coef_curve/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=100)
            ord_step += 1

        plt.clf()

        x2 = dict_ord_x2_a[str(legendre_order[ind][-1])]
        x2ndf = dict_ord_x2ndf_a[str(legendre_order[ind][-1])]
        # print(str(legendre_order[ind][-1]))
        # plt.scatter(energyList,x2ndf,s=8)
        # plt.title('%s channel - $\chi^{2}$ for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1]))
        # plt.savefig('legendre_out/chi2/%s/a%d/chi2.png'%(ch,legendre_order[ind][-1]),dpi=100)


    plt.clf()

    # """
    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    colLabels = ['a0','a2','a4','a6','a8','a10']
    colLabelss = ['a0_err','a2_err','a4_err','a6_err','a8_err','a10_err']
    for _ in range(0,leg_ord+1):
        stuff = dict_ord_a[str(legendre_order[_][-1])]
        stuff2 = dict_ord_err_a[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_+1])
        df1 =  pd.DataFrame(data=stuff2,index=energyList,columns=colLabelss[0:_+1])
        df = pd.concat([df, df1], axis=1, sort=False)
        df = df.assign(Chi2=pd.Series(dict_ord_x2_a[str(legendre_order[_][-1])],index=df.index).to_numpy())
        df = df.assign(Pval=pd.Series(dict_ord_pVal_a[str(legendre_order[_][-1])],index=df.index).to_numpy())
        df = df.assign(Chi2NDF=pd.Series(dict_ord_x2ndf_a[str(legendre_order[_][-1])],index=df.index).to_numpy())
        df.to_csv('legendre_out/DATA/%s/a%d/%s_a%dFit.csv'%(ch,legendre_order[_][-1],ch,legendre_order[_][-1]))
        df.to_excel('legendre_out/DATA/%s/a%d/%s_a%dFit.xlsx'%(ch,legendre_order[_][-1],ch,legendre_order[_][-1]))
    # """

    # Now I have all the 'a_i' coefficients for all the many order fits, need
    # to select/filter the necessary 'a0' terms from the fits. This will be
    # done through the p-value.
    # Read in all my a0s and p-values and compare each energies fit side by side
    # starting from lowest order to highest, select order which pvalue is maximal,
    # append the a0 term + continue and repeat for the next energy step
    # (NOT USED ANYMORE)

    a0_final = []
    a0_err_final = []
    p_final = []
    x2_final = []
    order = []

    # Only keep the a0s from the 4th order legendre polynomial fit
    for ind in range(len(ord_00)):
        a0_final.append(dict_ord_a[str(legendre_order[2][-1])][ind][0])
        a0_err_final.append(dict_ord_err_a[str(legendre_order[2][-1])][ind][0])
        p_final.append(dict_ord_pVal_a[str(legendre_order[2][-1])][ind])
        order.append(str(legendre_order[2][-1]))


    # Save the points
    df = pd.DataFrame(data=a0_final,index=energyList)
    df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).to_numpy())
    df = df.assign(Pval=pd.Series(p_final,index=df.index).to_numpy())
    df = df.assign(Order=pd.Series(order,index=df.index).to_numpy())
    df.to_excel('legendre_out/DATA/%s/%sFits.xlsx'%(ch,ch))

    # Save angle integrated cross section for rMatrix
    with open('rMatrix/24Mg_rMatrix_%s_angInt.dat'%ch,'w') as f:
        print(len(a0_final),len(a0_err_final),len(energyList))

        cross = np.array(a0_final)*4*np.pi
        cross_err = np.array(a0_err_final)*4*np.pi
        print(len(a0_final),len(a0_err_final),len(energyList))
        # exit()
        for loop in range(len(energyList)):
            if cross[loop] > 0 and cross_err[loop] > 0 and cross_err[loop] < cross[loop]:
                printOut= '%f \t %d \t %.8E \t %.8E \n' %(energyList[loop],0,cross[loop],cross_err[loop])
                f.write(printOut)
            else:
                print('Problem at Ep:',energyList[loop])


    # Now spline the 'a0' coefficients and plot as function of energy and overlay
    # the original 'a0' coefficents and make sure the splining is done well.

    solidAngle = 4*np.pi
    a0_final = np.array(a0_final) * solidAngle
    a0_err_final = np.array(a0_err_final) * solidAngle


    # Sort by energy, keeping others consistent!
    ind = energyList.argsort()
    energyList = energyList[ind]
    a0_final = a0_final[ind]
    a0_err_final = a0_err_final[ind]

    mask_ = ((a0_final>a0_err_final))
    energyList = energyList[mask_]
    a0_final = a0_final[mask_]
    a0_err_final = a0_err_final[mask_]

    plt.errorbar(energyList,a0_final,yerr=a0_err_final,c='k',fmt='.')
    plt.plot(energyList,a0_final,c='b')
    plt.ylabel('Yield (arb. units)')
    plt.xlabel('Lab Energy (MeV)')
    plt.title('Yield Curve')
    plt.savefig('O17_1_YieldCurve.png')
    # plt.show()



    with open('O17_1_angIt_YieldCurve.dat','w') as f:
        for loop in range(len(energyList)):
            printOut= '%f \t %d \t %.8E \t %.8E \n' %(energyList[loop],0,a0_final[loop],a0_err_final[loop])
            f.write(printOut)
