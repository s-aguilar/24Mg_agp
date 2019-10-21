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

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.

The 'plot_x' variable enables/disables the plots demonstrating the different
order legendre polynomial fits for either 'a' analytic or 'n' numeric fit

There is a 'stop' variable that can be turned off, it is used for
debugging purposes.
"""

plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12


# GetElement resturns (index,Z,A,Mass,'element')
_m1H = GetElement(1,1)[3]
_m4He = GetElement(2,4)[3]
_m24Mg = GetElement(12,24)[3]
_m27Al = GetElement(13,27)[3]



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
p1 = pd.read_table('rMatrix/rMatrix_p1s.dat',names=colNames)
p2 = pd.read_table('rMatrix/rMatrix_p2s.dat',names=colNames)
a1 = pd.read_table('rMatrix/rMatrix_a1s.dat',names=colNames)

dict_channels = {'p1':p1,'p2':p2,'a1':a1}


# If plot is 0 then no plotting, any other value it will generate the plots
# 'n' for numerical solution plots, 'a' for analytic solution plots
plot_n = 0
plot_a = 1


# """
# Perform the analysis over all the channels
channels = ['p1','p2']#,'a1']
for ch in channels:

    print(ch,' channel:\n-Working on Legendre fitting')
    angle = np.array([0,15,30,45,60,75,90])

    chan = dict_channels[ch]

    energy_chan = chan['Energy'].values
    angle_chan = chan['Angle'].values
    cross_chan = chan['Cross-section'].values
    crossErr_chan = chan['Error'].values


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
        masked_cross_chan = chan_temp['Cross-section'].values
        masked_err_chan = chan_temp['Error'].values

        # # Debugging diving by 0 errors ( only in p1 file, 90deg detector, lowest energy 3.99MeV, fixed manually)
        # for blah in range(len(masked_err_chan)):
        #     print(masked_nrg_chan[blah],'\t',masked_ang_chan[blah],'\t',masked_err_chan[blah]) #### TESTING IN/outputs

        temp_weights = [1/_**2 for _ in masked_err_chan]
        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED


        # Up to Legendre Polynomial "BLAH" in legendre_order[list] index starts from 0, can go up to 5
        leg_ord = 2


        # Loop through the number of terms used in Legendre expansion:
        # -- a0, a0+a2, ... , a0+a2+a4+a6+a8+a10
        for ind in range(0,leg_ord+1):
            coef,full = np.polynomial.legendre.legfit(np.cos(masked_ang_chan),masked_cross_chan,deg=legendre_order[ind],full=True,w=temp_weights)

            # print('Results from cookbook python minimization:\n',coef)
            # chi2_det(masked_cross_chan,masked_err_chan,ind)
            results = chi2_mat(masked_cross_chan,masked_err_chan,ind)

            # TURN ON/OFF if you want to debug values
            # stop = input('...Pause...\n\n')


            dict_ord[str(legendre_order[ind][-1])].append(coef)
            dict_ord_a[str(legendre_order[ind][-1])].append(results[0])
            dict_ord_err_a[str(legendre_order[ind][-1])].append(results[1])
            dict_ord_x2_a[str(legendre_order[ind][-1])].append(results[2])
            dict_ord_x2ndf_a[str(legendre_order[ind][-1])].append(results[3])
            dict_ord_pVal_a[str(legendre_order[ind][-1])].append(results[4])

            if plot_n:
                plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])

            elif plot_a:
                # Note: The legendre.legval function requires the coefficients
                # of all terms up to whatever order you want, must pad the even
                # terms with 0's
                temp_coef = convert(results[0])
                temp_coef_err = convert(results[1])

                # # if nrg == 4.30930:
                # if nrg < 4.962 and nrg > 4.96:
                #     print(ch,legendre_order[ind],nrg)
                #     # df = pd.DataFrame(data=np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),index=ang_knots,columns=['legVal'])
                #     df = pd.DataFrame(data=ang_knots,columns=['angle'])
                #     df = df.assign(legVal=pd.Series(np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),index=df.index).values)
                #     df.to_excel('legendre_out/DATA/analytically/%s/a%d/%s_a%d_legendrePointsAtSingleEnergy.xlsx'%(ch,legendre_order[ind][-1],ch,legendre_order[ind][-1]))

                # Scatter has bug when the y-values are small, the y-axis does not autoscale
                # plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])
                if ind ==2:
                    plt.plot(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),temp_coef),c=color_order[ind],label='%d$^{th}$ Order'%legendre_order[ind][-1])#: $\chi^{2}=%5.2f$'%(legendre_order[ind][-1],results[2]))

        # Angular distribution plots
        if plot_n or plot_a:
            plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
            plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
            plt.xlabel('Angle',fontsize=14)
            plt.ylabel('Differential Cross-Section (barns/sr)', fontsize=14)
            plt.title('$^{24}$Mg($\\alpha$,p$_{%s}\\gamma$) - %5.4f MeV' % (ch[-1],nrg),fontsize=20)
            plt.legend(prop={'size': 14})
            plt.savefig('legendre_out/fits/%s/%5.4fMeVFit.png'%(ch,nrg),dpi=300)
            plt.clf()


    # Convert lists into arrays
    for ind in range(0,leg_ord+1):
        dict_ord_a[str(legendre_order[ind][-1])] = np.array(dict_ord_a[str(legendre_order[ind][-1])])
        dict_ord_err_a[str(legendre_order[ind][-1])] = np.array(dict_ord_err_a[str(legendre_order[ind][-1])])


    # # Now we have all the coefficients for each order fit, lets plot them
    # # plot 'a' coefficients for each order as function of energy as well as
    # # the chi2ndf
    # ord_step = 0
    # print('-Working on plotting the coefficients')
    # for ind in range(0,leg_ord+1):
    #     for ord_step in range(0,leg_ord+1):
    #         plt.clf()
    #         if ind > ord_step: continue
    #         # a = [ (_[ind],) for _ in dict_ord[str(legendre_order[ord_step][-1])]]     # minimization solution
    #         a = [ (_[ind],) for _ in dict_ord_a[str(legendre_order[ord_step][-1])]]     # analytic solution
    #         plt.scatter(energyList,a,s=8)
    #         plt.title('%s channel - Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1],legendre_order[ord_step][-1]))
    #         # plt.savefig('legendre_out/coef_curve/numerically/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=100)
    #         plt.savefig('legendre_out/coef_curve/analytically/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=100)
    #         ord_step += 1
    #
    #     plt.clf()
    #
    #     x2 = dict_ord_x2_a[str(legendre_order[ind][-1])]
    #     x2ndf = dict_ord_x2ndf_a[str(legendre_order[ind][-1])]
    #     # print(str(legendre_order[ind][-1]))
    #     plt.scatter(energyList,x2ndf,s=8)
    #     plt.title('%s channel - $\chi^{2}$ for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1]))
    #     plt.savefig('legendre_out/chi2/%s/a%d/chi2.png'%(ch,legendre_order[ind][-1]),dpi=100)


    plt.clf()


    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    # colLabels = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','chi2','p-value','chi2ndf']
    colLabels = ['a0','a2','a4','a6','a8','a10']
    colLabelss = ['a0_err','a2_err','a4_err','a6_err','a8_err','a10_err']
    for _ in range(0,leg_ord+1):
        # stuff = dict_ord[str(legendre_order[_][-1])]
        # df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_*2+1])
        # df.to_csv('legendre_out/DATA/numerically/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        # df.to_excel('legendre_out/DATA/numerically/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))

        stuff = dict_ord_a[str(legendre_order[_][-1])]
        stuff2 = dict_ord_err_a[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_+1])
        df1 =  pd.DataFrame(data=stuff2,index=energyList,columns=colLabelss[0:_+1])
        df = pd.concat([df, df1], axis=1, sort=False)
        df = df.assign(Chi2=pd.Series(dict_ord_x2_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Pval=pd.Series(dict_ord_pVal_a[str(legendre_order[_][-1])],index=df.index).values)
        df = df.assign(Chi2NDF=pd.Series(dict_ord_x2ndf_a[str(legendre_order[_][-1])],index=df.index).values)
        df.to_csv('legendre_out/DATA/analytically/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        df.to_excel('legendre_out/DATA/analytically/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))


    # Now I have all the 'a_i' coefficients for all the many order fits, need
    # to select/filter the necessary 'a0' terms from the fits. This will be
    # done through the p-value.

    """
    Read in all my a0s and p-values and compare each energies fit side by side
    starting from lowest order to highest, select order which pvalue is maximal,
    append the a0 term + continue and repeat for the next energy step
    """
    a0_final = []
    a0_err_final = []
    p_final = []
    x2_final = []
    order = []
    PVAL = 0.08

    for ind in range(len(ord_00)):
        tester_p = [ord_00_pVal[ind],ord_02_pVal[ind],ord_04_pVal[ind]]
        tester_x2 = [ord_00_x2ndf[ind],ord_02_x2ndf[ind],ord_04_x2ndf[ind]]
        index = tester_p.index(max(tester_p))
        if tester_p[index] == 0:
            print('Testing with chi2')
            index = tester_x2.index(min(tester_x2))

        """
        ## 9/20/19 update, only get a0 term from the a4 fit to keep it all consistent
        """
        # a0_final.append(dict_ord_a[str(legendre_order[index][-1])][ind][0])
        # a0_err_final.append(dict_ord_err_a[str(legendre_order[index][-1])][ind][0])
        # p_final.append(dict_ord_pVal_a[str(legendre_order[index][-1])][ind])
        # order.append(str(legendre_order[index][-1]))
        a0_final.append(dict_ord_a[str(legendre_order[2][-1])][ind][0])
        a0_err_final.append(dict_ord_err_a[str(legendre_order[2][-1])][ind][0])
        p_final.append(dict_ord_pVal_a[str(legendre_order[2][-1])][ind])
        order.append(str(legendre_order[2][-1]))



    # # Save the points
    # df = pd.DataFrame(data=a0_final,index=energyList)
    # df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).values)
    # df = df.assign(Pval=pd.Series(p_final,index=df.index).values)
    # df = df.assign(Order=pd.Series(order,index=df.index).values)
    # df.to_excel('legendre_out/DATA/analytically/%s/%sFits.xlsx'%(ch,ch))


    # Now spline the 'a0' coefficients and plot as function of energy and overlay
    # the original 'a0' coefficents and make sure the splining is done well.

    solidAngle = 4*np.pi
    # temp_ord_00 = solidAngle*np.array([_[0] for _ in ord_00])
    a0_final = solidAngle*np.array(a0_final)
    a0_err_final = solidAngle*np.array(a0_err_final)


    # I want to spline from # <= 4.03514 and spline # > 4.26996
    e_knots_mask1 = energyList <= 4.03514
    spline1_mask = energyList <= 4.03514
    spline1 = CubicSpline(energyList[spline1_mask],a0_final[spline1_mask])
    y_val1 = spline1(energyList[e_knots_mask1])


    e_knots_mask2 = energyList >= 4.26996
    spline2_mask = energyList >= 4.26996
    spline2 = CubicSpline(energyList[spline2_mask],a0_final[spline2_mask])
    y_val2 = spline2(energyList[e_knots_mask2])


    # Plot the splines and the points
    plt.clf()
    plt.plot(energyList[e_knots_mask1],y_val1,c='b')
    plt.plot(energyList[e_knots_mask2],y_val2,c='b')
    plt.scatter(energyList,a0_final)
    plt.yscale('log')
    plt.title('%s channel - a0 vs E - Cross-section' % ch)
    plt.ylabel('Cross-Section (barns)')
    plt.xlabel('E$_\\alpha$ (MeV)')
    plt.savefig('legendre_out/excitationCurve/analytically/%s/%s_a0Curve.png' % (ch,ch),dpi=100)
    plt.clf()
    # plt.show()


    # Save the points, first merge the two splines
    yvals = np.ma.concatenate([y_val1,y_val2])
    xvals = np.ma.concatenate([energyList[e_knots_mask1],energyList[e_knots_mask2]])
    df = pd.DataFrame(data=a0_final,index=energyList,columns=['a0'])
    df = df.assign(a0err=pd.Series(a0_err_final,index=df.index).values)
    # df = df.assign(Pval=pd.Series(p_final,index=df.index).values)
    df = df.assign(fitOrder=pd.Series(order,index=df.index).values)
    df = df.assign(splineX=pd.Series(xvals,index=df.index).values)
    df = df.assign(splineY=pd.Series(yvals,index=df.index).values)
    df.to_excel('legendre_out/DATA/analytically/%s/%sFitPoints.xlsx'%(ch,ch))


    #### Note: Integrating both splines and adding
    # totalCross0 = spline.integrate(min(e_knots[e_knots_mask1]),max(e_knots[e_knots_mask1]))
    # totalCross0 += spline.integrate(min(e_knots[e_knots_mask2]),max(e_knots[e_knots_mask2]))
    # totalCross1 = integrate.quad(spline1,min(e_knots[e_knots_mask1]),max(e_knots[e_knots_mask1]))
    # totalCross1 += integrate.quad(spline2,min(e_knots[e_knots_mask2]),max(e_knots[e_knots_mask2]))
    totalCross2 = integrate.romberg(spline1,min(energyList[e_knots_mask1]),max(energyList[e_knots_mask1]),divmax=20)
    totalCross2 += integrate.romberg(spline2,min(energyList[e_knots_mask2]),max(energyList[e_knots_mask2]),divmax=20)

    # totalCross0 = totalCross0
    # totalCross1 = totalCross1[2]
    # totalCross2 = totalCross2

    print("The absolute cross section is:")
    # print(totalCross0)
    # print(totalCross1)
    print(totalCross2,"Barns")



    # Calculate the reaction rates
    fname = 'rate_Temps.dat'
    temperature = np.loadtxt( fname )

    mu = _m24Mg*_m4He/(_m27Al+_m1H)

    if ch == 'a1':
        mu = _m24Mg*_m4He/(_m24Mg+_m4He) # reduced mass


    rxnRate = []
    for T in temperature:

        E1 = energyList[spline1_mask]
        integrand1 = CubicSpline(E1,a0_final[spline1_mask]*np.exp(-11.604*E1/T)*E1)

        E2 = energyList[spline2_mask]
        integrand2 = CubicSpline(E2,a0_final[spline2_mask]*np.exp(-11.604*E2/T)*E2)

        intg = integrate.romberg(integrand1,min(energyList[e_knots_mask1]),max(energyList[e_knots_mask1]),divmax=20)
        intg += integrate.romberg(integrand2,min(energyList[e_knots_mask2]),max(energyList[e_knots_mask2]),divmax=20)

        rate =  3.7313E10* mu**(-.5) * T**(-1.5) * intg
        rxnRate.append(rate)

    rxnRate = np.array(rxnRate)
    plt.scatter(temperature,rxnRate)
    plt.plot(temperature,rxnRate)
    plt.yscale('log')
    plt.xlim(0,10)
    plt.ylim(1e-30,1e5)
    plt.ylabel('Reaction Rate (cm$^3$ mol$^{-1}$ s$^{-1}$)')
    plt.xlabel('Temperature (T9)')
    plt.savefig('%s_ReactionRate.png'%ch)
    plt.clf()
    df = pd.DataFrame(data=temperature,index=temperature,columns=['T9'])
    df = df.assign(Rate=pd.Series(rxnRate,index=df.index).values)
    df.to_excel('legendre_out/DATA/analytically/%s/a%d/%s_rates.xlsx'%(ch,legendre_order[0][-1],ch))




    # continue
    print('\n\n')
# """

print('DONE!')
