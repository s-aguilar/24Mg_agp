import os
import openpyxl # conda install -c anaconda openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from chi2 import chi2_det,chi2_mat
from scipy.interpolate import CubicSpline


"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.

The 'plot_x' variable enables/disables the plots demonstrating the different
order legendre polynomial fits for either 'a' analytic or 'n' numeric fit

Around line 110 there is a stop variable that can be turned off, it is used for
debugging purposes.
"""




# Read in the data into dataframe
colNames = ['Energy','Angle','Cross-section','Error']
p1 = pd.read_table('rMatrix/rMatrix_p1s.dat',names=colNames)
p2 = pd.read_table('rMatrix/rMatrix_p2s.dat',names=colNames)
a1 = pd.read_table('rMatrix/rMatrix_a1s.dat',names=colNames)

dict_channels = {'p1':p1,'p2':p2,'a1':a1}


# """ once everything is debugged put it all in the for loop here
channels = ['p1']#,'p2','a1']
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
    ord_00 = []
    ord_02 = []
    ord_err_00 = []
    ord_err_02 = []
    ord_00_x2 = []
    ord_02_x2 = []

    # Dict for analytic solutions
    dict_ord_a = {'0':ord_00,'2':ord_02}
    dict_ord_err_a = {'0':ord_err_00,'2':ord_err_02}
    dict_ord_x2_a = {'0':ord_00_x2,'2':ord_02_x2}

    # If plot is 0 then no plotting, any other value it will generate the plots
    # 'n' for numerical solution plots, 'a' for analytic solution plots
    plot_n = 0
    plot_a = 0


    # Passing list to a set does not preserve the original order (ascending) so
    # manually sort it back into a new list
    temp_nrg = set()
    energyList = [x for x in energy_chan if not (x in temp_nrg or temp_nrg.add(x))]


    # perform Legendre fit over the angular distribution at each energy
    for nrg in energyList:
        is_nrg = chan['Energy']==nrg      # create boolean mask
        chan_temp = chan[is_nrg]

        print(nrg)

        masked_nrg_chan = chan_temp['Energy'].values
        masked_ang_chan = np.radians(chan_temp['Angle'].values)
        masked_cross_chan = chan_temp['Cross-section'].values
        masked_err_chan = chan_temp['Error'].values

        # # Debugging diving by 0 errors ( only in p1 file, 90deg detector, lowest energy 3.99MeV, fixed manually)
        # for blah in range(len(masked_err_chan)):
        #     print(masked_nrg_chan[blah],'\t',masked_ang_chan[blah],'\t',masked_err_chan[blah]) #### TESTING IN/outputs

        temp_weights = [1/_**2 for _ in masked_err_chan]

        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED
        spline = CubicSpline(masked_ang_chan,masked_cross_chan) # not really needed but nice to look at
        y_val = spline(ang_knots)


        # Up to Legendre Polynomial BLAH in legendre_order[list] index starts from 0
        leg_ord = 1


        # Loop through the number of terms used in Legendre expansion:
        # -- a0, a0+a2, ... , a0+a2+a4+a6+a8+a10
        for ind in range(0,leg_ord+1): # TURN OFF THE HIGHER ORDER LEGENDERE POLYNOMIALS (was 6)
            coef,full = np.polynomial.legendre.legfit(np.cos(masked_ang_chan),masked_cross_chan,deg=legendre_order[ind],full=True,w=temp_weights)

            # print('Results from cookbook python minimization:\n',coef)
            # chi2_det(masked_cross_chan,masked_err_chan,ind)
            results = chi2_mat(masked_cross_chan,masked_err_chan,ind)

            # TURN ON/OFF if you want to debug values
            # stop = input('...Pause...\n\n')


            # THIS NEEDS TO BE THOUGHT OUT, WHEN TO SCALE BY 4PI
            # Append the coefficients of the fits to the lists and multiply by 4pi
            # coef = coef*(4*np.pi)

            dict_ord[str(legendre_order[ind][-1])].append(coef)
            dict_ord_a[str(legendre_order[ind][-1])].append(results[0])
            dict_ord_err_a[str(legendre_order[ind][-1])].append(results[1])
            dict_ord_x2_a[str(legendre_order[ind][-1])].append(results[2][1]) # replace the last index with [0] for chi2 or [1] for chi2ndf

            if plot_n:
                plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])

            elif plot_a:
                pass
        # Angular distribution plots
        if plot_n:
            plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
            plt.plot(ang_knots,y_val,c='k')                                                 # Spline of data points
            plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
            plt.title('%s - %5.4f MeV - Angular Distribution' % (ch,nrg))
            plt.legend()
            plt.show()

        elif plot_a:
            pass


    # Now we have all the coefficients for each order fit, lets plot them
    # plot 'a' coefficients for each order as function of energy as well as
    # the chi2ndf
    ord_step = 0
    print('-Working on plotting the coefficients')
    for ind in range(0,leg_ord+1):
        for ord_step in range(0,leg_ord+1):
            plt.clf()
            if ind > ord_step: continue
            # a = [ (_[ind],) for _ in dict_ord[str(legendre_order[ord_step][-1])]]     # minimization solution
            a = [ (_[ind],) for _ in dict_ord_a[str(legendre_order[ord_step][-1])]]     # analytic solution
            plt.scatter(energyList,a,s=8)
            plt.title('%s channel - Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1],legendre_order[ord_step][-1]))
            # plt.savefig('legendre_out/coef_curve/numerically/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=900)
            plt.savefig('legendre_out/coef_curve/analytically/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=900)
            ord_step += 1

        plt.clf()

        x2ndf = [ (_,) for _ in dict_ord_x2_a[str(legendre_order[ind][-1])] ]   # chi2
        plt.scatter(energyList,x2ndf,s=8)
        plt.title('%s channel - $\chi^{2}$ for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1]))
        plt.savefig('legendre_out/chi2/%s/a%d/chi2.png'%(ch,legendre_order[ind][-1]),dpi=900)


    plt.clf()


    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    colLabels = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10']
    for _ in range(0,leg_ord+1):
        stuff = dict_ord[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_*2+1])
        df.to_csv('legendre_out/CSV/numerically/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        df.to_excel('legendre_out/CSV/numerically/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))

        # WORK on these guys, input the chi2 values
        stuff = dict_ord[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_*2+1])
        df.to_csv('legendre_out/CSV/analytically/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        df.to_excel('legendre_out/CSV/analytically/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))


    # Now spline the 'a' coefficients and plot as function of energy and overlay
    # the original 'a' coefficents and make sure the splining is done well.

    #####Currently only doing the a0 of the order 0 fit to debug the process!!!!!!
    e_knots = np.linspace(min(energyList),max(energyList),1001)
    spline = CubicSpline(energyList,ord_0)
    y_val = spline(e_knots)

    # I want to spline from # <= 4.03514 and spline # > 4.03514
    spline1_mask = energyList <= 4.03514
    spline2_mask = energyList > 4.03514

    spline1 = CubicSpline(energyList[spline1_mask],ord_0[spline1_mask])




    # plt.clf()
    # plt.plot(e_knots,y_val)
    # a = [ (_[0],) for _ in ord_0]
    # plt.scatter(energyList,ord_0)
    # plt.title('%s channel - a0 vs E' % ch)
    # plt.savefig('legendre_out/excitationCurve/numerically/%s/a0Curve.png' % ch,dpi=900)
    # # plt.show()
    #
    # #### NOtE: NEED TO PROBABLY MAKE 2 SPLINES BECAUSE WE HAVE THAT LARGE step
    # # DOWN INTO LOW ENERGY AND THE SPLINE INTERPOLATION IN BETWEEN both
    # # regions IS NONSENSE. NEED TO ADDRESS THIS!
    # solidAngle = 4*np.pi
    # totalCross0 = solidAngle*spline.integrate(min(e_knots),max(e_knots))
    # totalCross1 = solidAngle*integrate.quad(spline,min(e_knots),max(e_knots))
    # totalCross2 = solidAngle*integrate.romberg(spline,min(e_knots),max(e_knots))
    #
    # print("The absolute cross section is:")
    # print(totalCross0)
    # print(totalCross1)
    # print(totalCross2)






    # continue
    print('\n\n')
# """

print('DONE!')
