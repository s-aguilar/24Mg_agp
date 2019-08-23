import os
import openpyxl # conda install -c anaconda openpyxl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import CubicSpline


"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.

The 'plot' variable enables/disables the plots demonstrating the different
order legendre polynomial fits
"""




# Read in the data into dataframe
colNames = ['Energy','Angle','Cross-section','Error']
p1 = pd.read_table('rMatrix/rMatrix_p1s.dat',names=colNames)
p2 = pd.read_table('rMatrix/rMatrix_p2s.dat',names=colNames)
a1 = pd.read_table('rMatrix/rMatrix_a1s.dat',names=colNames)

dict_channels = {'p1':p1,'p2':p2,'a1':a1}


# """ once everything is debugged put it all in the for loop here
channels = ['p1','p2','a1']
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

    ord_0 = []
    ord_2 = []
    ord_4 = []
    ord_6 = []
    ord_8 = []
    ord_10 = []

    dict_ord = {'0':ord_0,'2':ord_2,'4':ord_4,'6':ord_6,'8':ord_8,'10':ord_10}


    # If plot is 0 then no plotting, any other value it will generate the plots
    plot = 0


    # Passing list to a set does not preserve the original order (ascending) so
    # manually sort it back into a new list
    temp_nrg = set()
    energyList = [x for x in energy_chan if not (x in temp_nrg or temp_nrg.add(x))]


    # perform Legendre fit over the angular distribution at each energy
    for nrg in energyList:
        is_nrg = chan['Energy']==nrg      # create boolean mask
        chan_temp = chan[is_nrg]

        masked_nrg_chan = chan_temp['Energy'].values
        masked_ang_chan = np.radians(chan_temp['Angle'].values)
        # print(nrg, masked_ang_chan)
        masked_cross_chan = chan_temp['Cross-section'].values
        masked_err_chan = chan_temp['Error'].values
        for blah in masked_err_chan:
            print(blah)
        temp_weights = [1/(_*_) for _ in masked_err_chan]


        ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED
        spline = CubicSpline(masked_ang_chan,masked_cross_chan) # not really needed but nice to look at
        y_val = spline(ang_knots)

        for ind in range(0,6):
            coef,full = np.polynomial.legendre.legfit(np.cos(masked_ang_chan),masked_cross_chan,deg=legendre_order[ind],full=True,w=temp_weights)

            # THIS NEEDS TO BE THOUGHT OUT, WHEN TO SCALE BY 4PI
            # Append the coefficients of the fits to the lists and multiply by 4pi
            # coef = coef*(4*np.pi)

            dict_ord[str(legendre_order[ind][-1])].append(coef)

            if plot:
                plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])

        # Angular distribution plots
        if plot:
            plt.errorbar(masked_ang_chan,masked_cross_chan,yerr=masked_err_chan,c='k',fmt='o')    # Original data points
            plt.plot(ang_knots,y_val,c='k')                                                 # Spline of data points
            plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
            plt.title('%s - %5.4f MeV - Angular Distribution' % (ch,nrg))
            plt.legend()
            plt.show()

    """ Can turn this back on after
    # Now we have all the coefficients for each order fit, lets plot them
    ## plot 'a' coefficients for each order as function of energy
    ord_step = 0
    print('-Working on plotting the coefficients')
    for ind in range(0,6):
        # a0 = []
        # for _ in dict_ord[str(legendre_order[ind][-1])]:
        #     a0.append(_[0])
        for ord_step in range(0,6):
            plt.clf()
            if ind > ord_step: continue
            a = [ (_[ind],) for _ in dict_ord[str(legendre_order[ord_step][-1])]]
            plt.scatter(energyList,a,s=8)
            plt.title('%s channel - Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (ch,legendre_order[ind][-1],legendre_order[ord_step][-1]))
            # plt.show()
            plt.savefig('legendre_out/coef_curve/%s/a%d/a%dFit.png'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=900)
            ord_step += 1
    # """

    # Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
    print('-Writing coefficients into .csv and .xlsx files')
    colLabels = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10']
    for _ in range(0,6):
        stuff = dict_ord[str(legendre_order[_][-1])]
        df = pd.DataFrame(data=stuff,index=energyList,columns=colLabels[0:_*2+1])
        df.to_csv('legendre_out/CSV/%s/a%d/a%dFit.csv'%(ch,legendre_order[_][-1],legendre_order[_][-1]))
        df.to_excel('legendre_out/CSV/%s/a%d/a%dFit.xlsx'%(ch,legendre_order[_][-1],legendre_order[_][-1]))


    # Now spline the 'a' coefficients and plot as function of energy and overlay
    # the original 'a' coefficents and make sure the splining is done well.

    #####Currently only doing the a0 of the order 0 fit to debug the process!!!!!!
    e_knots = np.linspace(min(energyList),max(energyList),1001)
    spline = CubicSpline(energyList,ord_0)
    y_val = spline(e_knots)


    # TO DO:need to scale by 4pi somewhere!!!!!!!!

    plt.clf()
    plt.plot(e_knots,y_val)
    a = [ (_[0],) for _ in ord_0]
    plt.scatter(energyList,ord_0)
    plt.title('%s channel - a0 vs E' % ch)
    plt.savefig('legendre_out/excitationCurve/%s/a0Curve.png' % ch,dpi=900)
    # plt.show()

    #### NOtE: NEED TO PROBABLY MAKE 2 SPLINES BECAUSE WE HAVE THAT LARGE step
    # DOWN INTO LOW ENERGY AND THE SPLINE INTERPOLATION IN BETWEEN both
    # regions IS NONSENSE. NEED TO ADDRESS THIS!
    solidAngle = 4*np.pi
    totalCross0 = spline.integrate(min(e_knots),max(e_knots))
    totalCross1 = integrate.quad(spline,min(e_knots),max(e_knots))
    totalCross2 = integrate.romberg(spline,min(e_knots),max(e_knots))

    print("The absolute cross section is:")
    print(totalCross0)
    print(totalCross1)
    print(totalCross2)






    # continue
    print('\n\n')
# """




"""
angle = np.array([0,15,30,45,60,75,90])

# for ang in angle:
#     is_angle = p1['Angle']==ang     # create boolean mask
#     p1_temp = p1[is_angle]

energy_p1 = p1['Energy'].values
angle_p1 = p1['Angle'].values
cross_p1 = p1['Cross-section'].values
crossErr_p1 = p1['Error'].values


legendre_order = [[0],[0,2],[0,2,4],[0,2,4,6],[0,2,4,6,8],[0,2,4,6,8,10]]
color_order = ['b','g','r','c','m','y']

ord_0 = []
ord_2 = []
ord_4 = []
ord_6 = []
ord_8 = []
ord_10 = []

dict_ord = {'0':ord_0,'2':ord_2,'4':ord_4,'6':ord_6,'8':ord_8,'10':ord_10}


# If plot is 0 then no plotting, any other value it will plot
plot = 0


# Passing list to a set does not preserve the original order (ascending) so
# manually sort it back into a new list
temp_nrg = set()
energySet = [x for x in energy_p1 if not (x in temp_nrg or temp_nrg.add(x))]


# perform Legendre fit over the angular distribution at each energy
for nrg in energySet:
    is_nrg = p1['Energy']==nrg      # create boolean mask
    p1_temp = p1[is_nrg]

    masked_nrg_p1 = p1_temp['Energy'].values
    masked_ang_p1 = np.radians(p1_temp['Angle'].values)
    masked_cross_p1 = p1_temp['Cross-section'].values
    masked_err_p1 = p1_temp['Error'].values

    ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED
    spline = CubicSpline(masked_ang_p1,masked_cross_p1) # not really needed but nice to look at
    y_val = spline(ang_knots)

    for ind in range(0,6):
        coef,full = np.polynomial.legendre.legfit(np.cos(masked_ang_p1),masked_cross_p1,deg=legendre_order[ind],full=True)
        dict_ord[str(legendre_order[ind][-1])].append(coef)

        if plot:
            plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef),c=color_order[ind],s=8,label='%d'%legendre_order[ind][-1])

    # Angular distribution plots
    if plot:
        plt.errorbar(masked_ang_p1,masked_cross_p1,yerr=masked_err_p1,c='k',fmt='o')    # Original data points
        plt.plot(ang_knots,y_val,c='k')                                                 # Spline of data points
        plt.xticks(ticks=[0,np.pi/4,np.pi/2],labels=['0','$\pi/4$','$\pi/2$'])
        plt.title('%5.4f MeV - Angular Distribution' % nrg)
        plt.legend()
        plt.show()


# Now we have all the coefficients for each order fit, lets plot them
## plot 'a' coefficients for each order as function of energy
ord_step = 0
for ind in range(0,6):
    # a0 = []
    # for _ in dict_ord[str(legendre_order[ind][-1])]:
    #     a0.append(_[0])
    for ord_step in range(0,6):
        plt.clf()
        if ind > ord_step: continue
        a = [ (_[ind],) for _ in dict_ord[str(legendre_order[ord_step][-1])]]
        plt.scatter(energySet,a,s=8)
        plt.title('Legendre Polynomial a$_{%s}$ coefficients for an up to a$_{%s}$ fit' % (legendre_order[ind][-1],legendre_order[ord_step][-1]))
        # plt.show()
        plt.savefig('legendre_out/coef_curve/%s/a%d/a%dFit'%(ch,legendre_order[ind][-1],legendre_order[ord_step][-1]),dpi=900)
        ord_step += 1

# Now organize into pandas dataframe and write it as both 'csv' and 'xlsx' format
colLabels = ['a0','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10']
for _ in range(0,6):
    stuff = dict_ord[str(legendre_order[_][-1])]
    df = pd.DataFrame(data=stuff,index=energySet,columns=colLabels[0:_*2+1])
    df.to_csv('legendre_out/CSV/a%d/a%dFit.csv'%(legendre_order[_][-1],legendre_order[_][-1]))
    df.to_excel('legendre_out/CSV/a%d/a%dFit.xlsx'%(legendre_order[_][-1],legendre_order[_][-1]))


# Now spline the 'a' coefficients and plot as function of energy and overlay
# the original 'a' coefficents and make sure the splining is done well

# """






"""
#     print("0 order")
#     # coef0 = np.polynomial.legendre.Legendre.fit(np.cos(masked_ang_p1),masked_cross_p1,deg=0,domain=[0,90]) # Investigate the difference between both legendre fitting algos
#     # print (coef0)
"""

print('DONE!')
