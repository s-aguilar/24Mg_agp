import os
import numpy as np
# import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

"""
Reads in the cross-section data produced by 'crossSectionExp.py' file and
performs a fit using even terms in the legendre polynomial on the angular
distributions.
# """


# Read in the data into dataframe
colNames = ['Energy','Angle','Cross-section','Error']
p1 = pd.read_table('rMatrix/rMatrix_p2s.dat',names=colNames)
p2 = pd.read_table('rMatrix/rMatrix_p2s.dat',names=colNames)
a1 = pd.read_table('rMatrix/rMatrix_a1s.dat',names=colNames)

angle = np.array([0,15,30,45,60,75,90])

# for ang in angle:
#     is_angle = p1['Angle']==ang     # create boolean mask
#     p1_temp = p1[is_angle]

energy_p1 = p1['Energy'].values
angle_p1 = p1['Angle'].values
cross_p1 = p1['Cross-section'].values
crossErr_p1 = p1['Error'].values

p1_0_ord_coef = []
p1_2_ord_coef = []
p1_4_ord_coef = []
p1_6_ord_coef = []

for nrg in energy_p1:
    is_nrg = p1['Energy']==nrg      # create boolean mask
    p1_temp = p1[is_nrg]

    masked_nrg_p1 = p1_temp['Energy'].values
    masked_ang_p1 = np.radians(p1_temp['Angle'].values)
    masked_cross_p1 = p1_temp['Cross-section'].values
    masked_err_p1 = p1_temp['Error'].values

    # Angular distribution plots
    plt.errorbar(masked_ang_p1,masked_cross_p1,yerr=masked_err_p1,c='k',fmt='o')
    # plt.ylim(-0.00001,0.00001)
    # spline = CubicSpline(masked_ang_p1,masked_cross_p1)
    ang_knots = np.radians(np.linspace(0,90,91)) # RADIANSIFIED
    # y_val = spline(ang_knots)
    # plt.scatter(ang_knots,y_val,c='r')
    # plt.xticks(ticks=[0,np.pi()/4,np.pi()/2],labels=['0','$\pi/4$','$\pi/2$'])
    plt.title('%5.4f MeV - Angular Distribution'%nrg)
    # plt.show()

    print('====================================================================')
    print("0 order")
    # coef0 = np.polynomial.legendre.Legendre.fit(np.cos(masked_ang_p1),masked_cross_p1,deg=0,domain=[0,90])
    # print (coef0)

    coef0,full0 = np.polynomial.legendre.legfit(np.cos(masked_ang_p1),masked_cross_p1,deg=0,full=True)
    plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef0),c='b',s=8,label='0')
    # print(full0)
    print (coef0)

    print('\n0 + 2 order')
    # coef2 = np.polynomial.legendre.Legendre.fit(np.cos(masked_ang_p1),masked_cross_p1,deg=[0,2],domain=[0,90])
    # print(np.polynomial.legendre.legval(ang_knots,coef2))
    # print (coef2)

    coef2,full2 = np.polynomial.legendre.legfit(np.cos(masked_ang_p1),masked_cross_p1,deg=[0,2],full=True)
    plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef2),c='y',s=8,label='0+2')
    print (coef2)
    #
    print('\n0 + 2 + 4 order')
    # coef4,full4 = np.polynomial.legendre.legendre.legfit(masked_ang_p1,masked_cross_p1,deg=[0,2,4],domain=[0,90])
    # plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef4),c='y',s=8,label='0+2+4')
    # print (coef4)

    coef4,full4 = np.polynomial.legendre.legfit(np.cos(masked_ang_p1),masked_cross_p1,deg=[0,2,4],full=True)
    plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef4),c='m',s=8,label='0+2+4')
    print (coef4)

    #
    print('\n0 + 2 + 4 + 6 order')
    # coef6,full6 = np.polynomial.legendre.legendre.legfit(masked_ang_p1,masked_cross_p1,deg=[0,2,4,6],domain=[0,90])
    # plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef6),c='y',s=8,label='0+2+4+6')
    # # print (coef6)

    coef6,full6 = np.polynomial.legendre.legfit(np.cos(masked_ang_p1),masked_cross_p1,deg=[0,2,4,6],full=True)
    plt.scatter(ang_knots,np.polynomial.legendre.legval(np.cos(ang_knots),coef6),c='g',s=8,label='0+2+4+6')

    print (coef6)
    plt.legend()
    plt.show()
    print('====================================================================')











# energy_p2 = p2['Energy'].values
# angle_p2 = p2['Angle'].values
# cross_p2 = p2['Crosssection'].values
# crossErr_p2 = p2['Error'].values
#
# energy_a1 = a1['Energy'].values
# angle_a1 = a1['Angle'].values
# cross_a1 = a1['Crosssection'].values
# crossErr_a1 = a1['Error'].values


    # # Make the Cross-Section plot
    # plt.clf()
    # plt.plot(_a1Ealpha,_a1Cross)
    # plt.errorbar(_a1Ealpha,_a1Cross,yerr=_a1Cross_err,fmt='b.',markersize='2')
    # plt.yscale('log')
    # plt.ylim(1e-6,1)
    # plt.xlim(4,5.6)
    # plt.xlabel('$E_{\\alpha}$ (MeV)')
    # plt.ylabel('Differential Cross-Section (barns/sr)')
    # plt.title('a1 %s$^{\circ}$'%ang)
    # plt.savefig('crossSection/A1/a1_%s.png'%ang,dpi=300)
    # plt.clf()
    #
    # with open("rMatrix/rMatrix_a1.dat","a") as f:
    #     for loop in range(229):
    #         printOut= '%f \t %d \t %.8f \t %.8f \n' %(_a1Ealpha[loop],_Angle[loop],_a1Cross[loop],_a1Cross_err[loop])
    #         f.write(printOut)
# """

print('DONE!')
