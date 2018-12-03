import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def eff_Plotter(E_gam,eff,err,name):
    # eff = np.array([0.0003740410996,0.0002960081559,0.000251529233])

    labels = ["$10^{-4}$","$10^{-3}$","$10^{-2}$"]

    plt.figure(1)
    plt.yscale('log')
    plt.scatter(E_gam,eff,c='b')
    plt.errorbar(E_gam,eff,yerr=err,fmt='o')
    plt.xlim(400,1800)
    plt.xlabel('Energy (KeV)')
    plt.yticks((.0001,.001,.01),labels)
    plt.ylim(.0001,.001)
    plt.ylabel('Efficiency')
    plt.grid(b=True,which='both',axis='y',alpha=.5)
    plt.tight_layout()

    detName = "EfficiencyPlots/det_"+name
    plt.title(name)
    plt.savefig(detName,format="png",dpi=300)
    plt.show()
    return 0




e_gam = np.array([661.657,1173.228,1332.492])

h0_0eff = np.array([0.0004671717625,0.0003364504989,0.0003064110534])
h0_0err = np.array([0.00001036938748,0.000004127009793,0.000003541523082])
# eff_Plotter(e_gam,h0_0eff,h0_0err,"h0-0")

h0_1eff = np.array([0.0004930287339,0.000371630276,0.0003190202788])
h0_1err = np.array([0.00001068522709,0.000004359216644,0.000003627893029])
# eff_Plotter(e_gam,h0_1eff,h0_1err,"h0-1")

h0_2eff = np.array([0.0002510112315,0.0002001108139,0.0001661487954])
h0_2err = np.array([0.000009873989796,0.000003420824413,0.000002722836366])
# eff_Plotter(e_gam,h0_2eff,h0_2err,"h0-2")

h0_3eff = np.array([0.0003448294122,0.0002784331646,0.0002296181585])
h0_3err = np.array([0.00001019171629,0.00000384736078,0.000003153425562])
# eff_Plotter(e_gam,h0_3eff,h0_3err,"h0-3")

h0_4eff = np.array([0.000311673253,0.0002816223323,0.0002417831389])
h0_4err = np.array([0.000009880259111,0.000003900150691,0.000003185787076])
# eff_Plotter(e_gam,h0_4eff,h0_4err,"h0-4")

h0_5eff = np.array([0.0002757111184,0.0002575047535,0.0002135661294])
h0_5err = np.array([0.00001002530552,0.000003802882797,0.000003054041488])
# eff_Plotter(e_gam,h0_5eff,h0_5err,"h0-5")

h0_6eff = np.array([0.0004413567895,0.0003563727695,0.0003053017777])
h0_6err = np.array([0.0000107546765,0.000004264981687,0.000003560828451])
# eff_Plotter(e_gam,h0_6eff,h0_6err,"h0-6")

h0_7eff = np.array([0.0003012193211,0.0002452232686,0.0002129480826])
h0_7err = np.array([0.00001067013204,0.000003747510688,0.00000303661165])
# eff_Plotter(e_gam,h0_7eff,h0_7err,"h0-7")

h1_0eff = np.array([0.0003556162874,0.0002902669666,0.0002396092794])
h1_0err = np.array([0.00001070470457,0.000003916308529,0.000003187471617])
# eff_Plotter(e_gam,h1_0eff,h1_0err,"h1-0")

h1_1eff = np.array([0.0003605447001,0.0002879116657,0.0002483028074])
h1_1err = np.array([0.00001081755226,0.000003882732312,0.00000322217468])
# eff_Plotter(e_gam,h1_1eff,h1_1err,"h1-1")

h1_2eff = np.array([0.0002824673714,0.0002039279599,0.0001782744316])
h1_2err = np.array([0.0000108734283,0.000003590909533,0.000002832289485])
# eff_Plotter(e_gam,h1_2eff,h1_2err,"h1-2")

h1_3eff = np.array([0.00050444741,0.0003688269102,0.000295221273])
h1_3err = np.array([0.00001169708249,0.000004418729437,0.000003585118686])
# eff_Plotter(e_gam,h1_3eff,h1_3err,"h1-3")

h1_4eff = np.array([0.0004734569038,0.0003698246471,0.0003136748228])
h1_4err = np.array([0.0000110014319,0.000004289199345,0.000003598533893])
# eff_Plotter(e_gam,h1_4eff,h1_4err,"h1-4")


# peak position
h0_0pos = np.array([609.925,1070.58,1210.64])
h0_1pos = np.array([653.887,1158.89,1313.84])
h0_2pos = np.array([665.98,1170.56,1325.44])
h0_3pos = np.array([655.228,1154.35,1307.33])
h0_4pos = np.array([637.456,1124.03,1273.99])
h0_5pos = np.array([668.175,1181.92,1338.95])
h0_6pos = np.array([660.444,1165.27,1320.02])
h0_7pos = np.array([661.339,1172.13,1328.7])
h1_0pos = np.array([651.954,1156.26,1310.45])
h1_1pos = np.array([658.612,1169.4,1325.98])
h1_2pos = np.array([663.31,1178.97,1338.65])
h1_3pos = np.array([648.473,1151.9,1309.16])
h1_4pos = np.array([652.584,1157.09,1312.09])

# interpolate channel to Energy
det_interp = interp1d(e_gam,h1_4pos,fill_value="extrapolate")
energyToChannelp1 = det_interp(843.76)
energyToChannelp2 = det_interp(1014.56)
print("\n\n__________________________________________")
print("p1 at:",energyToChannelp1,"\tp2 at:",energyToChannelp2)
print("__________________________________________")
