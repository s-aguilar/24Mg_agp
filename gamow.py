import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from AtomicMassTable import GetElement


"""
INCOMPLETE
"""


_e = 1.60218E-19        # eV
# _e = 1.60218E-25        # MeV
_k = 8.617E-2           # MeV / GK
_hbarc = 197.327        # MeV fm

_K = 8.85418782E-12     # Permitivity of free space

_fineStruc = 1/137.036

# _hbarc = 197.327E-13    # MeV cm
_c = 2.99792458E10      # cm/s
_u = 931.494/_c**2      # MeV / c^2


def gamowEnergy(Z1,Z2,mu,T):

    E0 = (mu*_u/2)**(1/3) * (np.pi*(Z1*Z2*(1/_e**2))*(_k*T) / _hbarc) ** (2/3)
    E0_test = 0.122 * ((Z1*Z1)*(Z2*Z2)*mu)**(1/3) * T**(2/3)
    # print('fst:',_fineStruc,_e**2/_hbarc)
    return E0_test

def gamowWindow(Z1,Z2,mu,T):

    delta_E0 = 4*(gamowEnergy(Z1,Z2,mu,T)*_k*T/3)**(1/2)
    delta_E0_test = 0.2368*((Z1*Z1)*(Z2*Z2)*mu)**(1/6) * T**(5/6)

    print('compare windows:')
    print(delta_E0)
    print(delta_E0_test)
    return delta_E0_test


if __name__ == '__main__':

    # GetElement returns (index,Z,A,Mass,'element')
    _m1H = GetElement(1,1)[3]
    _m4He = GetElement(2,4)[3]
    _m24Mg = GetElement(12,24)[3]
    _m27Al = GetElement(13,27)[3]

    mu = _m24Mg*_m4He/(_m27Al+_m1H)
    z1 = GetElement(12,24)[1]       # 24Mg
    z2 = GetElement(2,4)[1]         # 4He
    T = 2       # GK


    print('\nGamow peak at:\n\t',gamowEnergy(z1,z2,mu,T),'MeV')
    print('\nGamow window of:\n\t',gamowWindow(z1,z2,mu,T),'MeV')
