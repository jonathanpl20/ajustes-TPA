import numpy as np


def fscanTH2(focal,beta,L,Pavg,Tp,wl,D,ds,Cf,alfa,R,N):
    freq = 11e3
    Leff = (1. - np.exp(-alfa*L))/alfa # Effective thickness
    w0 = 2.0*wl*focal*Cf/(np.pi*D) # beam waist
    z0 = np.pi*w0**2/wl # Rayleigh range
    w = w0*np.sqrt(1 + ((ds - focal)/z0)**2) # Beam waist
    I0 = 2*np.log(1 + np.sqrt(2))*Pavg/(Tp*freq*np.pi*w**2) # Peak intensity at sample
    B = beta*(1.0 - R)*I0*Leff

    # Transmitance
    T = []
    for cont in range(0,B.size):
        val = 0.0
        for m in range(0,N + 1):
            ter = 1.0
            for n in range(0,m + 1):
                if n == m:
                    ter = ter*1.0
                else:
                    ter = ter*(2.0*(m - n))/(2.0*(m - n) + 1.0)

            val = val + (-B[cont])**m/(m + 1.0)*ter

        T = np.append(T,val)
    return T
