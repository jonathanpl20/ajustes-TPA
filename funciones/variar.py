import numpy as np
import pandas as pd

from funciones.ajuste import ajuste
from funciones.focalcurr import focalcurr
from funciones.fscanTH2 import fscanTH2
from funciones.ruido import ruido
import os

def variar(par, TH, error_rel, cantidad_ruido, N=10):
    # Parametro a modificar ingresado como cadena: {L_c, Pavg_c,  Tp_c, wl_c , D_c, ds_c, Cf_c, alpha_c, R_c}
    # TH.
    # Error relativo a usar.
    # cantidad_ruido.
    # número de puntos

    ####################### PARAMETROS EXPERIMENTALES.
    df = pd.read_csv('Parametros_fscan.txt', delimiter="\t", index_col=0)

    corriente = np.arange(0, 202, 2) #current [mA]

    ######################## PARAMETROS QUE SE ASUMEN COMO LOS REALES
    focal = focalcurr('opt2', corriente) #focal distance

    N = 11  # number of terms in the sum of the theoretical formula
    # NN = int(df.loc['NN']['Valor']) # sample size
    beta = 500e-11  # TPA seed in (m/W )
    L = 0.70e-3  # thickness in (m).
    Pavg = 95e-3  # Average power in (W)
    Tp = 2e-9  # pulse width FWHM in (s)
    wl = 1063.5e-9  # central wavelength in (m)
    D = 1.5e-3  # beam diameter in (m)
    ds = 130e-3  # Distance EFTL-sample (m)
    Cf = 1.1  # beam correction coeff
    alfa = 975  # linear absorption in (1/m) (CdS@790 -2.64e-11, ZnSe@790 4.7720, CdSe@790 369.8)alfa_e = float(df.loc['alfa_e']['Valor'])*float(df.loc['alfa_e']['Valor2'])
    R = 0.33  # reflection percentage (CdS@790 = 0.15670, ZnSe@790 0.18164)
    freq = 11e3  # Frecuencia en HZ
    n = 2

    path = os.getcwd() + '\\Datos\data_simulada.txt'  # Ruta de archivo donde se guardaran los datos simulados

    ## Crea señal con ruido
    T = fscanTH2(focal, beta, L, Pavg, Tp, wl, D, ds, Cf, alfa, R, N)
    data = np.column_stack([corriente, ruido(T, cantidad_ruido)])  # corriente, señal simulada con ruido
    np.savetxt(path, data, fmt=['%d', '%.8f'])

    df = pd.read_csv('Parametros_fscan.txt', delimiter="\t", index_col=0)
    parametro = float(df.loc[par, 'Valor'])
    rango = np.linspace(parametro - error_rel * parametro, parametro + error_rel * parametro, N)
    betas2 = np.zeros(N)
    errorb2 = np.zeros(N)


    for j in range(len(rango)):

        df.loc[par, 'Valor'] = str(rango[j])
        df.loc['TH', 'Valor'] = str(TH)
        df.loc[par.replace('c', 'e'), 'Valor'] = str(0)

        yy = ajuste()

        if len(yy) == 2:
            betas2[j] = yy[0]
            errorb2[j] = yy[1]

        else:
            print('no')


    return betas2, errorb2, rango, parametro
