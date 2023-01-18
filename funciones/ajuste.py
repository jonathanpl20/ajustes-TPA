import pandas as pd
import numpy as np
from funciones.normalizacion import normalizacion
from funciones.focalcurr import focalcurr
from scipy.optimize import fmin
from scipy import stats


#Función de Ajuste
def ajuste():

    # Contedra valores de T PA
    TPA=[]
    df = pd.read_csv('Parametros_fscan.txt', delimiter = "\t",index_col=0)

    #ARCHIVO#
    nombreArch = df.loc['Carpeta']['Valor']+'/'+df.loc['Archivo']['Valor']
    nameFit = nombreArch + '-TFscan-FIT' + '.txt'

    f = open(nombreArch + '.txt', 'r')

    data = [] # Creates the list to store data
    f.readline() # Reads the data and store it in one line
    for linea in f:
        linea = linea.strip() # Separates rows
        column = linea.split() # Separates columns
        foco = float(column[0]) # focal distance
        X = float(column[1]) # X Lock-in value
        data = np.append(data,[foco , X] , axis=0) # saves the data in the list
    f.close() # close the file
    data = np.reshape(data,((int(data.size/2),2))) # reshapes the list into an array
    datos = data[:,1] # Experimental data, R lock-in value
    datos = normalizacion(datos,df.loc['Normalizacion']['Valor'])

    Xaxis = df.loc['Xaxis']['Valor'] # if current 'CURR', if focal distance 'FOC'
    Xaxis_u = df.loc['Xaxis']['Valor2']

    if Xaxis == 'FOC':
        focal = data[:,0]*Xaxis_u # Focal distance in (m)
    elif Xaxis == 'CURR':
        corr = data[:,0]*Xaxis_u # corriente en mA
        focal = focalcurr(df.loc['Focal-curr']['Valor'],corr)

    metric = df.loc['metric']['Valor'] # 'DOWN' revisa ubicación punto, dist puntos, y norma
    ECMs = [] # creates metric list
    Betas = [] # creates TPAs list
    Best = [] # creates best TPAs list
    N = int(df.loc['N']['Valor']) #number of terms in the sum of the theoretical formula
    #NN = int(df.loc['NN']['Valor']) # sample size
    NN = 50
    beta0 = float(df.loc['beta0']['Valor']) # TPA seed in (m/W = 1e11 cm/GW)
    ECM2 = float(df.loc['ECM2']['Valor'])
    TH = float(df.loc['TH']['Valor']) # Best fit is TH = 0


    for ii in range(NN):


        global L,alfa,Pavg,freq,Tp,wl,D,Cf,ds,R

        L_c = float(df.loc['L_c']['Valor'])*float(df.loc['L_c']['Valor2']) # thickness in (m).
        L_e = float(df.loc['L_e']['Valor'])*float(df.loc['L_e']['Valor2']) # thickness error

        Pavg_c = float(df.loc['Pavg_c']['Valor'])*float(df.loc['Pavg_c']['Valor2']) # Average power in (W)
        Pavg_e = float(df.loc['Pavg_e']['Valor'])*float(df.loc['Pavg_e']['Valor2']) # Average power error

        Tp_c = float(df.loc['Tp_c']['Valor'])*float(df.loc['Tp_c']['Valor2']) # pulse width FWHM in (s)
        Tp_e = float(df.loc['Tp_e']['Valor'])*float(df.loc['Tp_e']['Valor2']) # pulse width error

        wl_c = float(df.loc['wl_c']['Valor'])*float(df.loc['wl_c']['Valor2']) #central wavelength in (m)
        wl_e = float(df.loc['wl_e']['Valor'])*float(df.loc['wl_e']['Valor2']) # central wavelength error

        D_c = float(df.loc['D_c']['Valor'])*float(df.loc['D_c']['Valor2']) # beam diameter in (m)
        D_e = float(df.loc['D_e']['Valor'])*float(df.loc['D_e']['Valor2']) # beam diameter error

        ds_c = float(df.loc['ds_c']['Valor'])*float(df.loc['ds_c']['Valor2']) # Distance EFTL-sample (m)
        ds_e = float(df.loc['ds_e']['Valor'])*float(df.loc['ds_e']['Valor2']) # distance EFTL-sample error

        Cf_c = float(df.loc['Cf_c']['Valor']) # beam correction coeff
        Cf_e = float(df.loc['Cf_e']['Valor'])

        alfa_c = float(df.loc['alfa_c']['Valor']) # linear absorption in (1/m) (CdS@790 -2.64e-11, ZnSe@790 4.7720, CdSe@790 369.8)
        alfa_e = float(df.loc['alfa_e']['Valor'])

        R_c = float(df.loc['R_c']['Valor']) # reflection percentage (CdS@790 = 0.15670, ZnSe@790 0.18164)
        R_e = float(df.loc['R_e']['Valor'])

        freq = float(df.loc['freq_c']['Valor'])*float(df.loc['freq_c']['Valor2']) # Repetition rate in (Hz)

    # Randomly chosen parameters, normal distribution
        L = np.random.normal(L_c,L_e,1) # # thickness in (m).
        Pavg = np.random.normal(Pavg_c,Pavg_e,1) # Average power in (W)
        Tp = abs(np.random.normal(Tp_c,Tp_e,1)) # Pulse width FWHM in (s). sech
        wl = np.random.normal(wl_c,wl_e,1) #central wavelength in (m)
        D = np.random.normal(D_c,D_e,1) # beam diameter in (m)
        ds = np.random.normal(ds_c,ds_e,1) # EFTL-sample distance (m)
        Cf = np.random.normal(Cf_c,Cf_e,1) # Beam correction factor
        alfa =  abs(np.random.normal(alfa_c,alfa_e,1)) #linear absorption (1/m)
        R = np.random.normal(R_c,R_e,1) # Reflectance

    ########################
    # fscan TPA theoretical function with one variable parameter
        def fscanTH(focal,beta,N):
            Leff = (1. - np.exp(-alfa*L))/alfa # Effective thickness
            w0 = 2.0*wl*focal*Cf/(np.pi*D) # beam waist
            z0 = np.pi*w0**2/wl # Rayleigh range
            w = w0*np.sqrt(1 + ((ds - focal)/z0)**2) # Beam radius
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
    ########################

    ########################

        def metrica(beta):
            ########################global focal,datos
            T = fscanTH(focal,beta,N) # Transmittance
            if metric == 'DOWN':
                mdatos = np.amin(datos) # experimental data maximum value
                pmdatos = np.argmin(datos) # experimental data maximum value position
                mR = np.amin(T) # fitting data maximum value
                pmR = np.argmin(T) # fitting data maximum value position
                dist = np.sqrt((mdatos-mR)**2 + (pmdatos - pmR)**2) + (np.linalg.norm(T-datos))**2 # values distance
            return dist
     ########################


    # Fitting routine
        betaopt = fmin(metrica, beta0, full_output=False, xtol=1e-8, disp=False)
 #   betaopt, ECM, iter, funcalls, warnflag= fmin(metrica, beta0, full_output=True, xtol=1e-8, disp=False)

        Tplot = fscanTH(focal,betaopt,N) # fitted transmitance curve
        ECM = sum((Tplot-datos)**2/datos) # metric of the fitting
        Betas = np.append(Betas,betaopt) # TPAs values vector
        #print('Betas', Betas)
        ECMs = np.append(ECMs,ECM) # Metric values vector
        #print('ECMs', ECMs)

        if ECM < TH:
            miro = np.concatenate((np.array([ECM]),betaopt,L,Pavg,Tp,wl,D,ds,Cf,alfa,R))
            Best = np.append(Best,miro)

        elif ECM < ECM2:
            ECM2 = ECM2


    tamBest = int(len(Best)/11) # size of best-fit values vector
    Best = np.reshape(Best,(tamBest,11)) # Reshape into an array
    #print('Best', Best)
    ########################

    ########################

    #TPA con método 2.

    if len(Best[:,1])!=0:
# Best values average TPA
        Best_avg = np.average(Best[:,1],weights=1/(Best[:,0])**2)*1e11
        alpha = 0.5
        n = Best[:,0].size
        print('number of betas:', Best[:,1].size)
        print('betas:', Best[:,1])
        gdl = n - 1
        confi = 1. - alpha*1e-2
        aux = stats.t.interval(confi,gdl,loc=0,scale=1) # loc sirve para desplazar la distribución, scale para escalarla.
        valor_t = aux[1]
        error = valor_t*np.sqrt(np.average((Best[:,1] - Best_avg*1e-11)**2,weights=1/(Best[:,0])**2))*1e11

        #print('Best_avg', Best_avg)
    # Añade a TPA el TPA promedio producto de la condición T
        #print('Best_avg', Best_avg)
        TPA.append(Best_avg*1e-9)

        #print('error', error*1e-9/np.sqrt(len(Best[:,1] )))
        TPA.append(error*1e-9/np.sqrt(len(Best[:,1] )))

        betas_end = []
        betas_end.append(Best[:,1])
        #print('betas final:', betas_end[-1])
        #print('TPA', TPA)

    return TPA