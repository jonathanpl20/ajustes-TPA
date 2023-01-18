import numpy as np
from matplotlib import pyplot as plt
from funciones.variar import variar
from scipy import stats


def funcion(parametro, th, er_rel, cantidad_ruido, N):
    betas2, errorb2, rango, p = variar(parametro, th, er_rel, cantidad_ruido, N)
    beta = 500e-11

    # Ajuste \beta vs parametro
    slope, intercept, r_value, p_value, std_err = stats.linregress(rango, betas2)
    yy = slope * rango + intercept
    # Agregando barras de error al error relativo.
    x1 = abs((rango[0:int(N / 2)] - p) / p)

    erb1 = abs((betas2[0:int(N / 2)] - beta * 1e2)) / (beta * 1e2)  # error relativo del beta.

    x2 = abs((rango[int(N / 2):int(N) + 1] - p) / p)

    erb2 = abs((betas2[int(N / 2):int(N) + 1] - beta * 1e2)) / (beta * 1e2)  # error relativo del beta.

    ee1l = abs((betas2[0:int(N / 2)] - errorb2[0:int(N / 2)] - beta * 1e2)) / (
            beta * 1e2)  # Error asociado al error relativo del beta. 1

    ee1u = abs((betas2[0:int(N / 2)] + errorb2[0:int(N / 2)] - beta * 1e2)) / (
            beta * 1e2)  # Error asociado al error relativo del beta. 1

    eee1 = np.array([abs(ee1l - erb1), abs(ee1u + erb1)])

    ee2l = abs((betas2[int(N / 2):int(N) + 1] + errorb2[int(N / 2):int(N) + 1] - beta * 1e2)) / (
            beta * 1e2)  # Error asociado al error relativo del beta

    ee2u = abs((betas2[int(N / 2):int(N) + 1] + errorb2[int(N / 2):int(N) + 1] - beta * 1e2)) / (
            beta * 1e2)  # Error asociado al error relativo del beta 2.

    eee2 = np.array([abs(ee2l - erb2), abs(ee2u + erb2)])

    nombres = {'D_c': 'Diametro', 'L_c': 'Longitud', 'Pavg_c': 'Potencia', 'Tp_c': 'Ancho Temporal',
               'wl_c': 'Longitud de onda', 'ds_c': 'Distancia muestra-lente', 'Cf_c': 'Coeficiente',
               'alfa_c': 'Coeficiente de absorción', 'R_c': 'Coeficiente de reflexión'}

    font = {'family': 'Arial',
            'weight': 'bold',
            'size': 40}

    # print('erb1', eee1)
    # plt.rcParams['figure.fig-size'] = (15, 40)
    f = plt.figure(figsize=(12, 5), dpi=200)

    # print('p', p)
    # print('rango', rango)
    print('betas_f', betas2)
    # print('error valor real', errorb2)

    f = plt.figure(figsize=(12, 5), dpi=200)

    plt.subplot(1, 2, 1)
    plt.plot(rango, np.ones(len(rango)) * beta * 1e2, 'k--')
    plt.plot(rango, betas2, color='k', marker='o')
    plt.plot(rango, yy, color='g', marker='o', label='Ajuste: $R^{2}$=' + str(round(r_value ** 2, 4)))
    plt.errorbar(rango, betas2, yerr=errorb2, color='r', label='Valor real')
    plt.ylabel(r"$\beta$", size=10)
    plt.legend()
    plt.xlabel(nombres[parametro])

    plt.subplot(1, 2, 2)
    plt.plot(x1, erb1, 'ro', label='(' + nombres[parametro] + '$-e_{r}$' + nombres[parametro] + ')')
    plt.errorbar(x1, erb1, yerr=eee1, color='r', label='Valor real')
    plt.plot(x2, erb2, 'ko', label='(' + nombres[parametro] + '$+e_{r}$' + nombres[parametro] + ')')
    plt.errorbar(x2, erb2, yerr=eee2, color='k', label='Valor real')
    plt.ylabel("Error relativo en " r"$\beta$")
    plt.xlabel("Error relativo en el " + nombres[parametro])
    plt.legend()

    return betas2
