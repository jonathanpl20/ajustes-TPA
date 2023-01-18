# Experimental data normalization. Option 1
import numpy as np


def normalizacion(datos, opcion):
    if opcion == 'alas':
        datos = datos / np.mean(np.array([datos[0:5], datos[-6:-1]]))
    elif opcion == 'izq':
        datos = datos / np.mean(datos[0:5])
    elif opcion == 'der':
        datos = datos / np.mean(datos[-6:-1])
    return datos
