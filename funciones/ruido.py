import numpy as np


def ruido(señal, desv):  # Función para añadir ruido Gaussiano dada una desviación estandar desv
    # Recibe la señal y el porcentaje de ruido
    # Añade ruido generado por distribución gausiana hasta n desviaciones estandares.
    señal_ruido = np.zeros(len(señal))
    señal_ruido = [np.random.normal(señal[i], desv) for i in range(0, len(señal))]
    return señal_ruido
