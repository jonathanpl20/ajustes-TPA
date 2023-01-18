# Cambio de corriente a distancia focal
def focalcurr(ecu, corr):
    if ecu == 'opt1':
        focal = 1 / (0.045 * corr + 1.522)
    elif ecu == 'opt2':
        focal = (1 / (5.9e-5 * corr + 4.9e-3) + 1.7) * 1e-3

    return focal
