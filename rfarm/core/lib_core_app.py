"""
Library Features:

Name:           lib_core_app
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
                Mirko D'Andrea (mirko.dandrea@cimafoundation.org)
                Francesco Silvestro (francesco.silvestro@cimafoundation.org)
                Nicola Rebora (nicola.rebora@cimafoundation.org)
Date:           '20230524'
Version:        '3.5.2'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import numpy as np

from copy import deepcopy

import rfarm.core.lib_core_fx as lib_core_fx
from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to run RF model (nwp mode)
def exec_nwp(X, xscale, tscale, cssf, ctsf,
             f=None, celle3_rainfarm=None,
             sx=None, st=None):
    """
    x = rainfarm(X, sx, st, xscale, cssf, csst)
    INPUT
        X = matrice 3D di dimensioni nx * nx * nt, nx
        xscale = fattore di scala spaziale (numero intero)
        tscale = fattore di scala temporale (numero intero)
        cssf = scala spaziale affidabile (numero intero, 1 = risoluzione nativa)
        ctsf = scala temporale affidabile (numero intero, 1 = risoluzione nativa)
        f = ampiezze per il campo metagaussiano,  OPZIONALE
        celle3_rainfarm = matrice per l'interpolazione,   OPZIONALE
        sx = penza spettrale spaziale,    OPZIONALE
        st = penza spettrale temporale,    OPZIONALE

    OUTPUT
        x = matrice 3D con il campo disaggregato di dimensioni (nx*xscale) * (nx*xscale) * (nt*tscale)
        f = ampiezze per il campo metagaussiano calcolate
    """

    # info execution start
    log_stream.info(' ------> Run RF model (NWP Mode) ... ')

    # X and Y scale(s)
    yscale = deepcopy(xscale)

    # exponential of gaussian curve
    alfa = 1
    # X instantaneous rain data
    nx, ny, nt = X.shape

    #Inserito da Francesco 19/09/2023
    #Removes NAN values
    X[np.isnan(X)]=0.0
    '''
    plt.figure(1)
    plt.imshow(np.nansum(X,2), interpolation='none')
    plt.colorbar()
    plt.show()
    '''
    if cssf == 1 and ctsf == 1:
        pa = X
    else:
        pa = lib_core_fx.agg_xyt(X[:, :, 0:nt], nx / cssf, ny / cssf, nt / ctsf)

    '''
    # Debug
    from scipy.io import savemat
    from os.path import join
    path = '/home/fabio/Desktop/PyCharm_Workspace/RainFarm/Data/Data_Dynamic/Source/ecmwf0100/'
    file = 'rain_debug.mat'
    data = {'data': pa}
    savemat(join(path,file), data)
    import matplotlib.pylab as plt
    plt.figure(1); plt.imshow(pa[:, : ,0], interpolation='none'); plt.colorbar(); plt.show()
    '''

    # find spectral slopes
    sy = sx
    if sx is None and st is None:
        fxp, fyp, ftp = lib_core_fx.fft3d(X)

        kmin, kmax = 2, min(15, len(fxp) - 1)
        wmin, wmax = 2, min(9, len(ftp) - 1)

        sx, sy, st = lib_core_fx.fitallslopes(
            fxp, fyp, ftp,
            np.arange(kmin, kmax),
            np.arange(wmin, wmax),log_stream)

        # INIT: prepare f field for metagauss
        # np.random.rand('state',sum(100*clock))  ##ok<RAND> #
        # seme random differente per ogni run del modello

    # info spectral slopes
    log_stream.info(' -------> Slopes: sx=%f sy=%f st=%f' % (sx, sy, st))

    # initialize meta-gaussian field
    if f is None:
        f = lib_core_fx.initmetagauss(sx, st, nx * xscale, nt * tscale)
    # generate meta-gaussian field
    g = lib_core_fx.metagauss(f)

    # compute non-linear transform
    r = np.exp(alfa * g[0:nx * xscale, 0:ny * yscale, 0:nt * tscale])

    # aggregate field to be the same as pa
    ga = lib_core_fx.agg_xyt(r, nx / cssf, ny / cssf, nt / ctsf)
    ca = pa / ga
    if celle3_rainfarm is None:
        cai = lib_core_fx.interpola_xyt(ca, nx * xscale, ny * yscale, nt * tscale)
    else:
        cai = np.reshape(ca[celle3_rainfarm], (nx * xscale, ny * yscale, nt * tscale),
                         order='F')

    # define rainfarm output fields
    x = cai * r

    # info execution end
    log_stream.info(' ------> Run RF model (NWP Mode) ... OK')

    # return variable(s)
    return x, f

# -------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# method to run RF model (expert forecast mode)
def exec_expert_forecast(X, xscale, tscale, cssf, ctsf, pth=0.25,
                         f=None, celle3_rainfarm=None,
                         sx=None, st=None, nx=None, ny=None, nt=None):

    """
    x = rainfarm(X, sx, st, xscale, cssf, csst)
    INPUT
        X = matrice 3D di dimensioni nx * nx * nt, nx
        xscale = fattore di scala spaziale (numero intero)
        tscale = fattore di scala temporale (numero intero)
        cssf = scala spaziale affidabile (numero intero, 1 = risoluzione nativa)
        ctsf = scala temporale affidabile (numero intero, 1 = risoluzione nativa)
        f = ampiezze per il campo metagaussiano,  OPZIONALE
        celle3_rainfarm = matrice per l'interpolazione,   OPZIONALE
        sx = penza spettrale spaziale,    OPZIONALE
        st = penza spettrale temporale,    OPZIONALE

    OUTPUT
        x = matrice 3D con il campo disaggregato di dimensioni (nx*xscale) * (nx*xscale) * (nt*tscale)
        f = ampiezze per il campo metagaussiano calcolate
    """

    # info execution start
    log_stream.info(' ------> Run RF model (ExpertForecast Mode) ... ')

    # Nx and NY
    if X is None:
        X = np.float64(1.0)
        nx = nx
        ny = ny
        nt = nt
        cssf = nx
        ctsf = nt
    else:
        log_stream.error(' ===> ')
        raise NotImplemented('Case not implemented yet')

    # exponential of gaussian curve
    alfa = 1
    if cssf == 1 and ctsf == 1:
        pa = X
    else:
        if X.ndim == 0:
            pa = X
        elif X.ndim == 3:
            pa = lib_core_fx.agg_xyt(X[:, :, 0:nt], nx / cssf, ny / cssf, nt / ctsf)
        else:
            log_stream.error(' ===> X datasets dimensions are not supported')
            raise NotImplemented('Case not implemented yet')

    '''
    # DEBUG
    import scipy.io as sio
    Data = {}
    Data['rain'] = X
    sio.savemat('rain_debug.mat',Data)
    plt.figure(1); plt.imshow(X[:,:,0], interpolation='none'); plt.colorbar(); plt.show()
    '''

    # find spectral slope(s)
    if sx is None and st is None:

        fxp, fyp, ftp = lib_core_fx.fft3d(X)

        kmin, kmax = 3, min(15, len(fxp) - 1)
        wmin, wmax = 3, min(9, len(ftp) - 1)

        sx, sy, st = lib_core_fx.fitallslopes(
            fxp, fyp, ftp,
            np.arange(kmin, kmax + 1), np.arange(wmin, wmax + 1))

        # INIT: prepare f field for metagauss
        # np.random.rand('state',sum(100*clock))  ##ok<RAND> #
        # seme random differente per ogni run del modello
    else:
        sx, sy, st = sx, sx, st
    # info spectral slope(s)
    log_stream.info(' -------> Slopes: sx=%f sy=%f st=%f' % (sx, sy, st))

    # initialize meta-gaussian field
    if f is None:
        f = lib_core_fx.initmetagauss(sx, st, ny * xscale, nt * tscale)

    # generate meta-gaussian field
    g = lib_core_fx.metagauss(f)

    # compute non-linear transform
    g = np.exp(alfa * g[0:ny * xscale, 0:ny * xscale, 0:nt * tscale])

    # aggregate field to be the same as pa
    ga = lib_core_fx.agg_xyt(g, nx / cssf, ny / cssf, nt / ctsf)
    ca = pa / ga

    # define rainfarm output fields
    if ca.ndim == 0:
        x = ca * g
        x[np.where(x <= pth)] = 0.0
    elif ca.ndim == 3:
        if celle3_rainfarm is None:
            cai = lib_core_fx.interpola_xyt(ca, nx * xscale, ny * xscale, nt * tscale)
        else:
            cai = np.reshape(ca[celle3_rainfarm], (nx * xscale, ny * xscale, nt * tscale),
                             order='F')
        x = cai * g
    else:
        log_stream.error(' ===> CA datasets dimensions are not supported')
        raise NotImplemented('Case not implemented yet')

    # info execution end
    log_stream.info(' ------> Run RF model (ExpertForecast Mode) ... OK')

    # return variable(s)
    return x, f
    # --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------
