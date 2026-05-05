""" MOC.py - module for calculating the MOC 
             and derived indeces 
"""

import xgcm
import numpy as np
import xarray as xr

__all__ = [
    "calc_psi_rho",
    "levels_and_colorbarticks",
]


def calc_psi_rho(psi, psiGM, mask=None):
    """Load ty_trans_rho - and sum zonally. Also, if there is a ty_trans_rho_gm variable saved, 
    assume that GM is switched on and load that as well. 
    Most MOM simulations save transport with units of kg/s - convert to Sv."""

    #convert to volume transport
    #rho0 = 1025
    #psi = psi / rho0
    #psiGM = psiGM / rho0
    if mask is not None:
       psi = psi * mask
       psiGM = psiGM * mask
    psi = psi.cf.sum("longitude")
    psiGM = psiGM.cf.sum("longitude")
    """This is equivalent to the Ferret computation
       psi[i=@sum,k=@rsum]-psi[i=@sum,k=@sum]
       but I'm not sure why here I don't need
       to remove the vertical sum"""
    psi = psi.cf.cumsum("vertical")  ##- psi.cf.sum("vertical")
    psi = psi + psiGM

    psi.load()

    return psi


def levels_and_colorbarticks(max_value):
    """ Return the levels and the colorbarticks for the streamfunction plot.
    It may seem complicated but the truth is we just want to avoid the 0 contour
    so that the plot looks soothing to the eye"""

    levels =  np.hstack((np.arange(-max_value, 0, 2), np.flip(-np.arange(-max_value, 0, 2))))
    cbarticks = np.hstack((np.flip(-np.arange(2, max_value+4, 4)), np.arange(2, max_value+4, 4)))

    return levels, cbarticks

