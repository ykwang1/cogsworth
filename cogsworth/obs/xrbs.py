import numpy as np
import pandas as pd
import astropy.units as u
import astropy.constants as const

from cogsworth.classify import get_x_ray_lum

__all__ = ["get_xray_luminosity"]

def get_xray_luminosity(population=None, bcm=None):
    """
    Calculate the X-ray luminosity for RLOF/wind-fed XRBs and Be-XRBs.

    Parameters
    ----------
    population : :class:`~cogsworth.pop.Population`
        The population for which to compute X-ray luminosity (either supply this 
        after calculating a bcm or a bcm)
    bcm : :class:`~pandas.DataFrame`
        User-specified timestep table - must include these columns: [porb] and
        for each star it must have the columns: [mass, rad, kstar, RRLO, deltam]

    Returns
    -------
    xray_lums : :class:`~pandas.Series`
        Luminosity due to accretion in ergs/s/cm**2
    """

    c = const.c 
    G = const.G
    M_sun = const.M_sun
    R_sun = const.R_sun

    if population is not None:
        bcm = population.bcm

    # find which star is the one accreting mass
    CO_flag = {}
    CO_flag["1"] = bcm["kstar_1"].isin([13, 14]) & (bcm["RRLO_2"] >= 1)
    CO_flag["2"] = bcm["kstar_2"].isin([13, 14]) & (bcm["RRLO_1"] >= 1)

    # initialize all the variables for calculating X-ray lum
    deltam = np.zeros(len(bcm)) 
    mass_acc = np.zeros(len(bcm))
    rad_acc = np.zeros(len(bcm))
    porb = np.zeros(len(bcm))
    kstar = np.zeros(len(bcm))
    mass_don = np.zeros(len(bcm))
    RRLO_don = np.zeros(len(bcm))

    # fill in arrays with appropriate compact object and companion values
    for CO_id in "12":
        companion_id = "12".replace(CO_id, "")
        mask = CO_flag[CO_id]

        deltam[mask] = bcm.loc[mask][f"deltam_{CO_id}"].values
        mass_acc[mask] = bcm.loc[mask][f"mass_{CO_id}"].values
        rad_acc[mask] = bcm.loc[mask][f"rad_{CO_id}"].values
        porb[mask] = bcm.loc[mask]["porb"].values
        kstar[mask] = bcm.loc[mask][f"kstar_{CO_id}"].values
        mass_don[mask] = bcm.loc[mask][f"mass_{companion_id}"].values 
        RRLO_don[mask] = bcm.loc[mask][f"RRLO_{companion_id}"].values

    # Make arrays astropy quantities with the correct units
    deltam = deltam * M_sun / u.year
    mass_acc = mass_acc * M_sun
    rad_acc = rad_acc * R_sun
    porb = porb * u.day
    mass_don = mass_don * M_sun
    
    # call function that calculates xray luminosities using Misra+23, ignore BeXRBs
    xray_lums, _ = get_x_ray_lum(mass_acc, rad_acc, deltam, porb, kstar, mass_don, RRLO_don)

    return xray_lums
