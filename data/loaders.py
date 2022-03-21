# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import lightkurve as lk
from multiprocessing import Pool
import numpy as np
import pickle
import pandas as pd
import os 
#import resources

def load_lc(fp, fluxtype="PDC", mask=False):

    with open(fp, 'rb') as file:
        lc_list = pickle.load(file)

    fluxes = {"raw": lc_list[7], "corr": lc_list[8], "PDC": lc_list[9]}

    try:
        flux = fluxes[fluxtype]

    except KeyError:
        print("""
        The flux type must be 'raw', 'corr', or 'PDC'. Defaulting to 'PDC'.""")
        flux = fluxes["PDC"]

    finally:
        time = lc_list[6]
        flux_err = lc_list[10]
        quality = lc_list[11]

        if mask:
            mask = lc_list[11] == 0
            flux = flux[mask]
            time = time[mask]
            flux_err = flux_err[mask]
            quality = quality[mask]  # just 0's if masked

        # for meta information
        fluxes.update(
            {"TESS Magnitude": lc_list[3], "filename": fp.split("/")[-1]})
        lc = lk.lightcurve.TessLightCurve(
            time=time, flux=flux, flux_err=flux_err, targetid=lc_list[0],
            quality=quality, camera=lc_list[4], ccd=lc_list[5],
            ra=lc_list[1], dec=lc_list[2], label=f"TIC {lc_list[0]}",
            meta=fluxes
        )

    return lc

# -


