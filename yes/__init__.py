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
import pandas as pd
import numpy as np
import pkg_resources

def mount_drive(mount_point):
    try:
        # check if the lookup file for sector one exists
        ref = pd.read_csv(mount_point+"sector1lookup.csv")
        print("Disk already mounted")
    except FileNotFoundError:
        import os
        # If not, try to mount it.
        # presumes the disk is sdc
        os.system(f"sudo mount -o discard,ro /dev/sdc {mount_point}")
        print("Disk mounted")
# -




