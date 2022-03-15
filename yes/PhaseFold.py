#import and set up everything
import pandas as pd
import numpy as np
import sys
sys.path.append("../PhaseFold")
import os
import data
from data import loaders
import lightkurve as lk
import scipy.signal
import matplotlib.pyplot as plt
import math
from PIL import Image
import warnings
from ipywidgets.widgets import Button, Layout
from IPython.display import display
import functools


warnings.filterwarnings("ignore", category=RuntimeWarning) 

%load_ext autoreload
%autoreload 2
%matplotlib inline


# Mounting data
data_dir = "/mnt/disks/lcs/"
data.mount_drive(data_dir)
data_dir = "/mnt/disks/lcs/tess-goddard-lcs/"

def fold(ticid,sect):
    data_dir = "/mnt/disks/lcs/"
    data.mount_drive(data_dir)
    data_dir = "/mnt/disks/lcs/tess-goddard-lcs/"
    directory = f'S{sect}TIC{ticid}'
    parent_dir = f'/home/jupyter/PhaseFold/LightCurves' 
    path = os.path.join(parent_dir, directory)
    exist = os.path.isdir(path)
    if not exist:
        os.makedirs(path)
        
    global mini, maxi, midpt, lc
    ref = importsec(sect)
    tmcl = ["2_min_cadence" in fn for fn in ref.Filename]
    row = ref[ref.TIC_ID.isin([ticid])]
    string = str(row.Sector)
    s = string.split(" ")
    files = data_dir + ref.Filename.values
    filepath = files[int(s[0])]
    lc = loaders.load_lc(filepath)
    a = lc.scatter()
    pg = lc.normalize(unit='ppm').to_periodogram(oversample_factor=300)
    period = pg.period_at_max_power
    pg1 = lc.normalize(unit='ppm').to_periodogram(maximum_period = 2.1 * period.value, oversample_factor=100)
    fig = plt.figure(figsize = (12,5))
    mini = .7*period
    maxi = 1.3*period
    midpt = (mini + maxi)/2
    b = pg1.plot(view='period')
    midpt = redef()
    folded = lc.fold(midpt)
    cleanlightcurve = folded[folded.quality==0]
    c = cleanlightcurve.scatter(label=fr'Period = {midpt.value:.5f} d')
    
    b_save = Button (description = 'save', layout = Layout(width='100px'))
    def bsave(b_save):
        a.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png')
        #plt.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png')
        b.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png')
        c.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png')
        images = [Image.open(x) for x in [f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png']]
        widths, heights = zip(*(i.size for i in images))
        total_height = sum(heights)
        min_width = min(widths)
        new_im = Image.new('RGB', (min_width, total_height))
        y_offset = 0
        for im in images:
            new_im.paste(im, (0,y_offset))
            y_offset += im.size[1]
        #new_im.save(f'/home/jupyter/SPOcc/examples/Graphs/S{sect}TIC{ticid}.png')
        new_im.save(f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}.png')
    b_save.on_click(bsave)
    display(b_save)
    
# Importing the sector that the lightcurves are in
def importsec(sect):
    sec = sect
    ref = pd.read_csv(data_dir+"sector1lookup.csv")
    ref.head()
    ref[(ref.Camera==1)&(ref.CCD==1)&(ref.Magnitude>7)&(ref.Magnitude<9)].head()
    ref = pd.read_csv(data_dir+f"sector{sec}lookup.csv")
    return ref


# Folds and saves the original light curve, periodogram, folded light curve, and all 3 combined
# Into a folder in the LightCurves folder
def foldandsave(ticid,sect):
    global mini, maxi, midpt, lc
    directory = f'S{sect}TIC{ticid}'
    parent_dir = f'/home/jupyter/PhaseFold/LightCurves' 
    path = os.path.join(parent_dir, directory)
    exist = os.path.isdir(path)
    if not exist:
        os.makedirs(path)
    ref = importsec(sect)
    tmcl = ["2_min_cadence" in fn for fn in ref.Filename]
    row = ref[ref.TIC_ID.isin([ticid])]
    string = str(row.Sector)
    s = string.split(" ")
    files = data_dir + ref.Filename.values
    filepath = files[int(s[0])]
    lc = loaders.load_lc(filepath)
    a = lc.scatter()
    #a.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_LC.png')
    a.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png')
    pg = lc.normalize(unit='ppm').to_periodogram(oversample_factor=100)
    period = pg.period_at_max_power
    fig = plt.figure(figsize = (12,5))
    mini = .7*period
    maxi = 1.3*period
    midpt = (mini + maxi)/2
    xmax = 2.1 * period.value
    plt.xlim(0, xmax)
    plt.xlabel('Period [d]')
    plt.ylabel('Power [ppm]')
    plt.plot(pg.period, pg.power, color = 'black', markersize = 1)
    #plt.savefig(f'Individual_Images/S{sect}TIC{ticid}_Periodogram.png')
    plt.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png')
    midpt = redef()
    folded = lc.fold(midpt)
    cleanlightcurve = folded[folded.quality==0]
    c = cleanlightcurve.scatter(label=fr'Period = {midpt.value:.5f} d')
    #c.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_Folded.png')
    plt.close('all')
    c.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png')
    #images = [Image.open(x) for x in [f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Folded.png']]
    images = [Image.open(x) for x in [f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png']]
    widths, heights = zip(*(i.size for i in images))
    total_height = sum(heights)
    min_width = min(widths)
    new_im = Image.new('RGB', (min_width, total_height))
    y_offset = 0
    for im in images:
        new_im.paste(im, (0,y_offset))
        y_offset += im.size[1]
    #new_im.save(f'/home/jupyter/SPOcc/examples/Graphs/S{sect}TIC{ticid}.png')
    new_im.save(f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}.png')
    
#Saves the images as well as prints it out
def foldsaveprint(ticid,sect):
    global mini, maxi, midpt, lc
    directory = f'S{sect}TIC{ticid}'
    parent_dir = f'/home/jupyter/PhaseFold/LightCurves' 
    path = os.path.join(parent_dir, directory)
    exist = os.path.isdir(path)
    if not exist:
        os.makedirs(path)
    ref = importsec(sect)
    tmcl = ["2_min_cadence" in fn for fn in ref.Filename]
    row = ref[ref.TIC_ID.isin([ticid])]
    string = str(row.Sector)
    s = string.split(" ")
    files = data_dir + ref.Filename.values
    filepath = files[int(s[0])]
    lc = loaders.load_lc(filepath)
    a = lc.scatter()
    #a.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_LC.png')
    a.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png')
    pg = lc.normalize(unit='ppm').to_periodogram(oversample_factor=100)
    period = pg.period_at_max_power
    fig = plt.figure(figsize = (12,5))
    mini = .7*period
    maxi = 1.3*period
    midpt = (mini + maxi)/2
    xmax = 2.1 * period.value
    plt.xlim(0, xmax)
    plt.xlabel('Period [d]')
    plt.ylabel('Power [ppm]')
    plt.plot(pg.period, pg.power, color = 'black', markersize = 1)
    #plt.savefig(f'Individual_Images/S{sect}TIC{ticid}_Periodogram.png')
    plt.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png')
    midpt = redef()
    folded = lc.fold(midpt)
    cleanlightcurve = folded[folded.quality==0]
    c = cleanlightcurve.scatter(label=fr'Period = {midpt.value:.5f} d')
    #c.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_Folded.png')
    c.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png')
    #images = [Image.open(x) for x in [f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Folded.png']]
    images = [Image.open(x) for x in [f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png']]
    widths, heights = zip(*(i.size for i in images))
    total_height = sum(heights)
    min_width = min(widths)
    new_im = Image.new('RGB', (min_width, total_height))
    y_offset = 0
    for im in images:
        new_im.paste(im, (0,y_offset))
        y_offset += im.size[1]
    #new_im.save(f'/home/jupyter/SPOcc/examples/Graphs/S{sect}TIC{ticid}.png')
    new_im.save(f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}.png')

    
# Combines the seperate 3 pngs into one png
def combinefiles(ticid,sect):
    #images = [Image.open(x) for x in [f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/SPOcc/examples/Individual_Images/S{sect}TIC{ticid}_Folded.png']]
    images = [Image.open(x) for x in [f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png', f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png']]
    widths, heights = zip(*(i.size for i in images))

    total_height = sum(heights)
    max_width = max(widths)
    min_width = min(widths)

    new_im = Image.new('RGB', (min_width, total_height))

    y_offset = 0
    for im in images:
        new_im.paste(im, (0,y_offset))
        y_offset += im.size[1]

    #new_im.save(f'/home/jupyter/SPOcc/examples/Graphs/S{sect}TIC{ticid}.png')
    new_im.save(f'/home/jupyter/PhaseFold/LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}.png')

    
# Function to output original graph, periodogram, and phasefolded graph given TIC ID and their sector
def graphfoldprint(ticid,sect):
    global mini, maxi, midpt, lc
    ref = importsec(sect)
    tmcl = ["2_min_cadence" in fn for fn in ref.Filename]
    row = ref[ref.TIC_ID.isin([ticid])]
    string = str(row.Sector)
    s = string.split(" ")
    files = data_dir + ref.Filename.values
    filepath = files[int(s[0])]
    lc = loaders.load_lc(filepath)
    a = lc.scatter()
    pg = lc.normalize(unit='ppm').to_periodogram(oversample_factor=100)
    period = pg.period_at_max_power
    fig = plt.figure(figsize = (12,5))
    mini = .7*period
    maxi = 1.3*period
    midpt = (mini + maxi)/2
    xmax = 2.1 * period.value
    plt.xlim(0, xmax)
    plt.xlabel('Period [d]')
    plt.ylabel('Power [ppm]')
    plt.plot(pg.period, pg.power, color = 'black', markersize = 1)
    midpt = redef()
    folded = lc.fold(midpt)
    cleanlightcurve = folded[folded.quality==0]
    c = cleanlightcurve.scatter(label=fr'Period = {midpt.value:.5f} d')

# Function to save the original graph, periodogram, and phasefolded graph in the 
# LightCurves folder without combining them
def graphfoldsave(ticid,sect):
    data_dir = "/mnt/disks/lcs/"
    data.mount_drive(data_dir)
    data_dir = "/mnt/disks/lcs/tess-goddard-lcs/"
    directory = f'S{sect}TIC{ticid}'
    parent_dir = f'/home/jupyter/PhaseFold/LightCurves' 
    path = os.path.join(parent_dir, directory)
    exist = os.path.isdir(path)
    if not exist:
        os.makedirs(path)
        
    global mini, maxi, midpt, lc
    ref = importsec(sect)
    tmcl = ["2_min_cadence" in fn for fn in ref.Filename]
    row = ref[ref.TIC_ID.isin([ticid])]
    string = str(row.Sector)
    s = string.split(" ")
    files = data_dir + ref.Filename.values
    filepath = files[int(s[0])]
    lc = loaders.load_lc(filepath)
    a = lc.scatter()
    #a.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_LC.png')
    a.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_LC.png')
    pg = lc.normalize(unit='ppm').to_periodogram(oversample_factor=100)
    period = pg.period_at_max_power
    fig = plt.figure(figsize = (12,5))
    mini = 1.7*period
    maxi = 2.3*period
    midpt = (mini + maxi)/2
    xmax = 2.1 * period.value
    plt.xlim(0, xmax)
    plt.xlabel('Period [d]')
    plt.ylabel('Power [ppm]')
    plt.plot(pg.period, pg.power, color = 'black', markersize = 1)
    #plt.savefig(f'Individual_Images/S{sect}TIC{ticid}_Periodogram.png')
    plt.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Periodogram.png')
    midpt = redef()
    folded = lc.fold(midpt)
    c = folded.scatter(label=fr'Period = {midpt.value:.5f} d')
    #c.figure.savefig(f'Individual_Images/S{sect}TIC{ticid}_Folded.png')
    c.figure.savefig(f'LightCurves/S{sect}TIC{ticid}/S{sect}TIC{ticid}_Folded.png')

    
#cleanlightcurve - light curve without noise
#phasecurve - phasefolded light curve
#smoothcurve - median filtered version of phasecurve
#cleanlcmod - median filtered version of cleanlightcurve
#binary search within a range to find best period, A (min) B (max), use C(midpoint) and then compute A-C and B-C.
#then choose the smaller number and repeat. if same, B is best.

#adjusts min max and midpoint
def redef():
    global mini, maxi, midpt, lc
    while(mini.value + 0.0001 < maxi.value):
        minresstddev = calcresidualstddevmin(mini)
        midptresstddev = calcresidualstddevmidpt(midpt)
        maxresstddev = calcresidualstddevmax(maxi)
        if (minresstddev - midptresstddev) < (maxresstddev - midptresstddev):
            maxi = midpt
            midpt = (mini + maxi)/2
        else:
            mini = midpt
            midpt = (mini + maxi)/2
    
    #tests if a multiple of the midpt is better than the current one
    twomidptresstddev = calcresidualstddevmidpt(2*midpt)
    midptresstddev = calcresidualstddevmidpt(midpt)
    if (twomidptresstddev) < (midptresstddev):
        midpt = 2*midpt
    
    return midpt

#calculates the std dev of residuals of period given by argument (num)
def calcresidualstddevmin(num):
    lcfolded = lc.fold(num)
    cleanlightcurve = lcfolded[lcfolded.quality==0]
    phasecurve = lc.fold(num)[:]
    cleanlcmod = cleanlightcurve[:]
    cleanlcmod.flux = scipy.signal.medfilt(cleanlightcurve.flux, kernel_size=13)
    residual = cleanlcmod.flux - cleanlightcurve.flux
    ressqr = 0
    for x in range(len(cleanlcmod)):
        ressqr = ressqr + (residual[x] ** 2)
    minresstddev = math.sqrt((ressqr)/(len(cleanlcmod)-2))
    return minresstddev

def calcresidualstddevmax(num):
    lcfolded = lc.fold(num)
    cleanlightcurve = lcfolded[lcfolded.quality==0]
    phasecurve = lc.fold(num)[:]
    cleanlcmod = cleanlightcurve[:]
    cleanlcmod.flux = scipy.signal.medfilt(cleanlightcurve.flux, kernel_size=13)
    residual = cleanlcmod.flux - cleanlightcurve.flux
    ressqr = 0
    for x in range(len(cleanlcmod)):
        ressqr = ressqr + (residual[x] ** 2)
    maxresstddev = math.sqrt((ressqr)/(len(cleanlcmod)-2))
    return (maxresstddev)
    
def calcresidualstddevmidpt(num):
    lcfolded = lc.fold(num)
    cleanlightcurve = lcfolded[lcfolded.quality==0]
    phasecurve = lc.fold(num)[:]
    cleanlcmod = cleanlightcurve[:]
    cleanlcmod.flux = scipy.signal.medfilt(cleanlightcurve.flux, kernel_size=13)
    residual = cleanlcmod.flux - cleanlightcurve.flux
    ressqr = 0
    for x in range(len(cleanlcmod)):
        ressqr = ressqr + (residual[x] ** 2)
    midptresstddev = math.sqrt((ressqr)/(len(cleanlcmod)-2))
    return (midptresstddev)