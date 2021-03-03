# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:00:45 2021
main file of PyStat, v201
+ Time series generation

@author: MM42910
"""
# Imports
import numpy as np
import time, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pystats_functions as fc

# Initialisation
# Aruments reading
arguments = sys.argv
case_short = arguments[1]
csv_folder = arguments[2]
img_folder = arguments[3]
options = arguments[4]
n_avg = int(arguments[5])

##Arguments given
# case_short = 'A1'
# csv_folder = 'csv'
# img_folder = 'img'
# options = 'rR'
# n_avg = 3

# File detection -> to function generPath
case_name = fc.generPath(case_short, csv_folder, img_folder)
   
## Bin details
bin_width = 100;
bin_length = 1000;
bins = np.arange(bin_width/2, bin_length+bin_width/2, bin_width)

## Font size
SMALL_SIZE = 14
MEDIUM_SIZE = 20
BIGGER_SIZE = 25

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)

## Start timer
tic = time.perf_counter()
tic_loop = time.perf_counter()

# Reading loop and snapshots generation
#  function processData that takes into arguments (Csvname, options) and
#  print evolution+timer
timeseries = fc.processData(case_name, bins, csv_folder, img_folder, n_avg, 
                            options)

# End
#  print timer
toc = time.perf_counter()
print(case_name+f" completed in {toc - tic:0.4f} s")  