#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:56:36 2019

@author: anthony
"""

# Most recent edition for plotting the Multi-Frequency data AND the gating data
# Previous version was BasicMultiPlotv2_1

# This is looking for data files with the name **over**pqratio.dat
# The 'over' is the important bit

# One large running issue is that the fitting doesn't always play nice with the data
# In particular, when the denominator is even, there is almost no net motion
# Not only is it very small and hard to fit, it's also very complex
# The fit looks nothing like the data
# Ignore it, because we're just going to use v_max anyways

###############################################################################
# Section 0: Import the libraries
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.filedialog as fd
import scipy.optimize as opt
import os
import re

###############################################################################
# Section 1: Just about everything else
###############################################################################

# This one isn't quite as pretty as I had in previous versions
# That also means it's a lot less complicated though

root = tk.Tk()
root.withdraw()
# Start by opening a window that leads us to the data
# I'm assuming only one folder of data right now
startdir = fd.askdirectory(initialdir='/home/',title='Open the folder with data files')
# Change to that directory
os.chdir(startdir)
# Now pull a list of everything in that folder
startlist = os.listdir(startdir)
filelist = []
# Sort it out so that we only take the pieces in our naming scheme
for i in range(0,len(startlist)):
    if 'over' in startlist[i]:
        # Save all of those to a list of the data files only
        filelist.append(startlist[i])
        
# This is the fitting function we'll use to match the data
# This follows what Brown used
def fitfunc(x,A,B,phi0):
    return A*np.sin(x*B - (np.pi*.5) - phi0)
parray = []
qarray = []
datalist = []
#colorset=['r','b','m','g','c']
# Now we'll run through each file and do everything in one swing
for i in range(0,len(filelist)):
    # Load the file in a numpy array
    temparray=np.genfromtxt(filelist[i])
    # Search for the p and q values and save them as integers
    rs = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",filelist[i])
    p = int(rs[0])
    q = int(rs[1])
    parray.append(p)
    qarray.append(q)
    # Plug in the data to a fitting program
    fitvar,fiterr = opt.curve_fit(fitfunc,temparray[:,0],temparray[:,1],p0=(1,q,0))
    # Use the fitted parameters to make a fitted plot
    fittedcurve = fitfunc(temparray[:,0],fitvar[0],fitvar[1],fitvar[2])
    # Temp values for saving into a tuple
    tmpval1 = p/q
    tmpval2 = p*q
    # This tuple will get thrown into an array later
    # In order, it's p/q, then p*q, then max velocity
    temptuple = (tmpval1,tmpval2,max(temparray[:,1]))
    datalist.append(temptuple)
    # Finally plot one data file at a time
    plt.plot(temparray[:,0],temparray[:,1],'o',label='%s over %s' % (p,q))
    plt.plot(temparray[:,0],fittedcurve)
    ax = plt.gca()
    plt.xlabel(r'$Phase\ \phi\ (radians)$')
    ax.set_xticks([0,.5*np.pi,np.pi,1.5*np.pi,2*np.pi])
    ax.set_xticklabels(["$0$",r"$\frac{1}{2}\pi$",r"$\pi$",r"$\frac{3}{2}\pi$",r"$2\pi$"])
    plt.ylabel(r'$Average\ Velocity\ \langle v\rangle /v_r$')
    plt.title('Average velocity of multi-frequency ratchet')
    plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
    plt.show()
#plt.show()
# Take that information from the tuples/datalist and plot them
dataarray = np.array(datalist)
# This is the max velocity vs. p/q graph
plt.stem(dataarray[:,0],dataarray[:,2],'-b',alpha=0.3,markerfmt='bo',basefmt='w')
for i in range(0,len(qarray)):
    plt.text(dataarray[i,0]-.05,dataarray[i,2]+0.3,r'$\frac{%s}{%s}$' % (parray[i],qarray[i]))
plt.ylim(0,8)
plt.title(r'$Velocity\ versus\ p/q\ ratio$')
plt.xlabel(r'$p/q$')
plt.ylabel(r'$Velocity\ v_{max}$')
plt.show()
# This is the max velocity vs. p*q graph
plt.plot(dataarray[:,1],dataarray[:,2],'bo')
#for i in range(0,len(qarray)):
    #plt.text(dataarray[i,1]+0.02,dataarray[i,2]+0.02,r'%s/%s' % (parray[i],qarray[i]))
plt.xlabel(r'$p*q$')
plt.ylabel(r'$Velocity\ v_{max}$')
plt.title(r'$Velocity\ versus\ p*q\ ratio$')
plt.show()
    
root.quit()
