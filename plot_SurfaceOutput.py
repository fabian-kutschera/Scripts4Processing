#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:18:55 2023

@author: fkutschera
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import seissolxdmf 

description = '''Plot SeisSol surface output.'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("xdmfFilename", metavar=("file"), 
                    help="surface.xdmf output from your SeisSol simulation.")
parser.add_argument("-d", "--Data", type=str, required=False, nargs="+", metavar=("variable"), default=["u3"], 
                    help="Surface output variable to plot. Default: vertical coseismic displacement u3.")
parser.add_argument("--last", default=True, help="Plot last time step (default).")
parser.add_argument("--idt", help="Plot selected time step. Note: time step != time.", type=int)
parser.add_argument("--bounds", nargs="+", help="xmin xmax ymin ymax", type=float)
parser.add_argument("--proj", nargs=1, metavar=("projname"), help="project the data. projname: PROJ string describing the projection (ex epsg:32646 for UTM46N). Use geocent for geocentric cartesian")
parser.add_argument("--xlabel", default="x [m]", type=str, help="Choose xlabel")
parser.add_argument("--ylabel", default="y [m]", type=str, help="Choose ylabel")
parser.add_argument("--cmap", default="bwr", help="Matplotlib colormap")
parser.add_argument("--dpi", default=300, type=int, help="Dots per Inch")
parser.add_argument("-o", "--Output", help='Name of the output file', default='surface_output.png')
args = parser.parse_args()

if os.path.isfile(os.path.expanduser(args.xdmfFilename)):
    print("File selected: {}".format(os.path.expanduser(args.xdmfFilename)))
else:
    sys.exit("File '{}' does not exist.".format(os.path.expanduser(args.xdmfFilename)))

sx = seissolxdmf.seissolxdmf(args.xdmfFilename) # initiate class
xyz = sx.ReadGeometry() # load geometry array as a numpy array of shape
print(xyz)
connect = sx.ReadConnect() # load connectivity array as a numpy array of shape
# coords = xyz[connect].mean(axis=1)
# print(coords)
# x, y = coords.T[0], coords.T[1]
# print(x,y)

if (args.idt == True):
    ndt = args.idt
    print("Reading time step {}".format(ndt))
else:
    ndt = sx.ReadNdt()
    print("Reading last time step, i.e., {}".format(ndt))
    
    
# Plot data
surface_output = sx.ReadData("u1", ndt-1)
plt.tripcolor(xyz.T[0], xyz.T[1], connect, cmap="seismic", facecolors=surface_output, rasterized=True)
    
# if args.bounds:
#     idx = ((x > args.bounds[0]) & (x < args.bounds[1]) &
#           (y > args.bounds[2]) & (y < args.bounds[3]))
    
# if args.xlabel and args.ylabel:
#     plt.xlabel(args.xlabel)
#     plt.ylabel(args.ylabel)
# elif (args.proj == False):
#     plt.xlabel("x [m]")
#     plt.ylabel("y [m]")
    

print("Done.")

    
# output = sx.ReadData(args.Data, ndt)
# import pyproj
# plt.xlabel(args.xlabel)
# plt.ylabel(args.ylabel)