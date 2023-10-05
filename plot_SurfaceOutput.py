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
parser.add_argument("xdmfFilename", nargs="+", metavar=("file"), 
                    help="surface.xdmf output from your SeisSol simulation.")
parser.add_argument("-d", "--Data", type=str, required=False, nargs="+", metavar=("variable"), default=["u3"], 
                    help="Data to plot. Default: vertical coseismic displacement u3.")
parser.add_argument("--xlim", nargs="+", type=float, help="xmin xmax")
parser.add_argument("--ylim", nargs="+", type=float, help="ymin ymax")
parser.add_argument("--xlabel", default="x [m]", type=str, help="Choose xlabel")
parser.add_argument("--ylabel", default="y [m]", type=str, help="Choose ylabel")
parser.add_argument("--cmap", default="bwr", help="Matplotlib colormap")
parser.add_argument("--dpi", default=300, type=int, help="Dots per Inch")
parser.add_argument("-o", "--Output", help='Name of the output file', default='surface_output.png')
args = parser.parse_args()


import pyproj



plt.xlabel(args.xlabel)
plt.ylabel(args.ylabel)