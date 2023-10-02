import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot Energy ouput (default: Moment Rate)")
parser.add_argument("EnergyOutput", help="energy.csv file from your SeisSol simulation")
parser.add_argument("--Data", nargs="+", metavar=("variable"), default="MR", help="Data to plot. Default: Moment rate (MR). Other available options are Gravitational energy (GE), Acoustic energy (AE), Acoustic kinetic energy (AKE), Elastic kinetic energy (EKE), Elastic energy (EE), Total frictional work (TFW), Static frictional work (SFW), Seismic moment (SM), Plastic moment (PM)")

df = pd.read_csv(args.EnergyOutput)

M0 = df['seismic_moment'][len(df)-1]
Mw = 2/3 * (np.log10(M0) - 9.1)
print("M0 = {} and Mw = {}.".format(M0, Mw))
