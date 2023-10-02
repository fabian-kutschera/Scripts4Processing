import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

description = '''Plot Energy output(s). Default is Moment Rate (MR).)'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("EnergyOutput", help="energy.csv file(s) from your SeisSol simulation")
parser.add_argument("--Data", nargs="+", metavar=("variable"), default="MR", help="Data to plot. Default: Moment rate (MR). Other available options are Gravitational energy (GE), Acoustic energy (AE), Acoustic kinetic energy (AKE), Elastic kinetic energy (EKE), Elastic energy (EE), Total frictional work (TFW), Static frictional work (SFW), Seismic moment (SM), Plastic moment (PM)")
args = parser.parse_args()

df = pd.read_csv(args.EnergyOutput)

# print seismic moment and moment magnitude
M0 = df['seismic_moment'][len(df)-1]
Mw = 2/3 * (np.log10(M0) - 9.1)
print("M0 = {} and Mw = {}.".format(M0, Mw))

# plot chosen energy output
dt = df["time"][1]
moment_rate = np.gradient(df["seismic_moment"], dt)
plt.plot(df["time"], moment_rate, label="{:.2f}".format(Mw))
plt.legend()
plt.savefig("{:.2f}.png".format(Mw))
