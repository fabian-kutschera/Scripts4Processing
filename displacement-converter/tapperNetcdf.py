from netCDF4 import Dataset
import numpy as np
import sys
import argparse
from shutil import copyfile


def CreateHann(nx, nha):
    hann = np.zeros(nx)
    hann[:] = 1.0
    ix = np.arange(0, nha, 1)
    hann1 = 0.5 * (1 - np.cos(np.pi * ix / nha))
    hann2 = np.flip(hann1)
    hann[0:nha] = hann1
    hann[nx - nha : nx] = hann2
    return hann


parser = argparse.ArgumentParser(
    description="apply a Hanning window on a NetCDF file. This prevents sharp displacement discontinuities at the limits of the region of imposed displacements, which could generate spurious waves"
)
parser.add_argument("inputfile", help="netcdf input file")
parser.add_argument(
    "--nHann",
    nargs=1,
    metavar=("nHann"),
    default=[250],
    help="length of the Hann window",
    type=int,
)
args = parser.parse_args()


inputfile = args.inputfile
outputfile = "tapered_" + args.inputfile

print("creating copy of input file")
copyfile(inputfile, outputfile)

fh = Dataset(outputfile, mode="r+")
z = fh.variables["z"]
ndt, ny, nx = z.shape

print("running swapaxes")
z = np.swapaxes(z, 0, 2)

hannx = CreateHann(nx, args.nHann[0])
for i in range(nx):
    if i % 100 == 0:
        print(i)
    z[i, :, :] = z[i, :, :] * hannx[i]

print("running swapaxes")
z = np.swapaxes(z, 0, 1)
hanny = CreateHann(ny, args.nHann[0])
for i in range(ny):
    if i % 100 == 0:
        print(i)
    z[i, :, :] = z[i, :, :] * hanny[i]
print("running swapaxes")
z = np.swapaxes(z, 0, 1)
print("running swapaxes")
z = np.swapaxes(z, 0, 2)
print("write z")
fh.variables["z"][:, :, :] = z[:, :, :]
fh.close()
