from netCDF4 import Dataset
import numpy as np
import os
import argparse
from multiprocessing import Pool


def Writett3(kk):
    outputfile = inputfile + f"{kk}.tt3"
    ndt_chunk = len(my_chunks[kk])
    with open(outputfile, "w") as fid:
        if kk == 0:
            fid.write(f"{nx:7d} mx\n")
            fid.write(f"{ny:7d} my\n")
            fid.write(f"{ndt:7d} mt\n")
            fid.write(f"{xlower:14e} xlower\n")
            fid.write(f"{ylower:14e} ylower\n")
            fid.write(f"{t0:.14e} t0\n")
            fid.write(f"{dx:.14e} dx\n")
            fid.write(f"{dy:.14e} dy\n")
            fid.write(f"{dt:.14e} dt\n")
        for ik, k in enumerate(my_chunks[kk]):
            print(f"{kk}:{ik}/{ndt_chunk}")
            for j in range(ny):
                np.savetxt(fid, z[k, ny - j - 1, :], newline="\t", fmt="%.6e")
                fid.write("\n")
    print(f"done writing {outputfile}")


def chunk(xs, n):
    """Split the list, xs, into n evenly sized chunks"""
    L = len(xs)
    assert 0 < n <= L
    s, r = divmod(L, n)
    t = s + 1
    return [xs[p : p + t] for p in range(0, r * t, t)] + [
        xs[p : p + s] for p in range(r * t, L, s)
    ]


parser = argparse.ArgumentParser(
    description="convert a NetCDF displacement file to the tt3 format (GeoClaw)"
)
parser.add_argument("inputfile", help="input file (*.nc)")
parser.add_argument("--downsample", default=1, help="temporal downsample", type=int)
parser.add_argument(
    "--last", dest="last", action="store_true", help="use last time step only"
)
parser.add_argument(
    "--MP",
    nargs=1,
    metavar=("ncpu"),
    default=([1]),
    help="use np.pool to speed-up calculations",
    type=int,
)
args = parser.parse_args()


inputfile = args.inputfile

fh = Dataset(inputfile, mode="r+")
x = fh.variables["x"]
y = fh.variables["y"]
z = fh.variables["z"][:, :, :]
z = z[:: args.downsample, :, :]
ti = fh.variables["time"][:]
ti = ti[:: args.downsample]
print("ti=", ti)
t0 = ti[0]
dt = ti[1] - ti[0]
dx = x[1] - x[0]
dy = y[1] - y[0]
xlower = x[0]
ylower = y[0]

if args.last:
    z = z[-2:-1, :, :]
    print(z.shape)
    ti = [0]
    dt = 0
print("xlower, ylower", xlower, ylower)
ndt, ny, nx = z.shape
print("ndt, ny, nx", ndt, ny, nx)
print("done reading")
fh.close()


nproc = args.MP[0]
my_chunks = chunk(range(0, ndt), nproc)
with Pool(nproc) as p:
    p.map(Writett3, range(nproc))

fns = ""
for i in range(0, nproc):
    fns = fns + f"{inputfile}{i}.tt3 "
print("concatening files")
os.system(f"cat {fns} > cat{inputfile}.tt3")
print("removing temp files")
os.system(f"rm {fns}")
print(f"done writing cat{inputfile}.tt3")
