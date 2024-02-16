import argparse
import numpy as np
import os
import h5py
import pyproj
import seissolxdmf

parser = argparse.ArgumentParser(
    description="transform the geometry array of a SeisSol surface output file to the geocentric coordinate system"
)
parser.add_argument("xdmf_filename", help="seissol xdmf file")
parser.add_argument(
    "proj",
    nargs=1,
    metavar=("projname"),
    default=(""),
    help="string describing the origin projection",
)
args = parser.parse_args()

# initiate class
sx = seissolxdmf.seissolxdmf(args.xdmf_filename)
# load geometry array as a numpy array of shape ((nodes, 3))
geom = sx.ReadGeometry()

# set projection
lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
if args.proj[0] != "geocent":
    sProj = args.proj[0]
    myproj = pyproj.Proj(sProj)
else:
    myproj = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")

xyz = np.zeros_like(geom)
xyz[:, 0], xyz[:, 1], xyz[:, 2] = pyproj.transform(
    myproj, lla, geom[:, 0], geom[:, 1], geom[:, 2]
)

(
    dataLocation,
    data_prec,
    nElements,
    MemDimension,
) = sx.GetDataLocationPrecisionNElementsMemDimension("Geometry")
splitArgs = dataLocation.split(":")
isHdf5 = True if len(splitArgs) == 2 else False

mypath = os.path.dirname(args.xdmf_filename)
prefix = os.path.basename(args.xdmf_filename).split("-surface")[0]

if isHdf5:
    fname = os.path.join(mypath, prefix + "-surface_vertex_WGS84.h5")
    h5fv = h5py.File(fname, "w")
    h5fv.create_dataset("/mesh0/geometry", data=xyz)
    h5fv.close()
else:
    path2 = os.path.dirname(dataLocation)
    fname = os.path.join(mypath, path2, "geometry_WGS84.bin")
    output_file = open(fname, "wb")
    xyz.tofile(output_file)
    output_file.close()

print(f"done writing {fname}")

bn_fname = os.path.basename(fname)

fname_xdmf = os.path.join(mypath, "WGS84_" + os.path.basename(args.xdmf_filename))
with open(fname_xdmf, "w") as fout:
    with open(args.xdmf_filename) as fid:
        lines = fid.readlines()
    for line in lines:
        if isHdf5:
            text = line.replace(dataLocation, f"{bn_fname}:/mesh0/geometry")
        else:
            loca = os.path.join(path2, prefix + "-surface_vertex_WGS84.h5")
            text = line.replace(dataLocation, f"{loca}")
        fout.write(text)
print(f"done writing {fname_xdmf}")
