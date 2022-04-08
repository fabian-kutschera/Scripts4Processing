from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt
import matplotlib
from cmcrameri import cm
import argparse
import numpy as np

plt.rcParams["font.family"] = "sans-serif"
ps = 20
plt.rcParams.update({"font.size": ps})
plt.rcParams["font.family"] = "sans"
matplotlib.rc("xtick", labelsize=ps)
matplotlib.rc("ytick", labelsize=ps)


parser = argparse.ArgumentParser(description="plot xyz ground displacement (obs or synth)")
parser.add_argument("--surface", help="surface output filename (xdmf), for synthetics")
parser.add_argument("--trace", dest="trace", action="store_true", help="plot fault trace")
args = parser.parse_args()

fig = plt.figure()
fig.set_size_inches(10, 10)
ax1 = plt.gca()

x1, x2, y1, y2 = 119.3, 120.2, -1.0, 0.2

m = Basemap(epsg=23839, resolution="c", llcrnrlon=x1, urcrnrlon=x2, llcrnrlat=y1, urcrnrlat=y2)


bathy = "/home/ulrich/work/SulawesiWL/WithBatnas_v1.0/BATNAS_MERGED_TRIMMED_v1.0_smooth.nc"
from netCDF4 import Dataset

fh = Dataset(bathy, mode="r")
y = fh.variables["lat"][:]
x = fh.variables["lon"][:]
z = fh.variables["Band1"][:, :]
fh.close()
idx = np.where((x > x1) & (x < x2))[0]
idy = np.where((y > y1) & (y < y2))[0]

x = x[idx]
y = y[idy]
z = z[idy, :]
z = z[:, idx]

lon, lat = np.meshgrid(x, y)
X, Y = m(lon, lat)

# Coast line
plt.contour(X, Y, z, levels=[0], colors="k")

# Bathymetry
z[z < 0] = np.nan

plt.pcolormesh(X, Y, z, rasterized=True, cmap="gist_earth", vmin=-1000, vmax=2000)

if args.surface:
    import seissolxdmf

    xdmfFilename = args.surface
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    U = sx.ReadData("w", sx.ndt - 1)

    if True:
        XYZcenters0 = (xyz[connect[:, 0], :] + xyz[connect[:, 1], :] + xyz[connect[:, 2], :]) / 3.0
        Zcenters0 = abs(XYZcenters0[:, 2])
        U = U[Zcenters0 < 0.001]
        connect = connect[Zcenters0 < 0.001]

    import pyproj

    myproj = pyproj.Proj("EPSG:23839")
    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    xyzb = pyproj.transform(myproj, lla, xyz[:, 0], xyz[:, 1] + 1420000, xyz[:, 2], radians=False)

    lons = xyzb[0]
    lats = xyzb[1]
    x, y = m(lons, lats)
    # n=100000
    # connect = connect[0:n]
    # U=U[0:n]
    plt.tripcolor(x, y, connect, facecolors=U, cmap=cm.broc, rasterized=True)
    plt.clim(-1, 1)


if args.trace:
    # Print fault trace from pl file
    import pyproj

    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    myproj = pyproj.Proj("EPSG:23839")
    # fid=open('segtrace.pl')
    # fid=open('straighttrace.pl')
    fid = open("data/straighttrace_noBay.pl")
    # fid=open('seg_ext90long.pl')
    lines = fid.readlines()
    fid.close()
    init = 1
    vs = []
    for line in lines:
        if line in ["ILINE\n", lines[-1]]:
            if init == 1:
                init = 0
            else:
                vs = np.asarray(vs)
                trace = pyproj.transform(myproj, lla, vs[:, 0], vs[:, 1], radians=False)
                xx, yy = m(trace[0], trace[1])
                m.plot(xx, yy, linewidth=1.5, color="k")
            vs = []
        if line.startswith("VR"):
            myitems = line.split()[2:5]
            vs.append([float(f) for f in myitems])


m.drawmapscale(119.60, -0.92, 119.60, -1.4, 40.0, barstyle="fancy", labelstyle="simple", fillcolor1="w", fillcolor2="#555555", fontsize=18, zorder=5)

if False:
    x, y = m(x1 + 0.02, y1 + 0.3)
    plt.text(x, y, "vertical displacement (m)", va="center", fontsize=ps, rotation="vertical")
    # loc=4 means lower right
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    cbaxes = inset_axes(ax1, width=0.1, height="40%", loc=1, borderpad=5)
    cbar = plt.colorbar(cax=cbaxes, ticks=np.linspace(-1, 1, 3))
    cbar.ax.tick_params(which="major", labelsize=ps, length=10, width=2, direction="inout")
    cbar.ax.tick_params(which="minor", length=6, width=1.0)
    cbar.ax.minorticks_on()

# draw parallels and meridians
parallels = np.arange(-2, 1, 0.5)
meridians = np.arange(118, 121, 0.5)

m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=ps, linewidth=0.2)
m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=ps, linewidth=0.2)


plot_fname = "vel.svg"
plt.savefig(plot_fname, bbox_inches="tight", dpi=250)
print("done writing %s" % (plot_fname))
# plt.show()
