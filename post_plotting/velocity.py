from mpl_toolkits.basemap import Basemap, maskoceans
import matplotlib.pyplot as plt
import matplotlib
from cmcrameri import cm
import argparse
import numpy as np

#run velocity.py --surface /import/freenas-m-05-seissol/kutschera/HIWI/fully-coupled/HFFZ/output_o4/HFFZ_fullycp_o4-surface.xdmf

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
fig.set_size_inches(15, 10)
ax1 = plt.gca()

#x1, x2, y1, y2 = -25.0, -12.0, 62.0, 68.0 #Iceland
x1, x2, y1, y2 = -19.0, -16.8, 65.9, 66.4 #Iceland

m = Basemap(epsg=4326, resolution="c", llcrnrlon=x1, urcrnrlon=x2, llcrnrlat=y1, urcrnrlat=y2)
bathy = "/import/freenas-m-05-seissol/kutschera/HIWI/Scripts/Github/GMT/Iceland/input/GEBCO_08_Dec_2021_074bdff0ecf9/gebco_2021_n68.0_s62.0_w-25.0_e-12.0.nc"

from netCDF4 import Dataset

fh = Dataset(bathy, mode="r")
y = fh.variables["lat"][:]
x = fh.variables["lon"][:]
#z = fh.variables["Band1"][:, :]
z = fh.variables["elevation"][:, :]

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
z[z < 0] = 0

plt.pcolormesh(X, Y, z, rasterized=True, cmap="gist_earth", vmin=-1000, vmax=2000, shading="auto")

if args.surface:
    import seissolxdmf

    xdmfFilename = args.surface
    sx = seissolxdmf.seissolxdmf(xdmfFilename)
    xyz = sx.ReadGeometry()
    connect = sx.ReadConnect()
    U = sx.ReadData("u3", sx.ndt - 1)

    if True:
        XYZcenters0 = (xyz[connect[:, 0], :] + xyz[connect[:, 1], :] + xyz[connect[:, 2], :]) / 3.0
        Zcenters0 = abs(XYZcenters0[:, 2])
        U = U[Zcenters0 < 0.001]
        connect = connect[Zcenters0 < 0.001]

    import pyproj
    
    transformer = pyproj.Transformer.from_crs("epsg:32627", "epsg:4326", always_xy=True)
    lon, lat = transformer.transform(xyz[:,0], xyz[:,1])

    #myproj = pyproj.Proj("EPSG:23839")
    #myproj = pyproj.Proj("EPSG:32627")

    #lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")
    #xyzb = pyproj.transform(myproj, lla, xyz[:, 0], xyz[:, 1], xyz[:, 2], radians=False)
    #xyzb = pyproj.transform(myproj, lla, xyz[:, 0], xyz[:, 1] + 1420000, xyz[:, 2], radians=False)


    #lons = xyzb[0]
    #lats = xyzb[1]
    #x, y = m(lons, lats)
    x, y = m(lon, lat)

    # n=100000
    # connect = connect[0:n]
    # U=U[0:n]
    plt.tripcolor(x, y, connect, facecolors=U, cmap=cm.cork, rasterized=True)
    plt.xlim(x1,x2)
    plt.xticks(np.arange(x1,x2, step=0.5))
    plt.ylim(y1,y2)
    plt.clim(-0.6, 0.6)


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


# m.drawmapscale(-20.0, 63.0, -20.0, 63.0, 40.0, barstyle="fancy", labelstyle="simple", fillcolor1="w", fillcolor2="#555555", fontsize=18, zorder=5)

if True:
    #x, y = m(x1 + 0.02, y1 + 0.3)
    x, y = m(-17.0, 66.04)

    plt.text(x, y, "vertical displacement (m)", va="center", fontsize=ps, rotation="vertical")
    # loc=4 means lower right
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    #cbaxes = inset_axes(ax1, width=0.15, height="50%", loc=4, borderpad=3)
    cbaxes = inset_axes(ax1, width=0.15, height="50%", loc=4)
    cbar = plt.colorbar(cax=cbaxes)#, ticks=np.linspace(-1, 1, 3))
    cbar.ax.tick_params(which="major", labelsize=ps, length=10, width=2, direction="inout")
    cbar.ax.tick_params(which="minor", length=6, width=1.0)
    cbar.ax.minorticks_on()

# draw parallels and meridians
# parallels = np.arange(y1, y2, 0.5)
# meridians = np.arange(118, 121, 0.5)

# m.drawparallels(parallels)#, labels=[1, 0, 0, 0], fontsize=ps, linewidth=0.2)
# m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=ps, linewidth=0.2)


# plot_fname = "vel.svg"
# plt.savefig(plot_fname, bbox_inches="tight", dpi=250)
# print("done writing %s" % (plot_fname))
# plt.show()
