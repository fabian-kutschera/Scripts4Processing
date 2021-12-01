# =============================================================================
# Example how to run code from Console:
#     run plot_mao_cmt.py /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/TeleseismicDataRelated/PointSourceFile_1_1.h5
# 
# =============================================================================

import numpy as np
import argparse
import h5py
import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from obspy.imaging.mopad_wrapper import beach
from pyproj import Transformer
import cmt
from scalebar import scale_bar


def setup_map(ax, MapBoundaries, grid_spacing, gridlines_left=True):
    """Setup the map with cartopy"""
    ax.set_extent(MapBoundaries, crs=ccrs.PlateCarree())
    scale = "10m"
    ax.add_feature(cfeature.LAND.with_scale(scale))
    ax.add_feature(cfeature.OCEAN.with_scale(scale))
    ax.add_feature(cfeature.COASTLINE.with_scale(scale))
    ax.add_feature(cfeature.BORDERS.with_scale(scale), linestyle=":")
    locs = np.arange(-180, 180, grid_spacing)
    gl = ax.gridlines(draw_labels=True, ylocs=locs, xlocs=locs)
    gl.right_labels = False
    gl.top_labels = False
    gl.left_labels = gridlines_left


parser = argparse.ArgumentParser(
    description="draw map from multiple point source (h5)file"
)
parser.add_argument("filename", help="point source file (h5)")
parser.add_argument(
    "--proj",
    nargs=1,
    metavar=("projname"),
    help="transform to geocentric point source coordinates from projection projname",
)
parser.add_argument(
    "--x0y0proj",
    nargs=2,
    metavar=("x0", "y0"),
    default=[0, 0],
    help="offset of projection (e.g. used for centering UTM projection on faults",
)
parser.add_argument(
    "--dGrid",
    nargs=1,
    metavar=("dGrid"),
    default=[1.0],
    help="distance between consecutive parallel or meridians drawn",
    type=float,
)
parser.add_argument(
    "--beachSize",
    nargs=1,
    metavar=("dGrid"),
    default=[1.0],
    help="adjustement factor for beach ball sizes",
    type=float,
)

parser.add_argument(
    "--scalebarSize",
    nargs=1,
    metavar=("dGrid"),
    default=[100.0],
    help="size of scale bar in km",
    type=float,
)

parser.add_argument(
    "--MapBoundaries",
    nargs=4,
    metavar=("lonmin", "lonmax", "latmin", "latmax"),
    help="coordinates of map frame",
    type=float,
)
parser.add_argument(
    "--unicolor", dest="unicolor", action="store_true", help="beach balls all in blue"
)
args = parser.parse_args()

h5f = h5py.File(args.filename, "r")
aMomentTensor = h5f["MomentTensor"][:, :]
nsrc = aMomentTensor.shape[0]
xyz = h5f["xyz"][:, :]
h5f.close()
xyz[:, 0] += args.x0y0proj[0]
xyz[:, 1] += args.x0y0proj[1]

if args.proj:
    # epsg:4226 is geocentric (lat, lon)
    transformer = Transformer.from_crs(args.proj[0], "epsg:4326", always_xy=True)
    print("projecting back to WGS84")
    xyz[:, 0], xyz[:, 1] = transformer.transform(xyz[:, 0], xyz[:, 1])
    print((xyz[:, :]))

if not args.MapBoundaries:
    print("MapBoundaries not specified, inferring...")
    xmin = np.amin(xyz[:, 0])
    xmax = np.amax(xyz[:, 0])
    ymin = np.amin(xyz[:, 1])
    ymax = np.amax(xyz[:, 1])
    dx = max(ymax - ymin, xmax - xmin, 0.1)
    xmin = xmin - 0.2 * dx
    xmax = xmax + 0.2 * dx
    ymin = ymin - 0.2 * dx
    ymax = ymax + 0.2 * dx
    args.MapBoundaries = [xmin, xmax, ymin, ymax]
print(
    f"lon_min, lon_max {args.MapBoundaries[0:2]}, \
lat_min, lat_max {args.MapBoundaries[2:4]}"
)

dydx = (args.MapBoundaries[3] - args.MapBoundaries[2]) / (
    args.MapBoundaries[1] - args.MapBoundaries[0]
)

fig = plt.figure(figsize=(8, 8 * dydx), dpi=80)
ax = plt.axes(projection=ccrs.PlateCarree())
setup_map(ax, args.MapBoundaries, args.dGrid[0])
cmap = matplotlib.cm.get_cmap("Spectral")

for isrc in range(nsrc):
    M0all = cmt.compute_seismic_moment(aMomentTensor[isrc, :])
    Mw = 2.0 / 3.0 * np.log10(M0all) - 6.07
    print(f"{isrc:5d}: {Mw:.2f}", end=" ")
    if isrc % 10 == 9:
        print()
    MomentTensor = cmt.RTP2NED(aMomentTensor[isrc, :] / M0all)

    color = cmap(isrc / nsrc)
    if args.unicolor:
        color = "b"
    beach1 = beach(
        MomentTensor[:],
        xy=xyz[isrc, 0:2],
        width=Mw * args.beachSize[0] * 0.01,
        linewidth=0.2,
        mopad_basis="NED",
        facecolor=color,
    )
    ax.add_collection(beach1)

if not args.unicolor:
    norm = matplotlib.colors.Normalize(vmin=0, vmax=nsrc)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, shrink=0.7 / dydx, label="id source")

scale_bar(ax, (0.1, 0.1), args.scalebarSize[0])
fname = "BeachBallPlot.pdf"
plt.savefig(fname, bbox_inches="tight")
print(f"done writing {fname}")
