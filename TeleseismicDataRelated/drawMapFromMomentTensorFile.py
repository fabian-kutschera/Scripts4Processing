# =============================================================================
# Example how to run code from Console:
#    first run:
#       run computeMomentTensorSubSeisSol.py /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/complex_fault_geometry/Complex_Middle_M7.07/HFFtest-fault.xdmf --NZ 1 --NH 1 --invertSls --oneDvel /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/TeleseismicDataRelated/NIceland_1D.dat --proj '+proj=utm +zone=27'
#    then:
#       run drawMapFromMomentTensorFile.py ../output/PointSourceFile_1_1_Simple_Middle_M7.333.h5 --MapBoundaries -20 -16 65 67 --outName 'Simple_Middle_M7.333' --Title 'Centroid Simple Middle M7.33' 
# using output from Thomas new script: run drawMapFromMomentTensorFile.py ../../TuSeisSolScripts/TeleseismicDataRelated/PointSourceFile_1_1.h5 --MapBoundaries -20 -16 65 67
# =============================================================================

import numpy as np
import argparse
import h5py
from math import log10
parser = argparse.ArgumentParser(description='compute Moment Tensor from fault output file')
parser.add_argument('filename', help='Moment Tensor (h5)')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='project back to WGS84 from given projection string (ex: +init=EPSG:32646 (sumatra), +init=EPSG:3994 (NewZealand)')
parser.add_argument('--dGrid', nargs=1, metavar=('dGrid'), default = [1.0], help='distance between consecutive parallel or meridians drawn', type=float)
parser.add_argument('--beachSize', nargs=1, metavar=('dGrid'), default = [1.0], help='adjustement factor for beach ball sizes', type=float)
parser.add_argument('--MapBoundaries', nargs=4, metavar=('lonmin','lonmax','latmin','latmax'), help='coordinates of map frame', type=float)
parser.add_argument('--outName', nargs=1, metavar=('outName'), default = (''),  help='give name which will be appended to output file')
parser.add_argument('--Title', nargs=1, metavar=('Title'), default = False,  help='give title name')

args = parser.parse_args()

h5f = h5py.File(args.filename,'r')
aMomentTensor = h5f['MomentTensor'][:,:]
nsrc=aMomentTensor.shape[0]
xyz = h5f['xyz'][:,:]
h5f.close()

if args.proj:
   print('projecting back to WGS84')
   import pyproj
   myproj=pyproj.Proj(args.proj[0])
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   txyz = pyproj.transform(myproj, lla, xyz[:,0], xyz[:,1], xyz[:,2], radians=False)
   xyz[:,0] = txyz[0]
   xyz[:,1] = txyz[1]
   xyz[:,2] = txyz[2]
   print((xyz[:,:]))


if not args.MapBoundaries:
   print('MapBoundaries not specified, inferring...')
   xmin=np.amin(xyz[:,0])
   xmax=np.amax(xyz[:,0])
   ymin=np.amin(xyz[:,1])
   ymax=np.amax(xyz[:,1])
   dx = max(ymax-ymin, xmax-xmin,0.1)
   xmin = xmin-0.2*dx
   xmax = xmax+0.2*dx
   ymin = ymin-0.2*dx
   ymax = ymax+0.2*dx
else:
   xmin=args.MapBoundaries[0]
   xmax=args.MapBoundaries[1]
   ymin=args.MapBoundaries[2]
   ymax=args.MapBoundaries[3]

print(('lon_min, lon_max %f %f, lat_min, lat_max %f %f' %(xmin, xmax,ymin, ymax)))

###SETTING UP THE MAP BACKGROUND########
from mpl_toolkits.basemap import Basemap
#alternative if ModuleNotFoundError: No module named 'mpl_toolkits.basemap'
#import importlib
#importlib.import_module('mpl_toolkits.basemap').Basemap

import matplotlib.pyplot as plt
print('use source activate obspy')
#from obspy.imaging.beachball import beach
from obspy.imaging.mopad_wrapper import beach
fig, ax = plt.subplots() 
m = Basemap(projection='tmerc', resolution="h",
            llcrnrlon=xmin , urcrnrlon=xmax , llcrnrlat=ymin, urcrnrlat=ymax, lat_0 = 0.5*(ymin+ymax), lon_0 = 0.5*(xmin+xmax))
m.drawcoastlines()
#m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000, verbose= True)
m.drawmapboundary(fill_color='blue')
m.fillcontinents(color='wheat', lake_color='skyblue')
# draw parallels.
parallels = np.arange(-90,90,args.dGrid[0])
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8,linewidth=0.2)
# draw meridians
meridians = np.arange(0,360,args.dGrid[0])
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8,linewidth=0.2)

import matplotlib
cmap = matplotlib.cm.get_cmap('Spectral')

for isrc in range(nsrc):
   M0all = np.linalg.norm(aMomentTensor[isrc,:])
   Mw = 2./3.*log10(M0all)-6.07
   print(Mw)
   aMomentTensorRTP = aMomentTensor[isrc,:]/np.linalg.norm(aMomentTensor[isrc,:])
   ##convert RTP moment tensor to NED
   ## derived from source code of obspy
   ## https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
   signs = [1,1,1,-1,1,-1]
   indices = [1,2,0,5,3,4]
   aMomentTensorNED = np.array([sign * aMomentTensorRTP[ind] for sign, ind in zip(signs, indices)])

   beach1 = beach(aMomentTensorNED[:], xy=m(xyz[isrc,0],xyz[isrc,1]), width=max(1000, (Mw-5)*10000*args.beachSize[0]/2.), linewidth=0.2, mopad_basis = 'NED', facecolor=cmap(isrc/nsrc))
   #beach1.set_zorder(10)
   ax.add_collection(beach1)

# Draw a map scale
m.drawmapscale(
    xmin+0.2*(xmax-xmin), ymin+0.2*(ymax-ymin),
    xmin+0.2*(xmax-xmin), ymin+0.2*(ymax-ymin),
    30.,
    barstyle='fancy', labelstyle='simple',
    fillcolor1='w', fillcolor2='#555555',
    fontsize=8, units='km',
    zorder=5,
    linewidth=0.2)

if args.Title:
    plt.title('{}'.format(str(args.Title[0])))
    
if args.outName:
    plt.savefig('../output/BeachBallPlot_{}.png'.format(str(args.outName[0])), dpi=300, bbox_inches='tight')
else:
    plt.savefig('BeachBallPlot.png', dpi=300, bbox_inches='tight')

plt.show()
