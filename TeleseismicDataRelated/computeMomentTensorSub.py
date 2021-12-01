import numpy as np
import os
from computeMomentTensorRoutines import *

import argparse
import h5py

parser = argparse.ArgumentParser(description='compute Moment Tensor from fault output file')
parser.add_argument('filename', help='fault model filename (xdmf)')
parser.add_argument('--drawplot', dest='drawplot', action='store_true', help='plot moment tensors on a map')
parser.add_argument('--invertSls', dest='invertSls', action='store_true', help='invert Sls (problem in normal definition)')
parser.add_argument('--invertSld', dest='invertSld', action='store_true', help='invert Sld (problem in normal definition)')
args = parser.parse_args()

#read fault geometry and data from hdf5 file
prefix, ext = os.path.splitext(args.filename)

h5f = h5py.File(prefix+'.h5')
xyz = h5f['mesh0/geometry'][:,:]
connect0 = h5f['mesh0/connect'][:,:]
area0 = h5f['mesh0/area'][:]
slip0=h5f['mesh0/ASl'][:]
Sls = h5f['mesh0/Sls'][:]
Sld = h5f['mesh0/Sld'][:]

if args.invertSls:
   Sls=-Sls
if args.invertSld:
   Sld=-Sld
rake0=np.arctan2(Sld, Sls)

nElements = np.shape(area0)[0]
G = np.zeros((nElements))
G[:]=32e9
Garea0 = G*area0
h5f.close()

XYZcenters0 = (xyz[connect0[:,0],:] + xyz[connect0[:,1],:] + xyz[connect0[:,2],:])/3.

#compute various geometric object needed by Rcalc
strike0,dip0 = ComputeStrikeDip(xyz, connect0)

#for defining the geometry subset
xycenters0=XYZcenters0[:,0]+XYZcenters0[:,1]
xycmin = np.amin(xycenters0)
xycmax = np.amax(xycenters0)
dxyc = (xycmax-xycmin)

if args.drawplot:
   ###SETTING UP THE MAP BACKGROUND########
   from mpl_toolkits.basemap import Basemap
   import matplotlib.pyplot as plt
   print('use source activate obspy')
   #from obspy.imaging.beachball import beach
   from obspy.imaging.mopad_wrapper import beach
   fig, ax = plt.subplots() 
   m = Basemap(epsg = 3994, resolution="h",
               llcrnrlon=172.5, urcrnrlon=174.5, llcrnrlat=-42.7, urcrnrlat=-41.2)
   m.drawcoastlines()
   #m.arcgisimage(service='World_Shaded_Relief', xpixels = 1000, verbose= True)
   m.drawmapboundary(fill_color='aqua')
   m.fillcontinents(color='wheat', lake_color='skyblue')

   import pyproj
   myproj=pyproj.Proj('+init=EPSG:3994')
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

myrange=np.linspace(0.,1,6)
dv = myrange[1]
print(myrange)
#myrange = [0.8,1]

for v0 in myrange[0:-1]:
   ids = np.where((xycenters0>(xycmin+v0*dxyc)) & (xycenters0<(xycmin+(v0+dv)*dxyc)))
   ids = ids[0]
   XYZcenters = XYZcenters0[ids,:]
   xyzc = np.median(XYZcenters, axis=0)
   #project back to WGS84
   xyzll = pyproj.transform(myproj, lla, xyzc[0], xyzc[1], xyzc[2], radians=False)
   connect = connect0[ids,:]
   Garea = Garea0[0][ids]
   slip = slip0[0][ids]
   rake = rake0[0][ids]
   strike= strike0[ids]
   dip=dip0[ids]
   #print np.median(strike*180/3.14),np.median(dip*180/3.14),np.median(rake*180/3.14), np.median(Sls), np.median(Sld)

   FaceMomentTensor = ComputeFaceMomentTensorNED(strike, dip, rake, Garea, slip)
   aMomentTensor = computeMomentTensor(FaceMomentTensor)
   print(aMomentTensor)

   if args.drawplot:
      beach1 = beach(aMomentTensor, xy=m(xyzll[0],xyzll[1]), width=10000, linewidth=1, mopad_basis = 'NED')
      #beach1.set_zorder(10)
      ax.add_collection(beach1)


if args.drawplot:
   plt.savefig('test.png', bbox_inches='tight')
   plt.show()

