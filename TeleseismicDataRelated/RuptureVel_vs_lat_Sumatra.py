import numpy as np
import os
from computeMomentTensorRoutines import *
from numpy import linalg as LA
from functools import reduce

import argparse
import h5py
from submodules.pythonXdmfReader.pythonXdmfReader import *

parser = argparse.ArgumentParser(description='compute Moment Tensor from fault output file')
parser.add_argument('filename', help='fault model filename (xdmf)')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='project back to WGS84 from given projection string (ex: +init=EPSG:32646 (sumatra), +init=EPSG:3994 (NewZealand)')
parser.add_argument('--NZ', nargs=1, metavar=('nz'), default=([2]), help='number of point sources along z' ,type=int)
parser.add_argument('--NH', nargs=1, metavar=('nh'), default=([20]), help='number of point sources in the horizontal plane' ,type=int)
parser.add_argument('--vH', nargs=2, metavar=('vhx','vhy'), default=([1,1]), help='vector defining the slicing direction in the horizontal plane' ,type=float)

args = parser.parse_args()

if args.proj:
   import pyproj
   myproj=pyproj.Proj(args.proj[0])
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')


#read fault geometry and data from hdf5 file
prefix, ext = os.path.splitext(args.filename)

xyz0 = ReadGeometry(args.filename)
connect0 = ReadConnect(args.filename)
nElements = np.shape(connect0)[0]
ndt = ReadNdt(args.filename)

slip0,data_prec = LoadData(args.filename, 'ASl', nElements, ndt-1, True)
Vr0,data_prec = LoadData(args.filename, 'Vr', nElements, ndt-1, True)

print('computing area...')
cross0 = np.cross(xyz0[connect0[:,1],:] - xyz0[connect0[:,0],:], xyz0[connect0[:,2],:] - xyz0[connect0[:,0],:])
area0= 0.5*np.apply_along_axis(np.linalg.norm, 1, cross0)
print('computing barycenters...')
XYZcenters0 = (xyz0[connect0[:,0],:] + xyz0[connect0[:,1],:] + xyz0[connect0[:,2],:])/3.

IdnonZeroSlip = np.where(abs(slip0)>2.0)
#for defining the geometry subset
print(('%d slices along direction (horizontal): (%f,%f)' %(args.NH[0],args.vH[0],args.vH[1])))
xycenters0=args.vH[0]*XYZcenters0[:,0]+args.vH[1]*XYZcenters0[:,1]
xycenters00=args.vH[0]*XYZcenters0[IdnonZeroSlip,0]+args.vH[1]*XYZcenters0[IdnonZeroSlip,1]
xycmin = np.amin(xycenters00)
xycmax = np.amax(xycenters00)
dxyc = (xycmax-xycmin)
print(xycmin,xycmax)
myrange=np.linspace(0.,1,args.NH[0]+1)
dv = myrange[1]

#for defining the geometry subset
print(('%d slices along z' %(args.NZ[0])))
zcenters0=XYZcenters0[:,2]
zcenters00=XYZcenters0[IdnonZeroSlip,2]
zcmin = np.amin(zcenters00)
zcmax = np.amax(zcenters00)
print(zcmin,zcmax)
dzc = (zcmax-zcmin)
myZrange=np.linspace(0.,1,args.NZ[0]+1)
dz = myZrange[1]
print(myrange)

#write the momentTensor and MomentRate in a temp file
nsources = (myrange.shape[0]-1)*(myZrange.shape[0]-1)
print(nsources)
NormalizedMomentRate = np.zeros((nsources,ndt))
VrPer=np.zeros((nsources,3))
xyz=np.zeros((nsources,3))

M0allsrc=0
isrc=0
for v0 in myrange[0:-1]:
   for z0 in myZrange[0:-1]:
      idxys = np.where((xycenters0>(xycmin+v0*dxyc)) & (xycenters0<(xycmin+(v0+dv)*dxyc)))
      idxys = idxys[0]
      idzs = np.where((zcenters0>(zcmin+z0*dzc)) & (zcenters0<(zcmin+(z0+dz)*dzc)))
      idzs = idzs[0]
      ids = reduce(np.intersect1d, (idxys, idzs, IdnonZeroSlip))
      if ids.shape[0]==0:
          continue
      XYZcenters = XYZcenters0[ids,:]
      connect = connect0[ids,:]
      slip = slip0[ids]
      Vr=Vr0[ids]

      xyzc = [np.percentile(XYZcenters[:,0], 50), np.percentile(XYZcenters[:,1], 50), np.percentile(XYZcenters[:,2], 50)]
      print(xyzc)
      if args.proj:
         #project back to geocentric coordinates
         xyzc2 = pyproj.transform(myproj, lla, xyzc[0], xyzc[1], xyzc[2], radians=False)
         #z remains unchanged
         xyzc[0:2]=xyzc2[0:2]
      VrPer[isrc,0] = np.percentile(Vr, 50-34.1)
      VrPer[isrc,1] = np.percentile(Vr, 50)
      VrPer[isrc,2] = np.percentile(Vr, 50+34.1)
      print(xyzc[0], xyzc[1], np.percentile(Vr, 50-34.1), np.percentile(Vr, 50), np.percentile(Vr, 50+34.1) )
      xyz[isrc,:] = xyzc
      isrc=isrc+1

nz=3
zcmin = -30e3
zcmax = 0e3
dzc = (zcmax-zcmin)
myZrange=np.linspace(0.,1,args.NZ[0]+1)
dz = myZrange[1]
print(myrange)

VrVariations=np.zeros((args.NH[0],nz))
myZrange=np.linspace(0.,1,nz+1)
dz = myZrange[1]
xyz0=np.zeros((args.NH[0],nz,3))
for ii,v0 in enumerate(myrange[0:-1]):
   for jj,z0 in enumerate(myZrange[0:-1]):
      idxys = np.where((xycenters0>(xycmin+v0*dxyc)) & (xycenters0<(xycmin+(v0+dv)*dxyc)))
      idxys = idxys[0]
      idzs = np.where((zcenters0>(zcmin+z0*dzc)) & (zcenters0<(zcmin+(z0+dz)*dzc)))
      idzs = idzs[0]
      ids = reduce(np.intersect1d, (idxys, idzs, IdnonZeroSlip))
      if ids.shape[0]==0:
          continue
      XYZcenters = XYZcenters0[ids,:]
      connect = connect0[ids,:]
      slip = slip0[ids]
      Vr=Vr0[ids]

      xyzc = [np.percentile(XYZcenters[:,0], 50), np.percentile(XYZcenters[:,1], 50), np.percentile(XYZcenters[:,2], 50)]
      #print(xyzc)
      if args.proj:
         #project back to geocentric coordinates
         xyzc2 = pyproj.transform(myproj, lla, xyzc[0], xyzc[1], xyzc[2], radians=False)
         #z remains unchanged
         xyzc[0:2]=xyzc2[0:2]
      VrVariations[ii,jj] = np.percentile(Vr, 50)
      print(xyzc[0], xyzc[1], np.percentile(Vr, 50-34.1), np.percentile(Vr, 50), np.percentile(Vr, 50+34.1) )
      xyz0[ii,jj,:] = xyzc


import matplotlib.pyplot as plt
import matplotlib
modeltitle = 'dynamic rupture scenario'
modeltitle='this study,'
ps = 24
matplotlib.rcParams.update({'font.size': ps})
plt.rcParams['font.family']='sans'
matplotlib.rc('xtick', labelsize=ps)
matplotlib.rc('ytick', labelsize=ps)

fig = plt.figure(figsize=(16, 6), dpi=80)
ax = fig.add_subplot(111)

ids = np.where(VrPer[:,1]>0)[0]
xyz = xyz[ids,:]
VrPer = VrPer[ids,:]/1e3
VrVariations=VrVariations/1e3
plt.plot(xyz[:,1], VrPer[:,1], 'k', label='%s median' %(modeltitle))


#plt.xlim([0,600])

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlabel('latitude')
ax.set_ylabel('rupture speed (km/s)')

###Plotting Vallee rupture velo
#Line along which is computed the moment
a = np.loadtxt('data/lineVallee2.dat')
#Project line
import pyproj
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
myproj=pyproj.Proj('+init=EPSG:32646')
vs = pyproj.transform(lla, myproj, a[:,0], a[:,1], radians=False)
vs = np.asarray(vs)
#Compute curv. axis
diff = np.zeros(np.shape(vs))
diff[:,1:] = vs[:,1:]-vs[:, 0:-1]
faultlength = np.sqrt(np.square(diff[0,:])+np.square(diff[1,:]))/1e3
for i in range(1, np.shape(faultlength)[0]):
   faultlength[i]= faultlength[i]+faultlength[i-1]

from scipy import interpolate
f = interpolate.interp1d(faultlength, a[:,1], bounds_error=False, fill_value=3.3)
a = np.loadtxt('data/ruptureVelVallee.dat')
fa = [f(x) for x in a[:,0]]
plt.plot(fa, a[:,1], 'k-.', label='Vallee et al., 2007')

xi=[2,7,7.0001, 11,11.0001, 14]
yi=[2.7,2.7,2.5, 2.5,2,2 ]
plt.plot(xi,yi, 'k:', label='Guilbert et al., 2005')

cols='rgb'
print(myZrange)
for jj,z0 in enumerate(myZrange[0:-1]):
   zmin = zcmin+z0*dzc
   zmax = zcmin+(z0+dz)*dzc
   print(zmin, dzc, z0, dz, zmax)
   ids = np.where(VrVariations[:,jj]>0.9)[0]
   #added a min to avoid -0-10 in label
   plt.plot(xyz0[ids,jj,1], VrVariations[ids,jj], cols[jj]+'o', label='%s %.0f-%.0f km' %(modeltitle, min(0, -zmax/1e3), -zmin/1e3))

plt.legend(ncol=2)
plt.xlim(left=2, right=14.3)
x = np.arange(2,16,2)
labels = ['%d$^\circ$N' %(xi) for xi in list(x)]
y = np.arange(1,3.5,0.5)
ax.set_yticks(y)

ax.set_xticks(x)
ax.set_xticklabels(labels)

plt.ylim([0.5, 3.5])
plt.savefig('plots/RuptureVel.svg', bbox_inches='tight')



