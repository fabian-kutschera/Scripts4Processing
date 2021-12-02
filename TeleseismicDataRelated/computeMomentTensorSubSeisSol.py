# =============================================================================
# Example how to run code from Console:
#     run computeMomentTensorSubSeisSol.py /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/complex_fault_geometry/Complex_Middle_M7.07/HFFtest-fault.xdmf --NZ 1 --NH 1 --invertSls --oneDvel /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/TeleseismicDataRelated/NIceland_1D.dat --proj '+init=EPSG:32627'
#   or try (oneDvel)
#     run computeMomentTensorSubSeisSol.py /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/complex_fault_geometry/Complex_Middle_M7.07/HFFtest-fault.xdmf --NZ 1 --NH 1 --invertSls --oneDvel /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/TeleseismicDataRelated/NIceland_1D.dat --proj '+proj=utm +zone=27'
#   or try (asagiFile)
#     run computeMomentTensorSubSeisSol.py /import/freenas-m-05-seissol/kutschera/HIWI/SeisSol/complex_fault_geometry/Complex_Middle_M7.07/HFFtest-fault.xdmf --NZ 1 --NH 1 --invertSls --asagiFile /import/freenas-m-05-seissol/kutschera/HIWI/Scripts/TeleseismicDataRelated/vel_test.nc --proj '+proj=utm +zone=27'
# =============================================================================

import numpy as np
import os
from computeMomentTensorRoutines import *
from numpy import linalg as LA

import argparse
import h5py
from submodules.pythonXdmfReader.pythonXdmfReader import *

parser = argparse.ArgumentParser(description='compute Moment Tensor from fault output file')
parser.add_argument('filename', help='fault model filename (xdmf)')
#Note that it could be that we have a high SR sampled file output, and a coarse sampled (in time) file with the rest of the outputs
parser.add_argument('--STFfromSR', nargs=1, metavar=('xdmf SR File'), help='use SR to compute Source time function (high sampling rate required)')
parser.add_argument('--asagiFile', nargs=1, metavar=('asagiFile'), help='asagiFile containing the 3d material parameters (if not set use constant G=32e9)')
parser.add_argument('--oneDvel', nargs=1, metavar=('oneDvel'), help='1d velocity structure (col1: z col2: mu)')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='project back to WGS84 from given projection string (ex: +init=EPSG:32646 (sumatra), +init=EPSG:3994 (NewZealand)')
parser.add_argument('--NZ', nargs=1, metavar=('nz'), default=([2]), help='number of point sources along z' ,type=int)
parser.add_argument('--NH', nargs=1, metavar=('nh'), default=([20]), help='number of point sources in the horizontal plane' ,type=int)
parser.add_argument('--ndt', nargs=1, metavar=('ndt'), help='use a subset of time frames' ,type=int)
parser.add_argument('--vH', nargs=2, metavar=('vhx','vhy'), default=([1,1]), help='vector defining the slicing direction in the horizontal plane' ,type=float)
parser.add_argument('--static', dest='static', action='store_true', help='static fault output (only 1 dt)')
parser.add_argument('--invertSls', dest='invertSls', action='store_true', help='invert Sls (problem in normal definition)')
parser.add_argument('--invertSld', dest='invertSld', action='store_true', help='invert Sld (problem in normal definition)')

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

if not args.static:
   if args.ndt==[]:
      ndt=args.ndt[0]
      print('using ndt=', ndt)
   else:
      ndt = ReadNdt(args.filename)
else:
   ndt=1

#ndt, nElements
Sls,data_prec = LoadData(args.filename, 'Sls', nElements, ndt-1, True)
Sld,data_prec = LoadData(args.filename, 'Sld', nElements, ndt-1, True)
if args.invertSls:
   Sls=-Sls
if args.invertSld:
   Sld=-Sld
rake0=np.arctan2(Sld, Sls)

Sls = None
Sld = None

slip0,data_prec = LoadData(args.filename, 'ASl', nElements, ndt-1, True)

print('computing area...')
cross0 = np.cross(xyz0[connect0[:,1],:] - xyz0[connect0[:,0],:], xyz0[connect0[:,2],:] - xyz0[connect0[:,0],:])
area0= 0.5*np.apply_along_axis(np.linalg.norm, 1, cross0)
cross0 = None

print('computing barycenters...')
XYZcenters0 = (xyz0[connect0[:,0],:] + xyz0[connect0[:,1],:] + xyz0[connect0[:,2],:])/3.

if (not args.asagiFile) and (not args.oneDvel):
  print('using constant G')
  G = np.zeros((nElements))
  G[:]=32e9
elif args.asagiFile:
  print('loading G from asagi file')
  from getG import *
  G= getG(args.asagiFile[0], XYZcenters0)
else:
  print('loading 1d velocity file')
  depthmu = np.loadtxt(args.oneDvel[0])
  depthmu = depthmu[depthmu[:,0].argsort()]
  print(depthmu)
  G = np.interp(XYZcenters0[:,2], depthmu[:,0], depthmu[:,1])
  print(G)

Garea0= G*area0
G=None
area0=None

if not args.STFfromSR:
  print('computing moment rate function from absolute slip...')
  if not args.static:
    dt = ReadTimeStep(args.filename)
  else:
    dt=1.0
  FaceMomentRate0=np.zeros((ndt, nElements))
  #if too many elements this may generate a memory overflow, therefore the if
  if ndt*nElements<10e6:
     ASl,data_prec = LoadData(args.filename, 'ASl', nElements)
     FaceMomentRate0[1:,:] = np.diff(ASl, axis=0)/dt
  else:
     print('using a slower but more memory efficient implementation')
     slip1, data_prec = LoadData(args.filename, 'ASl', nElements, 0, True)
     for i in range(1,ndt):
        print(i)
        slip0 = np.copy(slip1)
        slip1 , data_prec = LoadData(args.filename, 'ASl', nElements, i, True)
        FaceMomentRate0[i,:] =(slip1-slip0)/dt
  FaceMomentRate0 = Garea0* FaceMomentRate0
  ASl = None
  FaceMoment0 = None

else:

  ndt = ReadNdt(args.STFfromSR[0])
  print('computing moment rate function from SR time histories')
  SRs,data_prec = LoadData(args.STFfromSR[0], 'SRs', nElements)
  SRd,data_prec = LoadData(args.STFfromSR[0], 'SRd', nElements)
  nElements1 = np.shape(SRs)[0]
  if (nElements1!=nElements):
     print(('main output file (%d) and SR file (%d) do not have the same number of elements' %(nElements, nElements1)))
  dt = ReadTimeStep(args.STFfromSR[0])
  FaceMomentRate0=np.zeros(SRs.shape)
  SR = np.sqrt(np.power(SRs,2)+np.power(SRd,2))
  SRs = None
  SRd = None
  FaceMomentRate0[:,:] = Garea0*SR
  SR = None

  xyz1 = ReadGeometry(args.STFfromSR[0])
  connect1 = ReadConnect(args.STFfromSR[0])
  XYZcenters1 = (xyz1[connect1[:,0],:] + xyz1[connect1[:,1],:] + xyz1[connect1[:,2],:])/3.
  xycenters1=args.vH[0]*XYZcenters1[:,0]+args.vH[1]*XYZcenters1[:,1]
  zcenters1=XYZcenters1[:,2]

#compute various geometric object needed by Rcalc
print('computing strike and dip angles...')
strike0,dip0 = ComputeStrikeDip(xyz0, connect0)

#Useful if large portion of the fault not ruptured
#idsRupture = np.where(slip0>0.01)[0]
#XYZcenters0 = XYZcenters0[idsRupture,:]

#for defining the geometry subset
print(('%d slices along direction (horizontal): (%f,%f)' %(args.NH[0],args.vH[0],args.vH[1])))
xycenters0=args.vH[0]*XYZcenters0[:,0]+args.vH[1]*XYZcenters0[:,1]
xycmin = np.amin(xycenters0)
xycmax = np.amax(xycenters0)
dxyc = (xycmax-xycmin)

myrange=np.linspace(0.,1,args.NH[0]+1)
dv = myrange[1]

#for defining the geometry subset
print(('%d slices along z' %(args.NZ[0])))
zcenters0=XYZcenters0[:,2]
zcmin = np.amin(zcenters0)
zcmax = np.amax(zcenters0)
dzc = (zcmax-zcmin)
myZrange=np.linspace(0.,1,args.NZ[0]+1)
dz = myZrange[1]
print(myrange)

#write the momentTensor and MomentRate in a temp file
nsources = (myrange.shape[0]-1)*(myZrange.shape[0]-1)
print(nsources)
NormalizedMomentRate = np.zeros((nsources,ndt))
MomentTensor=np.zeros((nsources,6))
xyz=np.zeros((nsources,3))

M0allsrc=0
isrc=0
for v0 in myrange[0:-1]:
   for z0 in myZrange[0:-1]:
      if not args.STFfromSR:
         idxys = np.where((xycenters0>(xycmin+v0*dxyc)) & (xycenters0<(xycmin+(v0+dv)*dxyc)))
         idxys = idxys[0]
         idzs = np.where((zcenters0>(zcmin+z0*dzc)) & (zcenters0<(zcmin+(z0+dz)*dzc)))
         idzs = idzs[0]
         ids = np.intersect1d(idxys, idzs)
      else:
         idxys = np.where((xycenters1>(xycmin+v0*dxyc)) & (xycenters1<(xycmin+(v0+dv)*dxyc)))
         idxys = idxys[0]
         idzs = np.where((zcenters1>(zcmin+z0*dzc)) & (zcenters1<(zcmin+(z0+dz)*dzc)))
         idzs = idzs[0]
         ids = np.intersect1d(idxys, idzs)
 
      XYZcenters = XYZcenters0[ids,:]
      connect = connect0[ids,:]
      Garea = Garea0[ids]
      slip = slip0[ids]
      rake = rake0[ids]
      strike= strike0[ids]
      dip=dip0[ids]

      FaceMomentRate = FaceMomentRate0[:,ids]
      FaceMoment = np.trapz(FaceMomentRate, axis=0)*dt
      MomentRate = np.sum(FaceMomentRate, axis=1)
      Mom = np.trapz(MomentRate)*dt

      if (abs(Mom)<1e-3) & (not args.static): 
         continue

      xyzc = np.average(XYZcenters, axis=0, weights=FaceMoment)
      MomentRate = MomentRate/Mom
      NormalizedMomentRate[isrc,:] = MomentRate
      nElements = np.shape(strike)[0]
      FaceMomentTensor = ComputeFaceMomentTensorNED(strike, dip, rake, Garea, slip)
      aMomentTensor = computeMomentTensor(FaceMomentTensor)
      M0all = LA.norm(aMomentTensor, 2)
      M0allsrc = M0allsrc + M0all
      MomentTensor[isrc,:] =  aMomentTensor
      if args.proj:
         #project back to geocentric coordinates
         xyzc2 = pyproj.transform(myproj, lla, xyzc[0], xyzc[1], xyzc[2], radians=False)
         #z remains unchanged
         xyzc[0:2]=xyzc2[0:2]
      xyz[isrc,:] = xyzc
      isrc=isrc+1

print('Mw(allsrc)=', 2./3.*log10(M0allsrc)-6.07)


##convert NED moment tensor to RTP
## derived from source code of obspy
## https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
signs = [1,1,1,1,-1,-1]
indices = [2,0,1,4,5,3]

aMomentTensorRTP = np.array([sign * MomentTensor[:,ind] for sign, ind in zip(signs, indices)])
MomentTensor = np.transpose(aMomentTensorRTP)


nsources = isrc
h5f = h5py.File('../output/PointSourceFile_%d_%d.h5' %(args.NH[0],args.NZ[0]),'w')
h5f.create_dataset('NormalizedMomentRate', (nsources,ndt), dtype='d')
h5f.create_dataset('xyz', (nsources,3), dtype='d')
h5f.create_dataset('MomentTensor', (nsources,6), dtype='d')
h5f.create_dataset('dt', (1,), dtype='d')
h5f['dt'][0]=dt
h5f['NormalizedMomentRate'][:,:] = NormalizedMomentRate[0:nsources,:]
h5f['MomentTensor'][:,:] = MomentTensor[0:nsources,:]
h5f['xyz'][:,:] = xyz[0:nsources,:]
if args.proj:
   h5f.attrs["CoordinatesConvention"] = np.string_("geographic")
else:
   h5f.attrs["CoordinatesConvention"] = np.string_("projected")
h5f.close()
