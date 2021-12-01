import numpy as np
import os
from computeMomentTensorRoutines import *
from io import BytesIO
import h5py
from math import pow
def createHatFunction(x,tc,td):
   for i, xi in enumerate(x):
      if xi<tc-0.5*td:
         y[i] = 0
      elif xi<tc:
         y[i] = (xi-(tc-0.5*td))/(0.5*td)
      elif xi<tc+0.5*td:
         y[i] = 1.0-(xi-tc)/(0.5*td)
      else:
         y[i] = 0.0
   return y

import argparse

parser = argparse.ArgumentParser(description='compute Moment Tensor from Multiple CMT Analysis inversions')
parser.add_argument('filename', help='Multiple CMT Analysis inversions file')
parser.add_argument('--tdconstant', nargs=1, metavar=('tdconstant'), default=([30]), help='use constant td for STF (half duration)' ,type=float)
parser.add_argument('--tdformula', dest='tdformula', action='store_true', help='compute td from M0 as in Tsai 2005')
parser.add_argument('--STFsampling', nargs=1, metavar=('STFsampling'), default=([0.5]), help='sampling rate of the STF (in s)' ,type=float)
args = parser.parse_args()

mydata = np.genfromtxt(args.filename, delimiter=" ")
print(mydata)
nsources = mydata.shape[0]
print('nb sources %d' %nsources)

xyz = np.zeros((nsources,3))
xyz[:,0]=mydata[:,1]
xyz[:,1]=mydata[:,0]
xyz[:,2]=-mydata[:,8]*1e3

"""
print('projecting back to WGS84 sphere')
import pyproj
myproj=pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='sphere', datum='WGS84')
txyz = pyproj.transform(myproj, lla, xyz[:,0], xyz[:,1], xyz[:,2], radians=False)
#this should only affect the latitude
xyz[:,0] = txyz[0]
xyz[:,1] = txyz[1]
"""

strike =np.radians(mydata[:,3])
dip =np.radians(mydata[:,4])
rake =np.radians(mydata[:,5])

M0=mydata[:,6]*1e20
Garea = np.zeros((nsources,))
Garea[:]=1
slip=M0

aMomentTensor = ComputeFaceMomentTensorNED(strike, dip, rake, Garea, slip)
aMomentTensor = np.transpose(aMomentTensor)

##convert NED moment tensor to RTP
## derived from source code of obspy 
## https://github.com/obspy/obspy/blob/master/obspy/imaging/source.py#L115
signs = [1,1,1,1,-1,-1]
indices = [2,0,1,4,5,3]

aMomentTensorRTP = np.array([sign * aMomentTensor[:,ind] for sign, ind in zip(signs, indices)])
MomentTensor = np.transpose(aMomentTensorRTP)

print((MomentTensor.shape))

#tc is the center of the hat
#td is the half width of the hat
tc= mydata[:,2]
td= np.zeros(tc.shape)

if args.tdformula:
   td = 2.2e-8*np.power(M0*1e7,1./3.)
else:
   td[:]=args.tdconstant[0]
   print('using constant td %f s' %td[0])

print('tc and td array:',(tc,td))

x=np.arange(0,np.amax(tc+td),args.STFsampling[0])
print(x)
ndt=np.shape(x)[0]
dt= x[1]
NormalizedMomentRate = np.zeros((nsources,ndt))

y=np.zeros(x.shape)
for i in range(nsources):
   NormalizedMomentRate[i,:]= createHatFunction(x, tc[i], td[i])
   Mom = np.trapz(NormalizedMomentRate[i,:])*dt
   NormalizedMomentRate[i,:]=NormalizedMomentRate[i,:]/Mom

from os.path import basename
prefix, ext = os.path.splitext(basename(args.filename))

h5f = h5py.File('PointSourceFile_%s.h5' %(prefix),'w')
h5f.create_dataset('NormalizedMomentRate', (nsources,ndt), dtype='d')
h5f.create_dataset('xyz', (nsources,3), dtype='d')
h5f.create_dataset('MomentTensor', (nsources,6), dtype='d')
h5f.create_dataset('dt', (1,), dtype='d')
h5f.attrs["CoordinatesConvention"] = np.string_("geographic")
h5f['dt'][0]=dt
h5f['NormalizedMomentRate'][:,:] = NormalizedMomentRate[0:nsources,:]
h5f['MomentTensor'][:,:] = MomentTensor[0:nsources,:]
h5f['xyz'][:,:] = xyz[0:nsources,:]

h5f.close()
