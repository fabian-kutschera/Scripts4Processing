import numpy as np
import argparse
import h5py
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='compute Moment Tensor from fault output file')
parser.add_argument('filename', help='Moment Tensor (h5)')
parser.add_argument('--idSTF', nargs='+', help='list of the STF to visualize (1st = 0); ' ,type=int)
parser.add_argument('--normalized', dest='normalized', action='store_true', help='show normalized stf (integral=1)')
parser.add_argument('--cumulative', dest='cumulative', action='store_true', help='show a cumulative plot of all STF')
parser.add_argument('--MRF', nargs=1, metavar=('MomentRateFile'), help='save MomentRate to file')
args = parser.parse_args()

h5f = h5py.File(args.filename,'r')
STF = h5f['NormalizedMomentRate'][:,:]
print((STF.shape))
dt = h5f['dt'][0]
nsources, ndt = STF.shape

if not args.normalized:
   aMomentTensor = h5f['MomentTensor'][:,:]
   Mom = np.linalg.norm(aMomentTensor, axis=1)
   c = np.diag(Mom)
   STF = np.dot(c,STF)
h5f.close()

if args.idSTF==None or args.cumulative:
   print('plotting all STF')
   args.idSTF=list(range(0,nsources))

time = np.linspace(0,(ndt-1)*dt,ndt)
print(time)
print((time.shape))
print((STF.shape))
cols = 'bgrcykb'

import matplotlib
cmap = matplotlib.cm.get_cmap('Spectral')

for ids,i in enumerate(args.idSTF):
   if ((args.cumulative) and (i>0)):
      STF[i,:] = STF[i-1,:] + STF[i,:]
   #plt.plot(time, STF[i,:], cols[ids%7])
   plt.plot(time, STF[i,:], color = cmap(ids/nsources))
plt.xlim([0,ndt*dt])
plt.show()

if args.MRF!=None:
   print(('saving moment rate function to %s' %args.MRF[0]))
   if not args.cumulative:
      for i in range(0,nsources):
         STF[i,:] = STF[i-1,:] + STF[i,:]
   np.savetxt(args.MRF[0], np.column_stack((time,STF[nsources-1,:]))) 
