import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator

def getG(asagiFile, xyz):
   #read Asagi file and get G at coordinates
   fh = Dataset(asagiFile, mode='r')
   x  = fh.variables['x']
   y  = fh.variables['y']
   z  = fh.variables['z']
   data = fh.variables['data']
   data = np.swapaxes(data,0,2)
   data = [test[1] for test in data.flatten()]
   data = np.array(data)
   data = data.reshape((x.shape[0],y.shape[0],z.shape[0]))
   if z[1]<z[0]:
      data = np.flip(data, axis=2)
      z=np.flip(z)
   fn = RegularGridInterpolator((x,y,z), data)
   return fn(xyz)

