#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 15:06:19 2021

@author: fkutschera and tulrich
"""

# =============================================================================
# Example how to run code from Console:
#     run momentrate.py
#
# =============================================================================


import h5py
import numpy as np
import matplotlib.pyplot as plt

labelSeisSol=['base scenario']
for i,fn in enumerate(['base']):
   #fname = f'PointSourceFile_1_1_{fn}.h5'
   #fname = f'PointSourceFile_20_2.h5'
   fname = f'PointSourceFile_1_1.h5'
   h5f = h5py.File(fname,'r')
   STF = h5f['NormalizedMomentRate'][:,:]
   aMomentTensor = h5f['MomentTensor'][:,:][0]
   # generate full Moment tensor to compute M0
   fullMomentTensor = np.zeros((3, 3))
   fullMomentTensor[0, :] = [aMomentTensor[0], aMomentTensor[3], aMomentTensor[4]]
   fullMomentTensor[1, :] = [aMomentTensor[3], aMomentTensor[1], aMomentTensor[5]]
   fullMomentTensor[2, :] = [aMomentTensor[4], aMomentTensor[5], aMomentTensor[2]]
   Mom = np.linalg.norm(fullMomentTensor) * np.sqrt(0.5)
   STF=Mom*STF
   dt = h5f['dt'][0]
   nsources, ndt = STF.shape
   time = np.linspace(0,(ndt-1)*dt,ndt)
   #print(time, STF[0,:])
   
# =============================================================================
#    a = np.loadtxt(f'PlastMoment_{fn}.dat')
#    plastSTF = np.interp(time, a[:,0], a[:,1])
#    #plt.plot(time, STF[0,:], 'b'+ls[i], label=labelSeisSol[i])
#    plt.plot(time, STF[0,:]+plastSTF, 'b'+ls[i], label=labelSeisSol[i])
# =============================================================================

   plt.plot(time, STF[0,:], label=labelSeisSol[i])# 'b'+ls[i], label=labelSeisSol[i])
   plt.xlim(0,20)
   plt.grid()
