import instaseis
import obspy
import h5py
from obspy.clients.fdsn.header import FDSNException
import os

def removeTopRightAxis(ax):
   #remove top right axis
   ax.spines['top'].set_visible(False)
   ax.spines['right'].set_visible(False)
   ax.get_xaxis().tick_bottom()
   ax.get_yaxis().tick_left()

def geographic2geocentric(lat):
   import numpy as np
   # geographic to geocentric
   # https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
   f = 1. / 298.257223563
   e2 = 2 * f - f ** 2
   return np.rad2deg(np.arctan((1 - e2) * np.tan(np.deg2rad(lat))))

def CreateInstaseisFiniteSource(db, filename, t1, myproj=''):
   #read HDF5 and create Finite Source for instaseis
   h5f = h5py.File(filename)
   NormalizedMomentRate = h5f['NormalizedMomentRate'][:,:]
   nsource, ndt = NormalizedMomentRate.shape
   xyz = h5f['xyz'][:,:]
   MomentTensor = h5f['MomentTensor'][:,:]
   dt = h5f['dt'][0]
   print(xyz)   
   if myproj!='':
      print('projecting back to geocentric')
      import pyproj
      myproj=pyproj.Proj(myproj)
      lla = pyproj.Proj(proj='latlong', ellps='sphere', datum='WGS84')
      txyz = pyproj.transform(myproj, lla, xyz[:,0], xyz[:,1], xyz[:,2], radians=False)
      xyz[:,0] = txyz[0]
      xyz[:,1] = txyz[1]
      #we use the same depth (else it is modified by the ellipsoid to sphere)
      #xyz[:,2] = txyz[2]
      print(xyz)
   else:
      print(filename)
      if h5f.attrs["CoordinatesConvention"]==b'geographic':
         xyz[:,1] = geographic2geocentric(xyz[:,1])
      elif h5f.attrs["CoordinatesConvention"]==b'geocentric':
         print('coordinates already in geocentric')
      else:
         print('coordinates are projected', h5f.attrs["CoordinatesConvention"])
         exit()
   h5f.close()
   lps= []
   for isrc in range(nsource):
      source = instaseis.Source( latitude=xyz[isrc,1], longitude=xyz[isrc,0], depth_in_m=-xyz[isrc,2],
         m_rr = MomentTensor[isrc,0], m_tt =MomentTensor[isrc,1], m_pp =MomentTensor[isrc,2],m_rt =MomentTensor[isrc,3],m_rp =MomentTensor[isrc,4],m_tp =MomentTensor[isrc,5],origin_time=t1, sliprate = NormalizedMomentRate[isrc,:], dt=dt)
      source.resample_sliprate(db.info.dt, int(ndt*dt/db.info.dt))
      lps.append(source)
   sources = instaseis.source.FiniteSource(pointsources=lps)
   return sources

def ImportObspyTraceGeonet(c, network,station,t1,t2):
   fname = 'observations/%s.mseed' %station
   if not os.path.exists('observations'):
      os.mkdir('observations')
   if os.path.isfile(fname):
      print('reading obs data from previous run...')
      st_obs = obspy.read(fname)
   else:
      print('requesting obs data...')
      # Fetch waveform from IRIS FDSN web service into a ObsPy stream object
      # and automatically attach correct response
      try:
         st_obs = c.get_waveforms(network=network, station=station, location='00',
                                     channel='LHZ', starttime=t1, endtime=t2,
                                     attach_response=True)
      except FDSNException:
         print(station + ":No data available for request.")
         return False
      st_obs.detrend('constant')
      # define a filter band to prevent amplifying noise during the deconvolution
      pre_filt = (0.0009, 0.001, 0.1, 0.3)
      st_obs.remove_response(output='DISP', pre_filt=pre_filt)
      st_obs.write(fname, format="MSEED")
      #gives similar results:
      #st_obs.remove_response(output="DISP")
   return st_obs
