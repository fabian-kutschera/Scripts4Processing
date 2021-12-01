import instaseis
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
from generateTeleseismicRoutines import *

def InitializeSeveralStationsFigure(nplots):
   ##for comparing signals on the same figure
   figall, axarr = plt.subplots(nplots,3, figsize=(10, 6), dpi=160, sharex=False, sharey=False)

   direction = ['EW', 'NS', 'UD']
   axarr[nplots-1,2].set_xlabel('time (s)')
   for j in range(3):
      axi = axarr[0,j]
      axi.set_title(direction[j])
      axi = axarr[nplots-1,j]
      for i in range(nplots):
         axi = axarr[i,j]
         if j>0:
            axi.get_shared_y_axes().join(axi, axarr[i,0])
            axi.set_yticks([])
            axi.spines['left'].set_visible(False)
         removeTopRightAxis(axi)
         if i<nplots-1:
            axi.spines['bottom'].set_visible(False)
            axi.set_xticks([])
   return [figall, axarr]


db = instaseis.open_db('/export/data/ulrich/Instaseis/10s_PREM_ANI_FORCES')
db = instaseis.open_db('/import/freenas-m-02-seismology/instaseis/PREMiso/')
c = Client("IRIS")
t1=UTCDateTime(2016, 11, 13, 11, 2, 56.34)

filename = 'PointSourceFile.h5'
sources0 = CreateInstaseisFiniteSource(db, filename, t1)


NetworkStationList = [['IU','FUNA'], ['IU','NWAO'], ['G','PPTF'],['G','CCD'],['IU','PMG']]

nstations = len(NetworkStationList)
figall, axarr = InitializeSeveralStationsFigure(nstations)

for ins, NetworkStation in enumerate(NetworkStationList):
   network=NetworkStation[0]
   station=NetworkStation[1]
   inventory = c.get_stations(network=network, station=station, level="channel", starttime=t1)
   stlat = inventory[0][0].latitude
   stlon = inventory[0][0].longitude
   stlat = geographic2geocentric(stlat)

   #create synthetic data with instaseis
   receiver = instaseis.Receiver(latitude=stlat, longitude=stlon, network=network, station=station)
   st = db.get_seismograms_finite_source(sources=sources0, receiver=receiver, kind='displacement')
   for i in range(3):
      st[i].stats.starttime=t1
   t2 = st[0].stats.endtime

   #retrieve teleseismic data
   st_obs = c.get_waveforms(network=network, station=station, location='00',
                               channel='LH*', starttime=t1, endtime=t2,
                               attach_response=True)
   st_obs.rotate(method="->ZNE",inventory=inventory)

   #st_obs.remove_response(output="DISP")
   #st_obs.detrend('constant')
   # define a filter band to prevent amplifying noise during the deconvolution
   pre_filt = (0.0009, 0.001, 0.1, 0.3)
   #st_obs.remove_response(output='DISP', pre_filt=pre_filt)
   st_obs.remove_response(output='DISP')

   for myst in [st, st_obs]:
      #myst.filter('bandpass', freqmin=0.00222222, freqmax=0.025, corners=4, zerophase=True)
      #myst.filter('bandpass', freqmin=0.00222222, freqmax=0.01, corners=4, zerophase=True)
      myst.filter('bandpass', freqmin=0.00222222, freqmax=0.01, corners=4, zerophase=False)
   components=['E','N','Z']
   for j,comp in enumerate(components):
      strace=st.select(component=comp)[0]
      otrace=st_obs.select(component=comp)[0]
      axarr[ins,j].plot(strace.times(reftime=t1),1e3*strace.data,'b')
      axarr[ins,j].plot(otrace.times(reftime=t1),1e3*otrace.data,'k')
      if j==0:
         axarr[ins,j].set_ylabel('%s.%s' %(network,station))
      axarr[ins,j].set_xlim(0,1800)

#plt.show()
figall.subplots_adjust(wspace=0.05)
figall.savefig('teleseismic_comparison_NZ.svg', bbox_inches='tight')


