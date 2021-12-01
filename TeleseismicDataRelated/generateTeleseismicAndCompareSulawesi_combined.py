import instaseis
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams["axes.axisbelow"] = False
from generateTeleseismicRoutines import *
from obspy.signal.trigger import recursive_sta_lta, trigger_onset
import numpy as np

def autoscale_y(ax,margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""


    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)

def InitializeSeveralStationsFigure(nplots, plot_time_axis_everywhere=False):
   ##for comparing signals on the same figure
   #figall, axarr = plt.subplots(nplots,3, figsize=(14, 9), dpi=160, sharex=False, sharey=False)
   figall, axarr = plt.subplots(nplots,3, figsize=(14, 12), dpi=160, sharex=False, sharey=False)

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
         axi.tick_params(axis='x', zorder=3)
         if i<nplots-1 and (not plot_time_axis_everywhere):
            axi.spines['bottom'].set_visible(False)
            axi.set_xticks([])
   return [figall, axarr]

def computerRMS(st,st_obs,t1,tmin, tmax):
   strace=st.select(component=comp)[0].copy()
   otrace=st_obs.select(component=comp)[0].copy()
   start_osTrace = max(strace.stats.starttime, otrace.stats.starttime)
   end_osTrace = min(strace.stats.endtime, otrace.stats.endtime)
   start_time_interp = max(t1+tmin, start_osTrace)
   npts_interp = min(int(end_osTrace-start_osTrace)+1, int(tmax)+1)
   strace.interpolate(sampling_rate=1.0,starttime=start_time_interp, npts=npts_interp)
   otrace.interpolate(sampling_rate=1.0,starttime=start_time_interp, npts=npts_interp)
   gof = np.std(strace.data-otrace.data)/np.std(otrace.data)
   del strace, otrace
   return gof

#db = instaseis.open_db('/export/data/ulrich/Instaseis/10s_PREM_ANI_FORCES')
#db = instaseis.open_db('/import/freenas-m-02-seismology/instaseis/PREMiso/')
db = instaseis.open_db('/import/freenas-m-05-seissol/ulrich/prem_a_10s_merge_compress2/')
#db = instaseis.open_db('syngine://prem_a_2s')
print(db)

c = Client("IRIS")
t1=UTCDateTime(2018, 9, 28, 10, 2, 43.60)

filename = 'PointSourceFile.h5'
sources0 = CreateInstaseisFiniteSource(db, filename, t1)


NetworkStationList = [['G','AIS'], ['II','PALK'], ['IU','TATO'],['IU','GUMO'],['IU','CTAO']]
#This is the network list for fig S12, S13
#NetworkStationList = [['G','ROCAM'],['II','COCO'],['II','NIL'],['IC','KMI'],['IC','ENH'],['IU','MAJO'],['IU','WAKE'],['IU','PMG'],['IU','MBWA'],['IU','NWAO'],]

nstations = len(NetworkStationList)
figall_P, axarr_P = InitializeSeveralStationsFigure(nstations, plot_time_axis_everywhere=True)
figall_whole, axarr_whole = InitializeSeveralStationsFigure(nstations)
components=['E','N','Z']

for ins, NetworkStation in enumerate(NetworkStationList):
   network=NetworkStation[0]
   station=NetworkStation[1]
   inventory = c.get_stations(network=network, station=station, level="channel", starttime=t1)
   print(inventory)
   stlat = inventory[0][0].latitude
   stlon = inventory[0][0].longitude
   stlat = geographic2geocentric(stlat)

   #create synthetic data with instaseis
   receiver = instaseis.Receiver(latitude=stlat, longitude=stlon, network=network, station=station)
   st0 = db.get_seismograms_finite_source(sources=sources0, receiver=receiver, kind='displacement')
   for i in range(3):
      st0[i].stats.starttime=t1
   t2 = st0[0].stats.endtime

   #retrieve teleseismic data
   st_obs0 = c.get_waveforms(network=network, station=station, location='00',
                               channel='LH*', starttime=t1, endtime=t2,
                               attach_response=True)
   st_obs0.rotate(method="->ZNE",inventory=inventory)

   # define a filter band to prevent amplifying noise during the deconvolution
   pre_filt = (0.0009, 0.001, 0.1, 0.3)
   #st_obs.remove_response(output='DISP', pre_filt=pre_filt)
   st_obs0.remove_response(output='DISP')

   #Pick Pwave arrival
   tr=st0.select(component='Z')[0]
   df = tr.stats.sampling_rate
   # Characteristic function and trigger onsets
   cft = recursive_sta_lta(tr.data, int(2.5 * df), int(20. * df))
   on_of = trigger_onset(cft, 6, 0.2)
   tmin = on_of[0, 0]/df - 20*df
   if station=='ROCAM':
      tmin=570
   tmax = tmin + 200
   #plt.plot(cft, 'k')
   #plt.show()

   ##########P wave plot #################################
   st_obs = st_obs0.copy()
   st = st0.copy()

   for myst in [st, st_obs]:
      myst.filter('bandpass', freqmin=0.00222222, freqmax=0.1, corners=4, zerophase=True)

   ymin0=1e10
   ymax0=-1e10
   for j,comp in enumerate(components):
      strace=st.select(component=comp)[0]
      otrace=st_obs.select(component=comp)[0]
      axarr_P[ins,j].plot(strace.times(reftime=t1),1e3*strace.data,'r')
      axarr_P[ins,j].plot(otrace.times(reftime=t1),1e3*otrace.data,'k')
      axarr_P[ins,j].set_xlim([tmin, tmax])
      autoscale_y(axarr_P[ins,j])
      ymin, ymax = axarr_P[ins,j].get_ylim()
      ymin0= min(ymin0, ymin)
      ymax0= max(ymax0, ymax)
      axarr_P[ins,j].set_title('')
   axarr_P[ins,0].set_ylabel('%s.%s' %(network,station))

   #Compute rMRS misfit and print it on plot
   for j,comp in enumerate(components):
      gof = computerRMS(st,st_obs,t1,tmin, tmax)
      axarr_P[ins,j].set_ylim([ymin0, ymax0])
      gofstring = 'rRMS %.2f' %(gof)
      axarr_P[ins,j].text(tmin+10, ymin0+0.15*(ymax0-ymin0), gofstring, fontsize=12)

   ##########whole signal  #################################
   st_obs = st_obs0.copy()
   st = st0.copy()

   for myst in [st, st_obs]:
      myst.filter('bandpass', freqmin=0.00222222, freqmax=0.015, corners=4, zerophase=True)

   for j,comp in enumerate(components):
      strace=st.select(component=comp)[0]
      otrace=st_obs.select(component=comp)[0]
      axarr_whole[ins,j].plot(strace.times(reftime=t1),1e3*strace.data,'r')
      axarr_whole[ins,j].plot(otrace.times(reftime=t1),1e3*otrace.data,'k')
      axarr_whole[ins,j].set_xlim([0, 1800])
   axarr_whole[ins,0].set_ylabel('%s.%s' %(network,station))

   #Compute rMRS misfit and print it on plot
   ymin0, ymax0 = axarr_whole[ins,0].get_ylim()
   for j,comp in enumerate(components):
      gof = computerRMS(st,st_obs,t1, 0., 1800.)
      gofstring = 'rRMS %.1f' %(gof)
      axarr_whole[ins,j].text(100, ymin0+0.15*(ymax0-ymin0), gofstring, fontsize=12)


figall_P.subplots_adjust(wspace=0.05)
figall_whole.align_ylabels(axarr_P[:, 0])
figall_P.savefig('teleseismic_comparison_Sulawesi_BodyWaves.svg', bbox_inches='tight', transparent=True)

figall_whole.subplots_adjust(wspace=0.05)
figall_whole.align_ylabels(axarr_whole[:, 0])
figall_whole.savefig('teleseismic_comparison_Sulawesi_whole_signal.svg', bbox_inches='tight')
