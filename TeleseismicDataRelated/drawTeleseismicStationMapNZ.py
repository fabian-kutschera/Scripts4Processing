
import obspy
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np

c = Client("IRIS")
t1=UTCDateTime(2016, 11, 13, 11, 2, 56.34)

#NetworkStationList = [['IU','FUNA'], ['IU','JOHN'], ['IU','NWAO'], ['G','PPTF'],['G','CCD'],['IU','PMG']]
NetworkStationList = [['IU','FUNA'], ['IU','NWAO'], ['G','PPTF'],['G','CCD'],['IU','PMG']]
nstations = len(NetworkStationList)

for ins, NetworkStation in enumerate(NetworkStationList):
   network=NetworkStation[0]
   station=NetworkStation[1]
   inventory = c.get_stations(network=network, station=station, level="station", starttime=t1)
   stlat = inventory[0][0].latitude
   stlon = inventory[0][0].longitude
   


from mpl_toolkits.basemap import Basemap
m = Basemap(projection='ortho',lat_0=-42,lon_0=173.5,resolution='c')

m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-90., 120., 30.))
m.drawmeridians(np.arange(0., 420., 60.))
m.drawmapboundary()

x=[]
y=[]
names=[]
for ins, NetworkStation in enumerate(NetworkStationList):
   network=NetworkStation[0]
   station=NetworkStation[1]
   inventory = c.get_stations(network=network, station=station, level="station", starttime=t1)
   stlat = inventory[0][0].latitude
   stlon = inventory[0][0].longitude
   xi, yi = m(stlon, stlat)
   x.append(xi)
   y.append(yi)
   name = '%s.%s' %(network, station)
   names.append(name)
   print(name, stlat, stlon)

print(x,y,names)

m.scatter(x, y, 200, color='r', marker="v", edgecolor="k", zorder=3)
for i in range(len(names)):
    plt.text(x[i], y[i], names[i], va="top", family="monospace", weight="bold", clip_on=True)
plt.savefig('TeleseismicStations.svg', bbox_inches='tight')
plt.show()

