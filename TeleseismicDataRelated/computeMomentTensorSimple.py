from math import pi,cos,sin
import numpy as np
def ComputeFaceMomentTensorNED(strike, dip, rake):
   cs = cos(strike)
   c2s = cos(2.*strike)
   cd = cos(dip)
   c2d = cos(2.*dip)
   cl = cos(rake)

   ss = sin(strike)
   s2s = sin(2.*strike)
   sd = sin(dip)
   s2d = sin(2.*dip)
   sl = sin(rake)

   M0 = 1.0
   MomentTensor = np.zeros((6))
   #0   1  2  3  4  5
   #xx,yy,zz,xy,xz,yz
   #http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/escidoc:65579/IS_3.8_rev1.pdf (eq 5)
   # with x y z : NED
   MomentTensor[0] = -M0*(sd*cl*s2s + s2d*sl*np.power(ss,2))
   MomentTensor[1] = M0*(sd*cl*s2s - s2d*sl*np.power(cs,2))
   MomentTensor[2] = M0*(s2d*sl)
   MomentTensor[3] = M0*(sd*cl*c2s + 0.5*s2d*sl*s2s)
   MomentTensor[4] = -M0*(cd*cl*cs + c2d*sl*ss)
   MomentTensor[5] = -M0*(cd*cl*ss-c2d*sl*cs)
   return MomentTensor

def ComputeFaceMomentTensorRTP(strike, dip, rake):
   # from obspy: https://docs.obspy.org/master/_modules/obspy/imaging/source.html
   # reoorder all moment tensors to NED and RTP convention
   # name : COMPONENT              : NED sign and index
   # NED  : NN, EE, DD, NE, ND, ED : [0, 1, 2, 3, 4, 5]
   # USE  : UU, SS, EE, US, UE, SE : [1, 2, 0, -5, 3, -4]
   # RTP  : RR, TT, PP, RT, RP, TP : [1, 2, 0, -5, 3, -4]
   # DSE  : DD, SS, EE, DS, DE, SE : [1, 2, 0, -5, -3, 4]
   #http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/escidoc:65579/IS_3.8_rev1.pdf (eq 5)
   #Mrr, Mtt, Mpp, Mrt, Mrp, Mtp in the spherical coordinate system defined by 
   #local upward vertical (r), North-South (t), and West-East (p) directions.
   mt = ComputeFaceMomentTensorNED(striker, dipr, raker)
   reorder = np.array([1, 2, 0, -5, 3, -4])
   mt_new = np.zeros((6))
   for i in range(6):
      n = reorder[i]
      if n>=0:
         mt_new[i] = mt[abs(n)]
      else:
         mt_new[i] = -mt[abs(n)]
   return mt_new

strike=-130.
dip=40.
rake=180

striker=strike*pi/180.
dipr=dip*pi/180.
raker=rake*pi/180.

mt = ComputeFaceMomentTensorNED(striker, dipr, raker)
print(mt)
from obspy.imaging.mopad_wrapper import beachball

beachball(mt, size=40, linewidth=2, facecolor='b',mopad_basis = 'NED', xy=(-70, 80))
print('ned',mt)
mt = [5.48666011e+19,-7.19747304e+19,1.71081293e+19,1.61560794e+19,-4.48619902e+19,-5.21000758e+19]


beachball(mt, size=40, linewidth=2, facecolor='b',mopad_basis = 'NED', xy=(-70, 80))
stop

mt0 = ComputeFaceMomentTensorRTP(striker, dipr, raker)
print('RT',mt0)
#beachball(mt0, size=40, linewidth=2, facecolor='b',mopad_basis = 'RT', xy=(-70, 80))
np1 = [strike, dip, rake]
#beachball(np1, mopad_basis = 'NED', xy=(0, 80))

import matplotlib.pyplot as plt
from obspy.imaging.beachball import beach
beach1 = beach(np1, xy=(-30, 0), width=30)
#beach2 = beach(mt, xy=(30, 0), width=30, mopad_basis = 'NED')
ax = plt.gca()
ax.add_collection(beach1) 
#ax.add_collection(beach2) 
ax.set_xlim((-60, 60))
ax.set_ylim((-30, 30))
ax.set_aspect("equal")
plt.show()
