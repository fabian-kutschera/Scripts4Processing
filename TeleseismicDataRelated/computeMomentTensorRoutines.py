from math import pi,cos,sin,log10
import numpy as np
def ComputeFaceMomentTensorNED(strike, dip, rake, Garea, slip):
   nElements = np.shape(strike)[0]

   cs = np.cos(strike)
   c2s = np.cos(2.*strike)
   cd = np.cos(dip)
   c2d = np.cos(2.*dip)
   cl = np.cos(rake)

   ss = np.sin(strike)
   s2s = np.sin(2.*strike)
   sd = np.sin(dip)
   s2d = np.sin(2.*dip)
   sl = np.sin(rake)

   M0 = Garea*slip
   M0all = np.sum(M0, axis=0)

   MomentTensor = np.zeros((6,nElements))
   #0   1  2  3  4  5
   #xx,yy,zz,xy,xz,yz
   #http://gfzpublic.gfz-potsdam.de/pubman/item/escidoc:65580/component/escidoc:65579/IS_3.8_rev1.pdf (eq 5)
   # with x y z : NED
   MomentTensor[0,:] = -M0*(sd*cl*s2s + s2d*sl*np.power(ss,2))
   MomentTensor[1,:] = M0*(sd*cl*s2s - s2d*sl*np.power(cs,2))
   MomentTensor[2,:] = M0*(s2d*sl)
   MomentTensor[3,:] = M0*(sd*cl*c2s + 0.5*s2d*sl*s2s)
   MomentTensor[4,:] = -M0*(cd*cl*cs + c2d*sl*ss)
   MomentTensor[5,:] = -M0*(cd*cl*ss-c2d*sl*cs)
   return MomentTensor

def computeMomentTensor(FaceMomentTensor):
   return np.sum(FaceMomentTensor, axis=1)

def ComputeStrikeDip(xyz, connect):
   nElements = np.shape(connect)[0]
   #compute triangle normals
   un = np.cross(xyz[connect[:,1],:]-xyz[connect[:,0],:],xyz[connect[:,2],:]-xyz[connect[:,0],:])
   norm=np.apply_along_axis(np.linalg.norm, 1, un)
   un = un/norm.reshape((nElements,1))
   un=un.T
   refVector=[1e-5, 0, 1]
   un[:,:] = un[:,:] * np.sign(np.dot(un.T, refVector))
   #compute Strike and dip direction
   us = np.zeros(un.shape)
   us[0,:] = -un[1,:]
   us[1,:] = un[0,:]
   norm=np.apply_along_axis(np.linalg.norm, 0, us)
   us = us/norm
   ud = np.cross(un.T, us.T).T
 
   strike = np.arctan2(us[0,:], us[1,:])
   idsNumericalError = np.where(ud[2,:]>1)
   ud[2,idsNumericalError]=1
   dip = np.arcsin(ud[2,:])
   return (strike, dip)

#### These latter two routines are simplified (routines just for testing)
def ComputeFaceMomentTensorNEDSimple(strike, dip, rake):
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

def ComputeFaceMomentTensorRTPSimple(strike, dip, rake):
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

