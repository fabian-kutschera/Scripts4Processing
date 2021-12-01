import argparse
import h5py

parser = argparse.ArgumentParser(description='compute nodal plate from Moment Tensor (h5) file')
parser.add_argument('filename', help='Moment Tensor (h5)')
args = parser.parse_args()

h5f = h5py.File(args.filename,'r')
aMomentTensor = h5f['MomentTensor'][:,:]

from obspy.imaging.beachball import mt2plane,MomentTensor
mt = MomentTensor(aMomentTensor[0,:], 1)

nodalplane = mt2plane(mt)
print(nodalplane.strike, nodalplane.dip, nodalplane.rake)
