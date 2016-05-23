from __future__ import division

from tools.convolution import convolution
from tools.extrema import extrema
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import sys

import argparse

parser = argparse.ArgumentParser(description='Compute the smoothing tree of a halo.')
parser.add_argument('--particles', type=str, help='Path to the particles file', default='tmp/536_around_clean')
parser.add_argument('--bins', type=int, help='Number of bins', default=64)
parser.add_argument('--cpus', type=int, help='Number of cpus to use', required=True)
parser.add_argument('--nsigmas', type=int, help='Number of sigmas to use', required=True)
parser.add_argument('--out', '-o', type=str, help='Name of the output file', required=True)
args = parser.parse_args()
print(args)


bins = args.bins
extrema.nbins = np.array([bins, bins, bins])
extrema.npeaks = bins**3
extrema.nprocs = args.cpus
NSIGMAS = args.nsigmas
sigmas = np.linspace(1, bins, NSIGMAS)

#df = pd.read_csv('tmp/536_around', delim_whitespace=True, engine='c', index_col=0)
# dfcleaned = df[df.x > 0.1][df.y > 0.01][df.z > -0.126]
# df = dfcleaned
try:
    df = pd.read_hdf(args.particles)
except:
    df = pd.read_csv(args.particles, delim_whitespace=True, engine='c', index_col=0)

pos = df[['x', 'y', 'z']].as_matrix().T

dens, edges = convolution.conv_density(pos, bins)

dfft = convolution.fft(dens)

def smooth(s):
    gaussian = convolution.kernel_gaussian3d(bins, s)
    gfft = convolution.fft(gaussian)
    cfft = convolution.conv_prod(gfft, dfft)
    conv = convolution.ifft(cfft)
    return conv

def extr(arr):
    _, npeak = extrema.extrema_compute(arr)
    index_bn, peakpos, eigvectors, eigvalues, peaktype = extrema.extrema_get(3, npeak)

    ind = np.unravel_index(index_bn, extrema.nbins, order='F')

    return dict(ind=ind, peakpos=peakpos, eigvectors=eigvectors, eigvalues=eigvalues, peaktype=peaktype)

def convert_to_df(res):

    return pd.DataFrame(dict(i=res['ind'][0], j=res['ind'][1], k=res['ind'][2],
                             x=res['peakpos'][0], y=res['peakpos'][1], z=res['peakpos'][2],
                             evx=res['eigvectors'][0], evy=res['eigvectors'][1], evz=res['eigvectors'][2],
                             ev=res['eigvalues'],
                             type=res['peaktype'])
    )

def extr_as_df(arr):
    res = extr(arr)

    return convert_to_df(res)


extrs = pd.DataFrame()

def doOne(sigma):
    # smooth
    smoothed = smooth(sigma)

    # get extrema
    tmp_extrs = extr(smoothed)
    tmp_extrs['sigma'] = sigma

    # append to general dataframe
    return tmp_extrs


for sigma in tqdm(sigmas):
    tmp_extrs = convert_to_df(doOne(sigma))
    tmp_extrs['sigma'] = sigma
    extrs     = extrs.append(tmp_extrs)


# store the result
HDF = pd.HDFStore(args.out)
HDF['extremas'] = extrs
HDF['dens']     = pd.Panel(dens)
HDF['edges']    = pd.DataFrame(edges.T, columns=['x', 'y', 'z'])
HDF.close()
