import matplotlib.pyplot as plt
import argparse
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from pint import UnitRegistry
ur = UnitRegistry()

parser = argparse.ArgumentParser(description='Compute the smoothing tree of a halo.')
parser.add_argument('--in', dest='infile', type=str, help='Path to the output file of compute_extrema (in hdf format)', required=True)
parser.add_argument('--out', '-o', type=str, help='Prefix of the outputs')
parser.add_argument('--infofile', type=str, help='Path to the information file (the one containing units of RAMSES, …)')

args = parser.parse_args()

# read the info file
infos = dict()
infos['headers'] = pd.read_csv(args.infofile, sep=' *= *', nrows=19, names=['key', 'value'], index_col='key').T
infos['domain'] = pd.read_csv(args.infofile, delim_whitespace=True, skiprows=20)


# read the center
from scipy.io import FortranFile
ff = FortranFile('data/halo_536-centers.bin')
ff.read_ints() # don't care
outputs = ff.read_ints()
centers = ff.read_reals().reshape(len(outputs), 3)
mins    = ff.read_reals().reshape(len(outputs), 3)
span    = ff.read_reals().reshape(len(outputs), 3)
maxs    = ff.read_reals().reshape(len(outputs), 3)

# create the output dir if required
# if not os.path.isdir(args.out):
#     os.mkdir(args.out)


HDF = pd.HDFStore(args.infile)
# read the data
df    = HDF['extremas']
dens  = HDF['dens']
edges = HDF['edges']

def to(v, unit):
    return infos['headers']['unit_l'].value * np.array(v) * (1*ur.cm).to(unit)

nbin = len(edges) - 1

# estimate positions
def realPos(S, edges):
    ''' Convert the sery into real positions in cm'''
    m, M = edges.min(), edges.max()
    return to(((S-1) / (nbin - 1) * M + m), 'Mpc')

# compute 'curvature' (product of eigen values of the hessian)
df['curvature'] = np.abs(df.evx * df.evy * df.evz)
df['rx'] = realPos(df.x, edges.x)
df['ry'] = realPos(df.y, edges.y)
df['rz'] = realPos(df.z, edges.z)

# some values have x,y,z above maximum, remove them
posMask = ((df.x >= 0) & (df.x < nbin) &
           (df.y >= 0) & (df.y < nbin) &
           (df.z >= 0) & (df.z < nbin))

# plots configuration
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.style.use(['ggplot'])
import matplotlib.pyplot as plt


fig = plt.figure()
############## plot sigma vs x
# plt.scatter('x', 'sigma', data=df[df.type == 4], c='r', marker='o', label='Maximum')
# plt.scatter('x', 'sigma', data=df[df.type == 3], c='b', marker='x', alpha=0.5, label='Point selle, 2 maximums')
# plt.xlim(0, len(edges))

# plt.legend(loc='best')
# plt.xlabel('x')
# plt.ylabel('$\sigma$')

############## plot sigma vs dist to center
# plt.scatter('dist', 'sigma', data=df[df.type == 4], c='r', marker='o', label='Maximum')
# plt.scatter('dist', 'sigma', data=df[df.type == 3], c='b', marker='.', alpha=0.5, label='Point selle, 2 maximums')

############## plot sigma vs x,y
ax = fig.add_subplot(111, projection='3d')

t = 4
# map evz onto 20-120
size_evz = 50*(df[df.type==t].evz.max()/df[df.type==t].evz) + 10
m = posMask & (df.type == t) & (df.sigma >= 5) & (df.i < nbin) & (df.j < nbin) & (df.k < nbin)
x, y = df[m].rx, df[m].ry
ax.scatter3D(x, y, df[m].sigma,
             s=size_evz,
             c='r', marker='o', label='Maximum')

t = 3
m = posMask & (df.type == t) & (df.sigma >= 5) & (df.i < nbin) & (df.j < nbin) & (df.k < nbin)
x, y = df[m].rx, df[m].ry
ax.scatter3D(x, y, df[m].sigma,
             c='b', marker='.', alpha=0.5, label='Point selle, 2 maximums')

t = 2
m = posMask & (df.type == t) & (df.sigma >= 5) & (df.i < nbin) & (df.j < nbin) & (df.k < nbin)
x, y = df[m].rx, df[m].ry
ax.scatter3D(x, y, df[m].sigma,
             c='g', marker='.', alpha=0.1, label='Point selle, 1 maximums')

# t = 1
# m = posMask & (df.type == t) & (df.sigma >= 5) & (df.i < nbin) & (df.j < nbin) & (df.k < nbin)
# x, y = df[m].rx, df[m].ry
# ax.scatter3D(x, y, df[m].sigma,
#              c='grey', marker='+', alpha=1, label='Minimum')

# ax.set_xlim(to(edges.x.min(), 'Mpc').m, to(edges.x.max(), 'Mpc').m)
# ax.set_ylim(to(edges.y.min(), 'Mpc').m, to(edges.y.max(), 'Mpc').m)
ax.set_xlabel('x [Mpc/ h]')
ax.set_ylabel('y [Mpc/ h]')
ax.set_zlabel('$\sigma$')

cMpc = to(centers, 'Mpc').m
ax.plot([cMpc[-1, 0], cMpc[0, 0]], [cMpc[-1, 1], cMpc[0, 1]], [df.sigma.min(), df.sigma.max()],
        label='Halo', c='black')

ax.legend()
plt.tight_layout()
# for azim in tqdm(np.linspace(0, 360)):
#     ax.view_init(elev=0, azim=azim)
#     ax.legend()

#     # plt.savefig('/tmp/plots/%s.pdf' % azim)
#     plt.savefig('/tmp/plots/%s.png' % azim, dpi=180)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter3D(df[m].rx, df[m].ry, df[m].rz, c=df[m].sigma, s=5)
cb = fig.colorbar(sc)

ax.set_xlabel('x [Mpc/h]')
ax.set_ylabel('y [Mpc/h]')
ax.set_zlabel('z [Mpc/h]')

ax.set_xlim()
cb.set_label('$\sigma$')
fig.tight_layout()

plt.show()
