# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import pandas as pd
import numpy as np
import os.path
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

try:
    from tqdm import tqdm
except:
    print('Missing tqdm library, install it (fallbacks to dummy function).')
    def tqdm(foo):
        return foo

parser = argparse.ArgumentParser(description='FIXME')
parser.add_argument('--galaxy-list', type=str, required=True)
parser.add_argument('--halo-list', type=str, required=True)
parser.add_argument('--association-list', type=str, required=True)
parser.add_argument('--ramses-output-start', type=str, required=True)
parser.add_argument('--DM-tree-bricks-start', '-dtbs', type=str, required=True)
parser.add_argument('--DM-tree-bricks-end', '-dtbe', type=str, required=True)

# for cpu in range(1, ncpu+1):
#     filename = 'part_{}.out{:0>5}'.format(nsim, cpu)
#     path = os.path.join(args.ramses_output_start, filename)
#     print('{:.2f}% ({}/{})'.format(100.*cpu/ncpu, cpu, ncpu))

#     f = FortranFile(path, 'r')
#     ncpu2 = f.read_ints()[0]
#     ndim2 = f.read_ints()[0]
#     npart = f.read_ints()[0]
#     f.read_ints()
#     nstar = f.read_ints()
#     f.read_ints()
#     f.read_ints()
#     f.read_ints()

#     x = f.read_reals(dtype=np.float64)
#     y = f.read_reals(dtype=np.float64)
#     z = f.read_reals(dtype=np.float64)
#     vx = f.read_reals(dtype=np.float64)
#     vy = f.read_reals(dtype=np.float64)
#     vz = f.read_reals(dtype=np.float64)
#     m = f.read_reals(dtype=np.float64)
#     part_id = f.read_ints()
#     birth = f.read_reals(dtype=np.float32)

#     particles.append(pd.DataFrame({
#         "x": x,
#         "y": y,
#         "z": z,
#         "vx": vx,
#         "vy": vy,
#         "vz": vz,
#         "m": m,
#         "id": part_id,
#         "birth": birth,
#     }, columns=['x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'id', 'birth']))

#     f.close()

def particles_in_halo(tree_brick, start=0, end=None, fun_filter=lambda x: True):
    ''' Open a tree bricks file and associate to each halo the corresponding particles.
    '''
    # Open file
    f = FortranFile(tree_brick, 'r')

    # Give a value to end, by default start + 1
    if end == None:
        end = start + 1

    # Read headers
    nbodies = f.read_ints()[0]
    f.read_reals(dtype=np.float32)
    aexp = f.read_reals(dtype=np.float32)
    f.read_reals(dtype=np.float32)
    age = f.read_reals(dtype=np.float32)
    nhalo, nsubhalo = f.read_ints()
    halo_tot = nhalo + nsubhalo

    halos = {}
    for i in tqdm(range(halo_tot)):
        parts = f.read_ints()[0]
        members = f.read_ints()
        this_id = f.read_ints()[0]
        if (start <= this_id and this_id < end and fun_filter(this_id)):
            for dm_particle_id in members:
                if not halos.has_key(this_id):
                    halos[this_id] = []

                halos[this_id].append(dm_particle_id)
        elif this_id >= end:
            break
        f.read_ints()

        # Irrelevant
        level, hosthalo, hostsub, nbsub, nextsub = f.read_ints()
        mstar = 1e11 * f.read_reals(dtype=np.float32)
        px, py, pz = f.read_reals(dtype=np.float32)
        f.read_reals(dtype=np.float32)
        f.read_reals(dtype=np.float32)
        rad = f.read_reals(dtype=np.float32)[0]
        f.read_reals(dtype=np.float32)
        f.read_reals(dtype=np.float32)
        rvir, mvir, tvir, cvel = f.read_reals(dtype=np.float32)
        f.read_reals(dtype=np.float32)

    f.close()
    return halos

def read_galaxy_list(listfile):
    galFile = FortranFile(listfile, 'r')
    ngal, columns = galFile.read_ints()
    _tmp = (galFile.read_reals(dtype=np.float32)).reshape((columns, ngal)).transpose()
    galaxies = pd.DataFrame(_tmp,
                            columns=['id', 'vt', 'dvz', 'dvr', 'dvtheta', 'mass', 'x', 'y', 'z'])
    galaxies.id.astype(int)

    galaxies['sigma'] = 1/3.*np.sqrt(galaxies.dvz**2 + galaxies.dvtheta**2 + galaxies.dvr**2)
    galaxies['sigmaoverv'] = galaxies.sigma / galaxies.vt
    galaxies['elliptic'] = galaxies.sigmaoverv > 1.5
    galaxies['spiral'] = galaxies.sigmaoverv < 0.8

    return galaxies

def read_halo_list(listfile):
    haloFile = FortranFile(listfile, 'r')
    nhalos, columns = haloFile.read_ints()
    _tmp = (haloFile.read_reals(dtype=np.float32)).reshape((columns, nhalos)).transpose()
    halos = pd.DataFrame(_tmp,
                         columns=['id', 'level', 'mass', 'vx', 'vy', 'vz', 'sigma'])
    halos.id.astype(int)
    halos.level.astype(int)
    return halos

def read_infos(path):
    with open(path, 'r') as f:
        ncpu = int(f.readline().replace('\n', '').split('=')[1])
        ndim = int(f.readline().replace('\n', '').split('=')[1])
        levelmin = int(f.readline().replace('\n', '').split('=')[1])
        levelmax = int(f.readline().replace('\n', '').split('=')[1])
        ngridmax = int(f.readline().replace('\n', '').split('=')[1])
        nstep_coarse = int(f.readline().replace('\n', '').split('=')[1])
        f.readline()
        boxlen = float(f.readline().replace('\n', '').split('=')[1])
        time = float(f.readline().replace('\n', '').split('=')[1])
        aexp= float(f.readline().replace('\n', '').split('=')[1])
        H0 = float(f.readline().replace('\n', '').split('=')[1])
        omega_m = float(f.readline().replace('\n', '').split('=')[1])
        omega_l = float(f.readline().replace('\n', '').split('=')[1])
        omega_k = float(f.readline().replace('\n', '').split('=')[1])
        omega_b = float(f.readline().replace('\n', '').split('=')[1])
        unit_l = float(f.readline().replace('\n', '').split('=')[1])
        unit_d = float(f.readline().replace('\n', '').split('=')[1])
        unit_t = float(f.readline().replace('\n', '').split('=')[1])
        f.readline()
        ordering_type = f.readline().replace('\n', '').split('=')[1].strip()

        headers = f.readline().replace('\n', '').split()

        infos = pd.read_csv(f, names=headers, delim_whitespace=True)

    return infos
def read_association(listfile):
    assocFile = FortranFile(listfile, 'r')
    nassoc, columns = assocFile.read_ints()
    _tmp = (assocFile.read_reals(dtype=np.float32)).reshape((columns, nassoc)).transpose()
    assoc = pd.DataFrame(_tmp,
                         columns=['halo_id', 'level', 'halo_mass', 'gal_id', 'gal_mass'])

    assoc.halo_id.astype(int)
    assoc.level.astype(int)
    assoc.gal_id.astype(int)
    return assoc

if __name__ == '__main__':
    args = parser.parse_args()
    print('Reading lists…')
    ipos = args.ramses_output_start.index('output_')
    nsim = args.ramses_output_start[ipos+7:ipos+13].replace('/', '')

    path = os.path.join(args.ramses_output_start, 'info_' + nsim + '.txt')
    infos = read_infos(path)

    galaxies = read_galaxy_list(args.galaxy_list)
    halo_list = read_halo_list(args.halo_list)
    associations = read_association(args.association_list)

    print('Reading tree…')
    halos_particles = particles_in_halo(args.DM_tree_bricks_start, start=181102)

    print('Finding CPUs containing the particles…')

