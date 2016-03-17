# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import pandas as pd
import numpy as np
import os.path
from scipy.io import FortranFile
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
from pybloomfilter import BloomFilter
from multiprocessing import Pool
try:
    from tqdm import tqdm
except:
    print('Missing tqdm library, install it (fallbacks to dummy function).')
    def tqdm(foo):
        return foo
from bintrees import FastAVLTree as Tree
parser = argparse.ArgumentParser(description='FIXME')
parser.add_argument('--galaxy-list', type=str, required=True)
parser.add_argument('--halo-list', type=str, required=True)
parser.add_argument('--association-list', type=str, required=True)
parser.add_argument('--ramses-output-start', type=str, required=True)
parser.add_argument('--ramses-output-end', type=str, required=True)
parser.add_argument('--DM-tree-bricks-start', '-dtbs', type=str, required=True)
parser.add_argument('--DM-tree-bricks-end', '-dtbe', type=str, required=True)
parser.add_argument('--process', '-p', type=int, default=1)

def hilbert3D(x, y, z, bit_length, npoint):
    x_bit_mask = np.zeros(bit_length, dtype=bool)
    y_bit_mask = np.zeros(bit_length, dtype=bool)
    z_bit_mask = np.zeros(bit_length, dtype=bool)
    i_bit_mask = np.zeros(3*bit_length, dtype=bool)
    order = np.zeros(npoint)
    if npoint == 1:
        x = [x]
        y = [y]
        z = [z]
    state_diagram = np.array([ 1, 2, 3, 2, 4, 5, 3, 5,
                               0, 1, 3, 2, 7, 6, 4, 5,
                               2, 6, 0, 7, 8, 8, 0, 7,
                               0, 7, 1, 6, 3, 4, 2, 5,
                               0, 9,10, 9, 1, 1,11,11,
                               0, 3, 7, 4, 1, 2, 6, 5,
                               6, 0, 6,11, 9, 0, 9, 8,
                               2, 3, 1, 0, 5, 4, 6, 7,
                              11,11, 0, 7, 5, 9, 0, 7,
                               4, 3, 5, 2, 7, 0, 6, 1,
                               4, 4, 8, 8, 0, 6,10, 6,
                               6, 5, 1, 2, 7, 4, 0, 3,
                               5, 7, 5, 3, 1, 1,11,11,
                               4, 7, 3, 0, 5, 6, 2, 1,
                               6, 1, 6,10, 9, 4, 9,10,
                               6, 7, 5, 4, 1, 0, 2, 3,
                              10, 3, 1, 1,10, 3, 5, 9,
                               2, 5, 3, 4, 1, 6, 0, 7,
                               4, 4, 8, 8, 2, 7, 2, 3,
                               2, 1, 5, 6, 3, 0, 4, 7,
                               7, 2,11, 2, 7, 5, 8, 5,
                               4, 5, 7, 6, 3, 2, 0, 1,
                              10, 3, 2, 6,10, 3, 4, 4,
                               6, 1, 7, 0, 5, 2, 4, 3]).reshape((12, 2, 8))
    for ipt in range(npoint):
        # convert to binary
        for i in range(bit_length):
            bit_mask = 2**i
            x_bit_mask[i] = (x[ipt] & bit_mask) > 0
            y_bit_mask[i] = (y[ipt] & bit_mask) > 0
            z_bit_mask[i] = (z[ipt] & bit_mask) > 0

        # Interleave bits
        for i in range(bit_length):
            i_bit_mask[3*i+2] = x_bit_mask[i]
            i_bit_mask[3*i+1] = y_bit_mask[i]
            i_bit_mask[3*i  ] = z_bit_mask[i]

        cstate = 0
        for i in range(bit_length-1, -1, -1):
            b2 = 1 if i_bit_mask[3*i+2] else 0
            b1 = 1 if i_bit_mask[3*i+1] else 0
            b0 = 1 if i_bit_mask[3*i  ] else 0
            sdigit = b2 * 4 + b1 * 2 + b0
            nstate = state_diagram[cstate, 0, sdigit]
            hdigit = state_diagram[cstate, 1, sdigit]
            i_bit_mask[3*i+2] = hdigit & 0b100 > 0
            i_bit_mask[3*i+1] = hdigit & 0b010 > 0
            i_bit_mask[3*i  ] = hdigit & 0b001 > 0
            cstate = nstate

        order[ipt] = 0
        for i in range(0, 3*bit_length):
            b0 = 1 if i_bit_mask[i] else 0
            order[ipt] = order[ipt] + b0*2**i

    return order

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
                         columns=['id', 'level', 'mass', 'x', 'y', 'z', 'rvir'])
    halos[['id', 'level']] = halos[['id', 'level']].astype(int)

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

    infos.DOMAIN = infos.DOMAIN.astype(np.int32)
    return infos

def read_association(listfile):
    assocFile = FortranFile(listfile, 'r')
    nassoc, columns = assocFile.read_ints()
    _tmp = (assocFile.read_reals(dtype=np.float32)).reshape((columns, nassoc)).transpose()
    assoc = pd.DataFrame(_tmp,
                         columns=['halo_id', 'level', 'halo_mass', 'gal_id', 'gal_mass'])

    assoc[['halo_id', 'level', 'gal_id']] =  assoc[['halo_id', 'level', 'gal_id']].astype(np.int32)
    return assoc

def read_output(path, header_only=True):
    f = FortranFile(path, 'r')
    ncpu = f.read_ints()
    dim = f.read_ints()
    nparts = f.read_ints()
    if header_only:
        f.close()
        return ncpu, dim, nparts
    f.read_ints()
    f.read_ints()
    f.read_ints()
    f.read_ints()
    f.read_ints()

    x = f.read_reals(dtype=np.float64)
    y = f.read_reals(dtype=np.float64)
    z = f.read_reals(dtype=np.float64)

    vx = f.read_reals(dtype=np.float64)
    vy = f.read_reals(dtype=np.float64)
    vz = f.read_reals(dtype=np.float64)

    m = f.read_reals(dtype=np.float64)

    part_ids = f.read_ints()

    birth = f.read_reals(dtype=np.float32)

    f.close()
    return  ncpu, dim, nparts, x, y, z, part_ids

def _process_one(data_file):
    ''' Process one output file to generate a bloom filter'''
    path, dump_name = os.path.split(data_file)
    _, parent_dir = os.path.split(path)

    # ensure the containing folder exists
    bf_dir_path = os.path.join('bloom_filters', parent_dir)
    if not os.path.isdir(bf_dir_path):
        os.mkdir(bf_dir_path)
    bf_file_path = os.path.join(bf_dir_path, dump_name)

    if not os.path.isfile(bf_file_path):
        ncpu, _, nparts, _, _, _, ids = read_output(data_file, header_only=False)
        bf = BloomFilter(nparts, 1./ncpu, bf_file_path)
        bf.update(ids)

    return bf_file_path

def _process_all(data_file_list, use_tqdm=True):
    '''Simple wrapper for a looper'''
    if not use_tqdm:
        iterator = data_file_list
    else:
        iterator = tqdm(data_file_list)

    return [_process_one(data_file) for data_file in iterator]

def build_bloom_filter(basepath):
    ''' Build a bloom filter for each outputs' dump, so that it is very fast to know
    whether a particle is part of a dump.'''

    allfiles = glob.glob(os.path.join(basepath,  'part*.out*'))
    allfiles.sort()

    if args.process > 1:
        p = Pool(args.process)
        # split the list of files into 10-files chunks
        chunk_l = 100
        splitted = [allfiles[i*chunk_l:(i+1)*chunk_l]
                    for i in range(int(len(allfiles) / (chunk_l*1.))+1)]
        return reduce(lambda prev, curr: prev + curr, p.map(_process_all, splitted))
    else:
        return _process_all(allfiles, use_tqdm=True)

def cpu_containing(particles, bloom_filters):
    ''' Iterate over all bloom filter and yield the one containing the particle'''
    for i in tqdm(range(len(bloom_filters))):
        bf = BloomFilter.open(bloom_filters[i])
        for p in particles:
            if p in bf:
                yield i
                break


def cpu_for(X0, X1, lvlmax, cpu_infos, ncpu=4096, ndim=3):
    (x0, y0, z0), (x1, y1, z1) = X0, X1
    idom = np.zeros((2,2,2), dtype=int)
    jdom = np.zeros((2,2,2), dtype=int)
    kdom = np.zeros((2,2,2), dtype=int)

    cpu_min = np.zeros(8, dtype=int)
    cpu_max = np.zeros(8, dtype=int)
    bounding_min = np.zeros(8, dtype=int)
    bounding_max = np.zeros(8, dtype=int)
    dmax = np.max(X1 - X0)

    # TODO: replace with log2
    lmin = int(min(np.ceil(np.log2(1/dmax)), lvlmax))

    bit_length = lmin - 1
    maxdom = 2**bit_length
    imin = imax = jmin = jmax = kmin = kmax = 0
    if bit_length > 0:
        imin = int(x0*maxdom)
        imax = imin + 1
        jmin = int(y0*maxdom)
        jmax = jmin + 1
        kmin = int(z0*maxdom)
        kmax = kmin + 1

    dkey = (1.* 2**(lvlmax + 1) / maxdom)**ndim
    ndom = 8 if bit_length > 0 else 1
    idom[:,:,0] = imin
    idom[:,:,1] = imax
    jdom[:,0,:] = jmin
    jdom[:,1,:] = jmax
    kdom[0,:,:] = kmin
    kdom[1,:,:] = kmax

    idom = idom.reshape(8)
    jdom = jdom.reshape(8)
    kdom = kdom.reshape(8)

    for i in range(ndom):
        if bit_length > 0:
            order_min = hilbert3D(idom[i], jdom[i], kdom[i], bit_length, 1)
        else:
            order_min = 0
        print(order_min)
        bounding_min[i] = order_min*dkey
        bounding_max[i] = (order_min + 1)*dkey

    for i in range(ndom):
        for impi in range(ncpu):
            cpu = cpu_infos.iloc[i]
            if (i > ndom - 5 and impi > ncpu - 5):
                print (cpu)
            if cpu.ind_min <= bounding_min[i] < cpu.ind_max:
                cpu_min[i] = impi
            if cpu.ind_min < bounding_max[i] <= cpu.ind_max:
                cpu_max[i] = impi

    return pd.DataFrame({'cpu_min':cpu_min, 'cpu_max': cpu_max})

if __name__ == '__main__':
    global args
    args = parser.parse_args()
    print('Reading lists…')
    ipos = args.ramses_output_start.index('output_')
    nsim = args.ramses_output_start[ipos+7:ipos+13].replace('/', '')

    path = os.path.join(args.ramses_output_start, 'info_' + nsim + '.txt')
    infos = read_infos(path)

    galaxies = read_galaxy_list(args.galaxy_list)
    halo_list = read_halo_list(args.halo_list)
    associations = read_association(args.association_list)

    # print('Reading tree…')
    # halos_particles = particles_in_halo(args.DM_tree_bricks_start, start=181102)

    print('Loading particles index…')
    bf_end = build_bloom_filter(args.ramses_output_end)
    bf_start = build_bloom_filter(args.ramses_output_start)

    # print('Finding cpus containing halos…')
    # for _id, line in associations[associations.gal_id > 0].iterrows():
    #     if _id > 10:
    #         break
    #     parts = (particles_in_halo(args.DM_tree_bricks_start, start=line.halo_id))[line.halo_id]
    #     cpus = list(cpu_containing(parts, bf_start))
    #     print('{} cpus containing {} particules of halo {}'.format(len(cpus),
    #                                                                len(parts),
    #                                                                line.halo_id))


    # Try to group CPUs


'''%run sort_galaxy.py --galaxy-list lists/list_kingal_00782.dat --halo-list lists/list_halo.dat.bin --association-list lists/associated_halogal_782.dat.bin --ramses-output-start /data52/Horizon-AGN/OUTPUT_DIR/output_00782/ -dtbs /data40b/Horizon-AGN/TREE_DM_raw/tree_bricks752 -dtbe /data40b/Horizon-AGN/TREE_DM_raw/tree_bricks032 --ramses-output-end /data52/Horizon-AGN/OUTPUT_DIR/output_00002/'''
