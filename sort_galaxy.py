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

from mpl_toolkits.mplot3d import Axes3D
import tools

parser = argparse.ArgumentParser(description='FIXME')
parser.add_argument('--galaxy-list', type=str,
                    default='lists/list_kingal_00782.dat')
parser.add_argument('--halo-list', type=str,
                    default='lists/list_halo.dat.bin')
parser.add_argument('--association-list', type=str,
                    default='lists/associated_halogal_782.dat.bin')
parser.add_argument('--ramses-output-start', type=str,
                    default='/data52/Horizon-AGN/OUTPUT_DIR/output_00032/')
parser.add_argument('--ramses-output-end', type=str,
                    default='/data52/Horizon-AGN/OUTPUT_DIR/output_00782/')
parser.add_argument('--DM-tree-bricks', '-dtb', type=str,
                    default='/data40b/Horizon-AGN/TREE_DM_raw/tree_bricks752')
parser.add_argument('--process', '-p', type=int, default=1)

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
    print(listfile)
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

    metadata = {
        "ncpu": ncpu,
        "ndim": ndim,
        "levelmin": levelmin,
        "levelmax": levelmax,
        "time": time,
        "aexp": aexp,
        "unit_l": unit_l,
        "unit_d": unit_d,
        "unit_t": unit_t
    }
    infos.metadata = metadata
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
    return ncpu, dim, nparts, x, y, z, part_ids

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

def cpu_containing(particles, bloom_filters, yieldAll=True):
    ''' Iterate over all bloom filter and yield the one containing the particle'''
    for cpu in tqdm(range(len(bloom_filters))):
        bf = BloomFilter.open(bloom_filters[cpu])
        yieldCPU = False
        cpu_contains = []
        for p in particles:
            if p in bf:
                yieldCPU = True
                cpu_contains.append(p)
        if yieldAll or yieldCPU:
            yield cpu+1, cpu_contains

def run_all_galaxies(galaxies, halo_list, associations, information, bf, ramses_output_dir):
    _tmp = associations.gal_id[associations.gal_id > 0]
    n_galaxies = _tmp.size
    _galaxies_treated = 0
    for _, gal_id in _tmp.iteritems():
        _galaxies_treated += 1

        print(('Get information of all particles of halo associated to galaxy {}'.format(gal_id) +
               ' ({}/{} - {:.2f}%)'.format(_galaxies_treated, n_galaxies,
                                           100.*_galaxies_treated/n_galaxies)))

        yield v0(gal_id, halo_list, associations, information, bf, ramses_output_dir)

def v0(gal_id, halo_list, associations, information, bf, ramses_output_dir):
    ramses_dump = int(ramses_output_dir.split('/')[-2].split('_')[-1])

    gal_halo = associations[associations.gal_id == gal_id]
    halo_id = int(gal_halo.halo_id)
    particles = pd.DataFrame({
        'id': particles_in_halo(args.DM_tree_bricks, start=halo_id)[halo_id]
    })

    particles_as_arr = particles.as_matrix().flatten()
    cpus = list(cpu_containing(particles_as_arr, bf))
    data = pd.DataFrame(columns=['x', 'y', 'z', 'ids'])
    cpu_visited = {}

    # create a set of found particles
    particles_found = set()
    particles_as_set = set(particles.id)

    # routine to filter and fill particles_found
    def filter_out(x, y, z, ids):
        # list of particles we're looking for
        particles_list = particles_as_set - particles_found

        for i in range(len(ids)):
            if ids[i] in particles_list:
                particles_found.add(ids[i])
                yield x[i], y[i], z[i], ids[i]

    # Sort CPUs by number of particles in them (first has most particles)
    def sortFun(x, y):
        lx, ly = len(x[1]), len(y[1])
        if (lx < ly):
            return 1
        elif (lx > ly):
            return -1
        else:
            return 0
    cpus.sort(sortFun)
    print('Read {} cpus containing the {} particles'.format(len(cpus), len(particles)))
    for cpu, part_in_cpu_raw in cpus:
        # remove already visited particles from part_in_cpu
        # part_in_cpu = set(part_in_cpu_raw) - particles_found
        part_remaining = particles_as_set - particles_found
        if len(part_remaining) == 0:
            print('\t\tFound everything at cpu {}!'.format(cpu))
            break

        print('\t~{} parts in cpu {} ({} remaining)'.format(len(part_in_cpu_raw),
                                                          cpu, len(part_remaining)))

        ncpu, dim, nparts, x, y, z, part_ids = read_output(
            os.path.join(ramses_output_dir,
                         'part_{:0>5}.out{:0>5}'.format(ramses_dump, cpu)), header_only=False)
        _tmp = list(filter_out(x, y, z, part_ids))
        # no elements found, go to next process
        if len(_tmp) == 0:
            continue
        arr = np.array(_tmp)

        data_tmp = pd.DataFrame(arr, columns=['x', 'y', 'z', 'ids'])
        data = data.append(data_tmp)

    # data_halo = data[[_d.ids in particles_as_set for k, _d in tqdm(data.iterrows())]]
    return gal_id, halo_id, data


if __name__ == '__main__':
    global args
    args = parser.parse_args()
    print('Reading lists…')
    def get_ramses_output(ramses_output_path):
        # get the part in "output_00xxx"
        output = os.path.dirname(ramses_output_path).split('/')[-1]
        num = output.split('_')[-1]
        return num
    nsim = get_ramses_output(args.ramses_output_start)
    path = os.path.join(args.ramses_output_start, 'info_' + nsim + '.txt')
    infos_start = read_infos(path)

    nsim = get_ramses_output(args.ramses_output_end)
    print(nsim)
    path = os.path.join(args.ramses_output_end, 'info_' + nsim + '.txt')
    infos_end = read_infos(path)

    galaxies = read_galaxy_list(args.galaxy_list)
    _tmp_rows, _tmp_cols = tools.io.read_list_header(args.halo_list)
    halo_list = tools.io.read_list_data(_tmp_rows, _tmp_cols)
    associations = read_association(args.association_list)

    # print('Reading tree…')
    # halos_particles = particles_in_halo(args.DM_tree_bricks_start, start=181102)

    print('Loading particles index…')
    bf_end = build_bloom_filter(args.ramses_output_end)
    bf_start = build_bloom_filter(args.ramses_output_start)

    genstart = run_all_galaxies(galaxies, halo_list, associations, infos_start, bf_start,
                  args.ramses_output_start)
    genend = run_all_galaxies(galaxies, halo_list, associations, infos_end, bf_end,
                args.ramses_output_end)

    def gen_plot(parts_start, parts_end, gal_id, halo_id):
        fig = plt.figure(figsize=(16, 9))
        fig.set
        ax = fig.add_subplot(121, projection='3d')
        ax.set_title('Start')
        ax.scatter3D(parts_start.x, parts_start.y, parts_start.z)
        ax = fig.add_subplot(122, projection='3d')
        ax.set_title('End')
        ax.scatter3D(parts_end.x, parts_end.y, parts_end.z)

        fig.suptitle('Halo {}, galaxy {}'.format(halo_id, gal_id))
        plt.show()

    def loop():
        gal_id, halo_id, parts_start = next(genstart)
        gal_id, halo_id, parts_end = next(genend)

        gen_plot(parts_start, parts_end, gal_id, halo_id)

    def plot(gal_id):
        _, _, parts_start = v0(gal_id, halo_list, associations, infos_start, bf_start,
                  args.ramses_output_start)
        _, _, parts_end = v0(gal_id, halo_list, associations, infos_end, bf_end,
                args.ramses_output_end)
        halo_halo = int(associations[associations.gal_id == gal_id].halo_id)
        gen_plot(parts_start, parts_end, gal_id, halo_id)
