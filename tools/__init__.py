import _tools
import pandas as pd
import os
from scipy.io import FortranFile
import numpy as np
import itertools

io = _tools.io
misc = _tools.misc

class Indexable(object):

    def __init__(self, it):
        self.it = iter(it)
        self.already_computed = []

    def __iter__(self):
        for elt in self.it:
            self.already_computed.append(elt)
            yield elt

    def __getitem__(self, index):
        try:
            max_idx = index.stop
        except AttributeError:
            max_idx = index
        n = max_idx - len(self.already_computed) + 1
        if n > 0:
            self.already_computed.extend(itertools.islice(self.it, n))
        return self.already_computed[index]

    def load_all(self):
        self.already_computed.extend(list(self.it))

def check_filename(fun):
    ''' Wrapper that takes the first argument of the function call and asserts
    it is an existing file '''
    def wrapped(filename, *args, **kwargs):
        if not os.path.isfile(filename):
            raise IOError('No such file or directory: \'{}\''.format(filename))
        else:
            return fun(filename, *args, **kwargs)
    return wrapped

@check_filename
def read_particles_wrap(filename):
    '''Parameters:
    -----------
    filename: input string

    Returns:
    --------
    ndim, nparts, nstar, data
    '''
    ndim, nparts, unit = _tools.io.read_particle_header(filename)
    nstar, [x, y, z], [vx, vy, vz], mass, ids, bd = _tools.io.read_particle_data(ndim, nparts, unit)
    data = pd.DataFrame({
        'x': x,
        'y': y,
        'z': z,
        'vx': vx,
        'vy': vy,
        'vz': vz,
        'mass': mass,
        'ids': ids,
        'birth_date': bd
    })
    return ndim, nparts, nstar, data

@check_filename
def read_brick_wrap(filename, dm_type=True, low_mem=False, preload=False):
    ''' Read a brick file.
    Params:
    -------
    filename: string, the brick file
    dm_type: boolean, True if dark matter brick file
    low_mem: boolean, if True, keep minimum data

    Return:
    -------
    res: pandas.DataFrame containing all the information
    '''

    def tmp():
        ff = FortranFile(filename)
        h = {}
        h["nbodies"] = ff.read_ints()
        h["massp"] = ff.read_ints()
        h["aexp"] = ff.read_reals(dtype=np.int32)
        h["omega_t"] = ff.read_reals(dtype=np.int32)
        h["age_univ"] = ff.read_reals(dtype=np.int32)
        h["n_halos"], h["n_subhalos"] = ff.read_ints()

        for i in range(h["n_halos"] + h["n_subhalos"]):
            infos = {
                "header": h
            }
            infos["nparts"] = ff.read_ints()
            infos["members"] = ff.read_ints()
            infos["idh"] = ff.read_ints()
            infos["timestep"] = ff.read_ints()
            infos["mylevel"], infos["hosthalo"], infos["hostsub"], infos["nbsub"], infos["nextsub"] = ff.read_ints()

            infos["mhalo"] = ff.read_reals(dtype=np.int32)
            infos["pos"] = ff.read_reals(dtype=np.int32)
            infos["speed"] = ff.read_reals(dtype=np.int32)
            infos["L"] = ff.read_reals(dtype=np.int32)
            infos["r"], infos["a"], infos["b"], infos["c"] = ff.read_reals(dtype=np.int32)
            infos["ek"], infos["ep"], infos["et"] = ff.read_reals(dtype=np.int32)
            infos["spin"] = ff.read_reals(dtype=np.int32)
            if not dm_type:
                ff.read_reals()
            infos["rvir"],infos["mvir"], infos["tvir"], infos["cvel"] = ff.read_reals(dtype=np.int32)
            ff.read_reals()
            if not dm_type:
                infos["npoints"] = ff.read_ints()
                infos["rdum"] = ff.read_reals(dtype=np.int32)
                infos["density"] = ff.read_reals(dtype=np.int32)

            if low_mem:
                yield {
                    "nparts": infos["nparts"],
                    "members": infos['members']
                }
            else:
                yield infos
        ff.close()
    ind = Indexable(tmp())
    if preload:
        ind.load_all()
    return Indexable(tmp())

@check_filename
def read_info_wrap(filename):
    _tools.io.read_info_headers(filename)

io.read_brick = read_brick_wrap
io.read_particles = read_particles_wrap
io.read_info = read_info_wrap

__all__ = [io, misc]
