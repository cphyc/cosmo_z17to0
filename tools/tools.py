import _tools
import os
from scipy.io import FortranFile
import numpy as np

io = _tools.io
misc = _tools.misc

class Indexable(object):
    def __init__(self, it):
        self.it = it
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
            nexts = [next(self.it) for _ in range(n)]
            self.already_computed.extend(nexts)
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
    wrapped.__doc__ = fun.__doc__
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
def read_brick(filename, dm_type=True, low_mem=False, preload=False, tqdm=lambda e: e):
    ''' Read a brick file.
    Params:
    -------
    filename: string, the brick file
    dm_type: boolean, True if dark matter brick file
    low_mem: don't save all data. If a list of keys is given, yield the keys for each halo,
    else, yield only the members and the number of particles

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

        for i in tqdm(range(h["n_halos"] + h["n_subhalos"])):
            infos = {
                "header": h
            }
            infos["nparts"] = ff.read_ints()
            infos["members"] = ff.read_ints()
            infos["idh"] = ff.read_ints()
            infos["timestep"] = ff.read_ints()
            infos["hlevel"], infos["hosthalo"], infos["hostsub"], infos["nbsub"], infos["nextsub"] = ff.read_ints()

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

            if low_mem != None:
                try:
                    keys = list(low_mem)
                except:
                    keys = ['nparts', 'members']

                tmp = {}
                for key in keys:
                    try:
                        tmp[key] = infos[key]
                    except KeyError:
                        print('Invalid key {}, can be any of', infos['keys'])
                yield tmp
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
