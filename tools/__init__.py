import _tools
import pandas as pd
import os
io = _tools.io
misc = _tools.misc

def check_filename(fun):
    ''' Wrapper that takes the first argument of the function call and asserts
    it is an existing file '''
    def wrapped(filename, *args, **kwargs):
        if not os.path.isfile(filename):
            raise IOError('No such file or directory: \'{}\''.format(filename))
        else:
            fun(filename, *args, **kwargs)
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
    if not os.isdir(filename):
        raise IOError('No such file or directory: \'{}\''.format(filename))
    ndim, nparts = _tools.io.read_particle_header(filename)
    nstar, x, y, z, vx, vy, vz, mass, ids, bd = _tools.io.read_particle_data(ndim, nparts)
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
def read_brick_wrap(filename, dm_type=True):
    ''' Read a brick file.
    Params:
    -------
    filename: string, the brick file
    dm_type: boolean, True if dark matter brick file

    Return:
    -------
    res: pandas.DataFrame containing all the information
    '''
    if _tools.io.assert_infos() > 0:
        return
    nbodies, aexp, age_univ, nb_of_halos, nb_of_subhalos = _tools.io.read_brick_header(filename)
    nDM = nb_of_halos + nb_of_subhalos
    res = pd.DataFrame(_tools.io.read_brick_data(nDM, dm_type),
                       columns=["mdm", "xdm", "ydm", "zdm", "rvirdm", "mvirdm",
                                "tvirdm", "hlevel", "lxdm", "lydm", "lzdm", "iddm"])
    return res

@check_filename
def read_info_wrap(filename):
    _tools.io.read_info_headers(filename)

io.read_brick = read_brick_wrap
io.read_particles = read_particles_wrap
io.read_info = read_info_wrap

__all__ = [io, misc]
