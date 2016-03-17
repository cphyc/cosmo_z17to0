from __future__ import print_function

import tools
import numpy as np
import pandas as pd

tests = []
def test(fun):
    def wrap (*args, **kwargs):
        print('Testing "{}":'.format(fun.__doc__), end=' ')
        r = fun(*args, **kwargs)
        r = [] if r is None else r
        r_as_list = list(r)
        if len(r_as_list) == 0:
            print('all good!')
        else:
            print('{} errors'.format(len(r_as_list)))

            for e in r_as_list:
                print('\tE:', e)

    tests.append(wrap)

def run_tests(*args, **kwargs):
    while len(tests) > 0:
        test = tests.pop()

        test(*args, **kwargs)

@test
def Hilbert3D():
    '''Hilbert 3D'''
    ret = tools.misc.hilbert3d([1, 1], [2, 1], [3, 5], 1)
    if ret.shape != (2, ):
        yield 'Shape mismatch'

    for a, b in zip(ret, np.array([ 6.,  5.])):
        if a != b:
            yield '{} != {}'.format(a, b)


@test
def get_cpu_list():
    '''Cpu_list'''
    # Test
    # remove last because it's a missing col
    bound_key = np.array(pd.read_csv('bound_key',
                                     delim_whitespace=True,
                                     header=None)).flatten()[:-1]

    x0, x1 = (np.array([0, 0.1, 0.2]), np.array([0.1, 0.2, 0.3]))
    lvlmax = 17
    ret = tools.misc.get_cpu_list(x0, x1, lvlmax, bound_key)
    cpu_list = (np.array(pd.read_csv('cpu_list',
                                     delim_whitespace=True,
                                     header=None))).flatten()[:-2].astype(int)

    for a, b in zip(cpu_list, ret):
        if (a != b):
            yield [a, b]

@test
def bricks():
    '''Bricks'''
    pass


run_tests()


