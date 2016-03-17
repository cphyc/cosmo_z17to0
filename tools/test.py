from tools import tools
import numpy as np
import pandas as pd
# Test Hilbert 3D
ret = tools.hilbert3d([1, 1], [2, 1], [3, 5], 1)
assert(ret.shape == (2, ))

# Test
# remove last because it's a missing col
bound_key = np.array(pd.read_csv('bound_key',
                                 delim_whitespace=True,
                                 header=None)).flatten()[:-1]
x0, x1 = (np.array([0, 0.1, 0.2]), np.array([0.1, 0.2, 0.3]))
lvlmax = 17
print(x0, x1, lvlmax, bound_key.shape)
ret = tools.get_cpu_list(x0, x1, lvlmax, bound_key)
