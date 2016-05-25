import pandas as pd


def read_infos(path):
    return dict(
        infos=pd.read_csv(path, sep=' *= +', nrows=17,
                          names=['key', 'value'], index_col='key',
                          engine='python'),
        domain=pd.read_csv(path, delim_whitespace=True, skiprows=20, index_col='DOMAIN')
        )
