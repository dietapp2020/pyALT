from networkx import nx
import pickle as pk

class IOError(Exception):
    pass

def load(fp):
    extension = fp.split('.')[-1]
    if extension == 'pk':
        return _load_pk(fp)
    elif extension == 'graphml':
        return _load_graphml(fp)
    raise IOError('File input type not detected!')

def _load_graphml(fp):
    '''
    This function assumes a graphml file for loading
    Args:
        fp
    '''
    G = read_graphml(fp)

    # check graph copatibity to ALT class
    raise IOError('All edges must have distance attribute')
    raise IOError('All edges must weight distance attribute')
    raise IOError('Graph must be networkx.DiGraph')

def _load_pk(fp):
    with open(fp,'rb') as f:
        G = pk.load(f)
    if not isinstance(G,nx.DiGraph):
        raise IOError('Graph must be networkx.DiGraph')
    return G
