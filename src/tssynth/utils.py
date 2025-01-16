import tempfile, os
from .solar_abundances import periodic_table

def mkdtemp():
    """
    Create a temporary directory.
    This function can be modified if you want to use a specific path,
    e.g. for debugging purposes or on a cluster scratch space.
    """
    
    return tempfile.mkdtemp(dir=os.environ['TWD_BASE'])

def parse_XFe_dict(XFedict):
    """
    """
    keys = list(XFedict.keys())
    values = list(XFedict.values())
    new_keys = []
    new_keys = [key if isinstance(key, int) else element_to_atomic_number(key) for key in keys]
    return dict(zip(new_keys, values))

def element_to_atomic_number(element):
    """
    Convert an element symbol to its atomic number.
    """
    try:
        return periodic_table.index(element)
    except ValueError as e:
        print(e)
        print(f"Element {element} not found in periodic table, returning None.")
        return None

