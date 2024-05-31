from importlib.util import find_spec

from sage.all import *

def check_package_exists(package_name):
    return False if find_spec(package_name) is None else True