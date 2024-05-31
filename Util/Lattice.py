from re import findall
from subprocess import check_output
import logging

from sage.all import *

def flatter(L: Matrix) -> Matrix:
    # compile GitHub - keeganryan/flatter: Fast lattice reduction and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in L) + "]]"
    try:
        ret = check_output(["flatter"], input=z.encode())
        return matrix(L.nrows(), L.ncols(), map(int, findall(b"-?\\d+", ret)))
    except:
        logging.info("Flatter not found, use FpLLL.")
        return Matrix(L).LLL()
    
