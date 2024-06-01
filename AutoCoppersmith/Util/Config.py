from sage.all import *

import logging

class RFConfig:
    method_GROEBNER = 1
    method_CRT = 2
    method_VARIETY = 3

    SIGNPOSITIVE = 1
    SIGNNEGATIVE = -1
    def __init__(self,method = method_GROEBNER,gbLimitNum = 5,fgbThreads = 2,crtRootSign = []) -> None:
        self.method = method
        self.gbLimitNum = gbLimitNum
        self.fgbThreads = fgbThreads
        self.crtRootSign = crtRootSign

class ETConfig:
    def __init__(self, Ts = []) -> None:
        self.Ts = Ts

class ULConfig:
    # Unravelled Linearization
    def __init__(self, qr, unqr: list) -> None:
        self.qr = qr
        self.unqr = unqr

        if len(unqr) != qr.parent().ngens():
            logging.error("The length of unqr must equal the number of vars in qr.")
