from sage.all import *


class RFConfig:
    method_GROEBNER = 1
    method_CRT = 2
    method_VARIETY = 3

    def __init__(self,method = method_GROEBNER,tryTimes = 20,gbLimitNum = 5,fgbThreads = 2) -> None:
        self.method = method
        self.tryTimes = tryTimes
        self.gbLimitNum = gbLimitNum
        self.fgbThreads = fgbThreads

class ETConfig:
    def __init__(self, Ts = []) -> None:
        self.Ts = Ts