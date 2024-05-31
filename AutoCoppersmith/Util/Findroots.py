import logging
from AutoCoppersmith.Util.Config import RFConfig
from AutoCoppersmith.Util.Misc import check_package_exists
from tqdm import tqdm

from sage.all import *

if check_package_exists("fgb_sage"):
    import fgb_sage
    FGB_EXIST = False
else:
    FGB_EXIST = False
    

class rootsFinder:
    def __init__(self,rfconfig: RFConfig) -> None:
        self.rfconfig = RFConfig() if rfconfig == None else rfconfig()            

    def find_roots(self,R,Hs,bounds = []) -> None:

        method = self.rfconfig.method

        if method == RFConfig.method_GROEBNER:
            roots = self.__find_roots_groebner(R,Hs)
        elif method == RFConfig.method_CRT:
            pass
        elif method == RFConfig.method_VARIETY:
            roots = self.__find_roots_variety(R,Hs)
        return roots
        
    def __find_roots_groebner(self,R,Hs):
        roots = []
        vars = R.gens()
        if self.rfconfig.gbLimitNum > len(Hs) - 1:
            self.rfconfig.gbLimitNum = -1
            logging.warning("Length of Hs less than gbLimitNum.")
        logging.debug("Start Findroots.")
        for i in tqdm(range(len(Hs) - 1, self.rfconfig.gbLimitNum,- 1)):
            Hs_ = deepcopy(Hs[:i])
            for _ in range(self.rfconfig.tryTimes):
                I = Sequence(Hs_,R.change_ring(QQ)).ideal()
                if FGB_EXIST:
                    I = ideal(fgb_sage.groebner_basis(I,threads = self.rfconfig.fgbThreads, verbosity = 0))
                else:
                    I = ideal(I.groebner_basis())
                if I.dimension() == 0:
                    for root in I.variety(ring = ZZ):
                        root = tuple(root[var] for var in vars)
                        if 0 not in root :
                            roots.append(root)
                if roots != []:
                    return roots
                else:
                    shuffle(Hs_)
                
    def __find_roots_variety(self,R,Hs):
        roots = []
        vars = R.gens()
        H = Sequence([], R.change_ring(QQ))
        for h in filter(None, Hs):
            H.append(h)
            I = H.ideal()
            if I.dimension() == -1:
                H.pop()
            elif I.dimension() == 0:
                for root in I.variety(ring = ZZ):
                    root = tuple(root[var] for var in vars)
                    if 0 not in root:
                        roots.append(root)
        return roots