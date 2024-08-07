from sage.all import *

import logging
from tqdm import tqdm
from itertools import product
from AutoCoppersmith.Util.Config import RFConfig
from AutoCoppersmith.Util.Misc import check_package_exists

if check_package_exists("fgb_sage"):
    import fgb_sage
    FGB_EXIST = True
else:
    logging.info(msg="Can not find fgb_sage")
    FGB_EXIST = False
    

class rootsFinder:
    """
        Main class for root finding.

        Methods:
            find_roots: Main method for root finding.
    """
    def __init__(self,rfconfig: RFConfig) -> None:
        """
            Initialize rootsFinder with configuration provided by RFConfig.

            :param rfconfig: Configuration for root finding details.
        """
        self.rfconfig = RFConfig() if rfconfig == None else rfconfig           


    def find_roots(self,R,Hs,bounds = []) -> None:
        """
            Main method for root finding.

            :param R:       The ring where roots are located.
            :param Hs:      Polynomials obtained after LLL.
            :param bounds:  Upper bounds of roots.
        """
        method = self.rfconfig.method

        if method == RFConfig.method_GROEBNER:
            roots = self.__find_roots_groebner(R,Hs)
        elif method == RFConfig.method_CRT:
            roots = self.__find_roots_CRT(R,Hs,bounds)
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
            I = Sequence(Hs[:i],R.change_ring(QQ)).ideal()
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


    def __find_roots_CRT(self,R,Hs,bounds):
        roots = []
        vars = R.gens()


        k = len(vars)
        Hs = [h.change_ring(ZZ) for h in Hs]
        maxBound = max(bounds)
        remainders = [[] for _ in range(k)]
        modulus = 1
        moduli = []
        ALERT = True
        
        P = min(Integer(maxBound),Integer(2 ** 25))
        while modulus < maxBound and len(Hs) > 0:
            P = Primes().next(P)
            R = R.change_ring(GF(P))
            solutions = []
            I = R * Hs
            # fgb_sage should add
            I = Ideal(I.groebner_basis())
            if I.dimension() == 0:
                solutions = I.variety()

                if len(solutions) == 1:
                    solutions = [int(solutions[0][var]) for var in vars ]
                    for remainder,solution in zip(remainders,solutions):
                        remainder.append(solution)
                    moduli.append(P)
                    modulus *= P
                    continue

                if len(solutions) > 1 and ALERT:
                    ALERT = False
                    logging.debug("Find more than one roots! The result may goes wrong. ")
            Hs.pop(-1)


        if len(Hs) != 0:
            res = list(int(crt(remainders[i], moduli)) for i in range(k))
            for root in res:
                roots.append([root,root - modulus])
            roots = [root for root in product(*roots)]
            if self.rfconfig.crtRootSign != []:
                i = 0
                while i != len(roots):
                    _ = False
                    for r,sign in zip(roots[i], self.rfconfig.crtRootSign):
                        if r * sign < 0:
                            _ = True
                    if _:
                        roots.pop(i)
                    else:
                        i += 1
        return roots

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
