import logging
import sys
from AutoCoppersmith.Util.Lattice import flatter
from AutoCoppersmith.Util.Findroots import rootsFinder

from sage.all import *

LOG_FORMAT = " %(levelname)s - %(message)s"
sys.set_int_max_str_digits(0)
class Coppersmith:
    logging_level_INFO = logging.INFO
    logging_level_DEBUG = logging.DEBUG

    code_level_ATTACK = 1
    code_level_EXP = 2
    def __init__(self,beta = 1, rfconfig = None, etconfig = None, code_level = code_level_ATTACK,logging_level = logging_level_DEBUG) -> None:
            
        self.beta = beta
        self.rfconfig = rfconfig
        self.etconfig = etconfig
        self.code_level = code_level

        logging.basicConfig(level = logging_level, format = LOG_FORMAT )

        if self.beta > 1:
            logging.error("Invalid beta")
            exit()

    def __digCalculate(self) -> None:
        B = deepcopy(self.B)
        factors = [monomial(*self.bounds) for monomial in self.monomials]
        for i, factor in enumerate(factors):
            B.rescale_col(i, 1/factor)
        if B.nrows() == B.ncols():
            P,_,_ = B.LU()
            B = B / P.change_ring(QQ)
            dig = {
                "mm": self.m * B.nrows(), 
                "modulus": 0
                }
            Vars = [str(var) for var in self.vars]
            for Var in Vars:
                dig[Var] = 0
            for _ in range(B.nrows()):
                for Var,exp in zip(Vars,list(self.monomials[_].exponents()[0])):
                    dig[Var] += exp

                Mp = B[_,_] 
                while Mp != 1:
                    assert Mp % self.modulus == 0
                    dig["modulus"] += 1
                    Mp = Mp // self.modulus
        self.dig = dig
        logging.debug("The det of B: {}".format(dig))

    def __ConstructM(self) -> None:
        # The origin method of constructing M in AutoCoppersmith
        i = self.m // self.n
        self.M = []

        jjs = []
        for _ in range((i + 1) ** self.n):
            _ = Integer(_).digits(i + 1)
            _ = _ + [0] * (self.n - len(_))
            jjs.append(tuple(_))

        P = lambda jj,polys:prod([poly ** j for j,poly in zip(jj,polys)])
        for jj in jjs:
            self.M += P(jj,self.polys).monomials()

        self.M = set(self.M)

        # extend strategy
        if self.etconfig != None:
            if self.etconfig.Ts == [] or len(self.etconfig.Ts) != self.k:
                self.etconfig.Ts = [0] * self.k
                logging.warning("Invalid ETConfig.")
           

            self.M = list(self.M)
            baseM = deepcopy(self.M)
            for v,t in zip(self.vars,self.etconfig.Ts):
                for i in range(floor(self.m * t) + 1):
                    self.M = list(set(list(v ** i * vector(baseM)) + self.M))
        
            self.M = set(self.M)

        logging.debug("ConstructM done.")
        logging.debug("The num of monomials: {}".format(len(self.M)))
        logging.debug("The monomials: {}".format(self.M))

    def __ConstructF(self) -> None:
        iis = []
        for _ in range((self.m + 1) ** self.n):
            _ = Integer(_).digits(self.m + 1)
            _ = _ + [0] * (self.m - (len(_)))
            if sum(_) <= self.m:
                iis.append(_)

        self.labels = []
        self.F = Sequence([],self.R)
        for mon in self.M:
            candidatePoly = 1
            label = -1
            for ii in iis:
                LM = prod([poly.lm() ** i for i,poly in zip(ii,self.polys)])
                if mon % LM == 0:
                    P = prod([poly ** i for i,poly in zip(ii,self.polys)]) * (mon // LM) * (self.u ** (self.m - sum(ii)) * self.modulus ** (self.t - sum(ii)))
                    if set(P.monomials()).issubset(set(self.M)) and  P.lm() == mon and sum(ii) > label:
                        candidatePoly = P
                        label = sum(ii)
            self.labels.append(self.t - label)
            self.F.append(candidatePoly)

        logging.debug("ConstructF done.")
        logging.debug("The num of polynomials: {}".format(len(self.F)))

    def __ConstructB(self) -> None:
        B, monomials = self.F.coefficient_matrix()
        monomials = vector(monomials)

        # Scale bounds
        factors = [monomial(*self.bounds) for monomial in monomials]
        for i, factor in enumerate(factors):
                B.rescale_col(i, factor)

        B = B.change_ring(QQ)

        self.B = B 
        self.monomials = monomials
        if self.beta == 1 and self.code_level == self.code_level_EXP:
            self.__digCalculate()
        logging.debug("Finial matrix: {}".format(B.parent()))

        logging.debug("Start reduced.")
        B = flatter(B)
        logging.debug("Reduced done.")

        B = B.change_ring(QQ)
        factors = [monomial(*self.bounds) for monomial in monomials]
        for i, factor in enumerate(factors):
            B.rescale_col(i, 1/factor)

        self.B = B 
        self.monomials = monomials
        logging.debug("ConstructB done.")

    def small_roots(self,fs: list,bounds: list,i: int,u = 1,ROOTS = [],HsFilter = []) -> list:
        """
        Small_roots:
        :param fs: Origin polys
        :param bounds: Bounds of roots
        :param i: m = i * n, where m is power of modulus, n is the number of origin polys
        :param ROOTS: The correct roots
        :param HsFilter: Select Hs
        :param u: when beta < 1, mod u * N ^ beta
        """
        self.f0 = fs[0]
        self.n = len(fs)
        self.R = self.f0.parent()
        self.k = self.R.ngens()
        self.modulus = self.f0.base_ring().cardinality()
        self.m = i * self.n
        self.t = self.beta * self.m
        self.bounds = bounds
        self.u = u

        if self.beta == 1 and self.u != 1:
            self.u = 1
            logging.warning("beta = 1 then u will not use")


        # monic
        polys = []
        for f in fs:
            LC = f.coefficients().pop(0)
            GCD = gcd(LC, self.modulus)
            f /= LC // GCD
            polys.append(f.change_ring(ZZ))

        self.R = polys[0].parent()
        self.polys = polys
        self.vars = self.R.gens()

        self.__ConstructM()
        self.__ConstructF()
        self.__ConstructB()
        self.Hs = self.B * self.monomials


        if ROOTS != []:
            logging.info("Check roots: ")
            for row in self.Hs:
                row = row.change_ring(QQ)
                if row(*ROOTS) == 0:
                    print("+",end=" ")
                else:
                    print("x",end=" ")
            print()

        if HsFilter != []:
            self.Hs = [self.Hs[i] for i in HsFilter]


        rf = rootsFinder(self.rfconfig)
        roots = rf.find_roots(self.R,self.Hs,bounds)
        return roots
