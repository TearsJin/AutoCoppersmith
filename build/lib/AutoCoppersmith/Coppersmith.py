from sage.all import *

from AutoCoppersmith.Util.Lattice import flatter
from AutoCoppersmith.Util.Findroots import rootsFinder
from AutoCoppersmith.Util.Config import ETConfig
from AutoCoppersmith.Util.loggingConfig import logging
from time import time
import sys

LOG_FORMAT = " %(levelname)s - %(message)s"
sys.set_int_max_str_digits(0)
class Coppersmith:
    '''
        Coppersmith Main Class

        Attributes:
            logging_level_INFO:  Sets the logging level to Info
            logging_level_DEBUG: Sets the logging level to Debug

            mode_MODN:           Sets the equation to be solved in the context of MOD N
            mode_MODP:           Sets the equation to be solved in the context of MOD P
            mode_MODUP:          Sets the equation to be solved in the context of MOD U * P

            code_level_ATTACK:   Sets the information level output during root-finding, outputs less information in ATTACK mode
            code_level_EXP:      Sets the information level output during root-finding, outputs more information in EXP mode
        
        Methods:
            Small_roots:         Main method
    '''
    logging_level_INFO = logging.INFO
    logging_level_DEBUG = logging.DEBUG

    mode_MODN = 1
    mode_MODP = 2
    mode_MODUP = 3

    code_level_ATTACK = 4
    code_level_EXP = 5
    def __init__(self,mode = mode_MODN,modulus = None,beta = 1, rfconfig = None, etconfig = None, ulconfig = None, code_level = code_level_ATTACK,logging_level = logging_level_DEBUG) -> None:
        """
            :param mode:            Sets the context in which the equation is solved
            :param modulus:         Sets the modulus of the equation, required only when mode is set to mode_MODUP
            :param beta:            Sets P as N^beta when mode is set to mode_MODP
            :param rfconfig:        Configuration for root-finding
            :param etconfig:        Configuration for extension strategy
            :param ulconfig:        Configuration for splitting and linearization
            :param code_level:      Code level, used to set the amount of detail obtained during the code execution
            :param logging_level:   Logging level, used to set the amount of logging obtained during the code execution
        """
        self.mode = mode
        self.beta = beta
        self.rfconfig = rfconfig
        self.etconfig = etconfig
        self.ulconfig = ulconfig
        self.code_level = code_level

        logging.basicConfig(level = logging_level, format = LOG_FORMAT )
        if self.beta > 1:
            logging.error("Invalid beta")
            exit()

        if self.mode == self.mode_MODUP:
            if modulus == None:
                logging.error("Mode MODUP need modulus.")
                exit()
            self.modulus = modulus

    def __exgcd(self,mons: list):
        if len(mons) == 2:
            return gcd(*mons)
        return self.__exgcd((gcd(*mons[:2]),) + mons[2:])
        
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

        # Extend strategy

        if self.etconfig != None: 
            if self.etconfig.Ts == []:
                # default extend strategy
                Ts = [0]
                Lm = prod(self.polys).lm()
                lmexps = Lm.exponents()[0]
                lms = sum([a * b for a,b in zip(lmexps,[float(log(b,self.modulus)) for b in self.bounds])])
                for _ in range(1,self.k):
                    Ts.append(ceil((1 - lms) / float(log(self.bounds[_],self.modulus))))
                self.etconfig = ETConfig(Ts = Ts)

        else:
            self.etconfig = ETConfig([0] * self.k)

            
        if len(self.etconfig.Ts) != self.k:
            self.etconfig.Ts = [0] * self.k
            logging.warning("Invalid ETConfig.")
        

        self.M = list(self.M)
        baseM = deepcopy(self.M)
        for v,t in zip(self.vars,self.etconfig.Ts):
            for i in range(floor(self.m * t) + 1):
                self.M = list(set(list(v ** i * vector(baseM)) + self.M))
    
        self.M = set(self.M)

        # Unravelled Linearization
        if self.ulconfig != None:
            qr = self.R.quotient(self.ulconfig.qr)
            M = []
            for mon in self.M:
                M += qr(mon).lift().monomials()
            self.M = set(M)

        logging.debug("ConstructM done.")
        logging.debug("The num of monomials: {}".format(len(self.M)))

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
                    P = prod([poly ** i for i,poly in zip(ii,self.polys)]) * (mon // LM) * (self.u ** (self.m - sum(ii)) * self.modulus ** (max(0,self.t - sum(ii))))
                    if self.ulconfig != None:
                        qr = self.R.quotient(self.ulconfig.qr)
                        P = qr(P).lift()
                    if set(P.monomials()).issubset(set(self.M)) and  P.lm() == mon and sum(ii) > label:
                        candidatePoly = P
                        label = sum(ii)
            self.labels.append(self.t - label)
            self.F.append(candidatePoly)

        logging.debug("ConstructF done.")
        logging.debug("The num of polynomials: {}".format(len(self.F)))

    def __ConstructB(self) -> None:

        try:
            # for sagemath 10.3
            B, monomials = self.F.coefficients_monomials()
        except:
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

    def small_roots(self, fs: list, bounds: list, i: int, u = 1, ROOTS = [], HsFilter = []) -> list:
        """
            Small_roots:
                :param fs:          The original system of equations for which the roots need to be found
                :param bounds:      Upper bounds for all the roots
                :param i:           Sets the size of the Coppersmith lattice, related to the power of the modulus, m = i * n, where n is the number of fs
                :param ROOTS:       Correct solutions, can be provided during experiments to quickly determine if the solutions exist in the lattice after LLL
                :param HsFilter:    Filters the polynomials obtained after LLL, for example, setting it to [0,1,2,3] will select the first 4 polynomials
                :param u:           Needs to be provided when mode is set to mode_MODUP
        """
        self.startTime = time()
        self.f0 = fs[0]
        self.n = len(fs)
        self.R = self.f0.parent()
        self.k = self.R.ngens()
        self.m = i * self.n
        self.t = ceil(self.beta * self.m)
        self.bounds = bounds
        self.u = u

        if self.mode == self.mode_MODN or self.mode == self.mode_MODP:
            if self.u != 1:
                self.u = 1
                logging.warning("Since beta = 1, u will not use.")
            self.modulus = self.f0.base_ring().cardinality()
           
        # monic
        polys = []
        for f in fs:
            LC = f.coefficients().pop(0)
            GCD = self.R(gcd(LC, self.modulus))
            if self.mode == self.mode_MODUP:
                f = (f * inverse_mod(LC // GCD,self.u * self.modulus)) % (self.u * self.modulus)
            else:
                f /= LC // GCD
            polys.append(f.change_ring(ZZ))
            
        self.polys = polys
        self.R = polys[0].parent()
        self.vars = self.R.gens()

        self.__ConstructM()
        self.__ConstructF()
        self.__ConstructB()
        self.Hs = self.B * self.monomials
        
        if self.ulconfig != None:
            Hs = []
            for h in self.Hs:
                Hs.append(h(*self.ulconfig.unqr))
            self.Hs = Hs
            self.R = self.ulconfig.unqr[0].parent()
            self.bounds = self.ulconfig.bounds
        

        if ROOTS != []:
            logging.info("Check roots: ")

            for row in self.Hs:
                row = row.change_ring(ZZ)
                if row(*ROOTS) == 0:
                    print("+",end=" ")
                else:
                    print("x",end=" ")
            print()

        if HsFilter != []:
            self.Hs = [self.Hs[i] for i in HsFilter]

        rf = rootsFinder(self.rfconfig)
        roots = rf.find_roots(self.R,self.Hs,bounds)
        logging.debug("Run time: {}s".format(str(time() - self.startTime)))
        return roots
