from sage.all import *

import logging

class RFConfig:
    """
        Class for setting root-finding details

        Attributes:
            method_Groebner:    Used to set the root-finding method to Groebner
            method_CRT:         Used to set the root-finding method to CRT
            method_VARIETY:     Used to set the root-finding method to variety

            SIGNPOSITIVE:       Used for the sign of the final root in CRT root-finding, positive
            SIGNNEGATIVE:       Used for the sign of the final root in CRT root-finding, negative
    """

    method_GROEBNER = 1
    method_CRT = 2
    method_VARIETY = 3

    SIGNPOSITIVE = 1
    SIGNNEGATIVE = -1
    def __init__(self,method = method_GROEBNER,gbLimitNum = 5,fgbThreads = 2,crtRootSign = []) -> None:
        """
            Initialization of RFConfig

            :param method:      Sets the root-finding method, options are method_GROEBNER, method_CRT, method_VARIETY
            :param gbLimitNum:  Sets the minimum number of polynomials used in GROEBNER root-finding
            :param fgbThreads:  Sets the number of threads for fgb_sage root-finding
            :param crtRootSign: Sets the sign of the final root in CRT root-finding; if not set, all cases will be output
        """

        self.method = method
        self.gbLimitNum = gbLimitNum
        self.fgbThreads = fgbThreads
        self.crtRootSign = crtRootSign

class ETConfig:
    """
        Class for setting details of extension strategy.
        
        Attributes:
            Ts: The additional shifts in the extension strategy, for example, when there are three variables, setting [0, 0.2, 0.7] will shift the second variable by 0.2 times, and the third variable by 0.7 * m times.
    """
    def __init__(self, Ts = []) -> None:
        self.Ts = Ts

class ULConfig:
    """
        Class for setting details of Unravelled Linearization.
    """
    def __init__(self, qr, unqr: list, bounds: list) -> None:
        """
            Initialize the Unravelled Linearization

                :param qr:      Equations between new variables in the Unravelled Linearization.
                :param unqr:    Relationships between old variables and new variables.
                :param bounds:  The bounds of old variables

            Example:
                PR1.<u1,u2,u3> = PolynomialRing(ZZ,3)
                PR2.<x,y> = PolynomialRing(ZZ,2)
                qr = u2 ** 2 - u1 - u3            #  u2 ** 2 = u1 + u3        
                unqr = [x ** 2 - y, x , y]          #  u1 = x ** 2 - y, u2 = x, u3 = y 
                ulconfig = ULConfig(qr,unqr)
        """
        self.qr = qr
        self.unqr = unqr
        self.bounds = bounds

        if len(unqr) != qr.parent().ngens():
            logging.error("The length of unqr must equal the number of vars in qr.")
