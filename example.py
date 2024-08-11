from Crypto.Util.number import getPrime
from AutoCoppersmith.Coppersmith import *
from AutoCoppersmith.Util.Config import ETConfig


from sage.all import *



Nbits = 1024
p,q = getPrime(Nbits // 2),getPrime(Nbits // 2)
N = p * q

bounds = (floor(N ** 0.21),)
roots = tuple(randrange(bound) for bound in bounds)

R = PolynomialRing(Zmod(N),["x"],1)
x = R.gens()[0]
monomials = [x, x ** 2, x ** 3]
f1 = sum(randrange(N) * monomial for monomial in monomials)
f1 -= f1(*roots)


f2 = sum(randrange(N) * monomial for monomial in monomials)
f2 -= f2(*roots)

Cp = Coppersmith(beta = 1,logging_level = Coppersmith.logging_level_DEBUG)
print(Cp.small_roots([f1],bounds,3))
print(Cp.small_roots([f1,f2],bounds,3))


bounds = (floor(N ** 0.1), floor(N ** 0.32))
roots = tuple(randrange(bound) for bound in bounds)
R = PolynomialRing(Zmod(N),["x","y"],2)
x,y = R.gens()
monomials = [x, y, x * y, x ** 2]
f = sum(randrange(N)*monomial for monomial in monomials)
f -= f(*roots)
etcfg = ETConfig()
Cp = Coppersmith(beta = 1, etconfig = etcfg, logging_level = Coppersmith.logging_level_DEBUG)
print(Cp.small_roots([f],bounds,2))


