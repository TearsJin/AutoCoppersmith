from sage.all import sqrt, Matrix, Integer, Zmod, PolynomialRing, inverse_mod
from Crypto.Util.number import getPrime
from AutoCoppersmith.Coppersmith import Coppersmith

Nbits = 200
beta = 0.5
p, q = getPrime(int(Nbits * beta)), getPrime(Nbits - int(Nbits * beta))
N = p * q
s = - (p + q)

A = N + 1
R = int(sqrt(N))

delta = 0.27
d = getPrime(int(delta * Nbits))
e = int(inverse_mod(d, (p - 1) * (q - 1)))
k = (e * d - 1) // (p - 1) // (q - 1)

M = Matrix([[R, e], [0, -A]])
M = M.LLL()
l1, l2 = Integer(M[0][0] / R), M[1][0] / R
print(Integer(d).nbits() - l1.nbits())
for x in range(2 ** 20):
    if (d - x * l1) % l2 == 0:
        z0 = (d - x * l1) // l2
        print(x, z0)
        break
    if (d + x * l1) % l2 == 0:
        z0 = (d + x * l1) // l2
        print(-x, z0)
        break

PR = PolynomialRing(Zmod(e * l1), ["x", "y", "z"], 3)
x, y, z = PR.gens()
f = x + int(A) * y - int(e * l2) * z

Cp = Coppersmith(beta=1, logging_level=Coppersmith.logging_level_DEBUG)
roots = Cp.small_roots(
        [f],
        [int(N ** (1/2 + delta)), int(N ** delta), int(N ** (delta - 1/4))],
        3,
        ROOTS=[k * s + 1, k, z0]
        )
print(roots, k, z0)
