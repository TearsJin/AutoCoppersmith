# from Crypto.Util.number import getPrime
# from sage.all import *
# from AutoCoppersmith.Coppersmith import Coppersmith
# from AutoCoppersmith.Util.Config import RFConfig

# p = getPrime(256)
# q = getPrime(256)

# N = p * q

# a0 = getPrime(64)
# a1 = getPrime(512)
# k = 0x1111111111111111

# X = 63
# Y = 63
# x0,y0 = getPrime(X),getPrime(Y)
# X,Y = 2 ** X,2 ** Y

# PR = PolynomialRing(Zmod(k * N),["x","y"],2)
# x,y = PR.gens()

# f = x ** 2 + a0 * x - y + a1 - (x0 ** 2 + a0 * x0 - y0 + a1 - k * p)

# rfcfg = RFConfig(method = RFConfig.method_CRT)
# Cp = Coppersmith(Coppersmith.mode_MODP,beta = float((256 + 64) / (512 + 64)),rfconfig = rfcfg)
# print(Cp.small_roots([f],[X,Y],i = 13,ROOTS = [x0,y0]))

# PR = PolynomialRing(ZZ,["x","y"],2)
# x,y = PR.gens()

# f = x ** 2 + a0 * x - y + a1 - (x0 ** 2 + a0 * x0 - y0 + a1 - k * p)
# rfcfg = RFConfig(method = RFConfig.method_CRT)
# Cp = Coppersmith(Coppersmith.mode_MODP,beta = 0.5,rfconfig = rfcfg)
# print(Cp.small_roots([f],[X,Y],i = 13,u = k,ROOTS = [x0,y0]))


from AutoCoppersmith.Coppersmith import Coppersmith
from AutoCoppersmith.Util.Config import RFConfig

from sage.all import *


N = 16560864354170001754058994512672943318940563370976087798466444727881483330942761871352258869577528553572663755957451586119530166974019489639824958689951617276161495644698610967505623963141115301453407833169441870214457028613126422544192376164562636961183685768882768312923714718974394046529737415829239869985138300035260197768109884692196300418128283001411459636932202117661934374282660653623201000258881545142769739613476895287578432668782523878551310276636940681367529703509645615758511763411234629347145718846540646434928674409510317321022775260360951885592032817087179680011170316690231298208905101051835044960739
k = 1384368362335504207966549759103327813579924965137177857370654240988463990646382435845333383
X = 2 ** 790
rfcfg = RFConfig(method = RFConfig.method_GROEBNER,crtRootSign = [RFConfig.SIGNPOSITIVE])
PR = PolynomialRing(ZZ,"x",1)
x = PR.gens()[0]
f = 0x8150ef932c24cabd * x  + 156836185692941550685700830474694810801481665972206761284102272745592205857267403669744401464803337850596191972703703054530778045202934271320509563462849102899119175874347103045300490340884691560845949337977390221147286384470287571194570435220134004423720066330220713360755490696975534869765259187925645069511484082883584773826540218677265291260254674449374287218468318302260595563533738636720497174
Cp = Coppersmith(mode = Coppersmith.mode_MODUP,modulus = N,beta = 0.5,rfconfig = rfcfg)
print(Cp.small_roots(fs = [f],bounds = [X],i = 25,u = k))


# Nbits = 2048
# R1 = 300
# R2 = 790
# win = 0
# beta = float((R1 + Nbits // 2) / (R1 + Nbits))
# rfcfg = RFConfig(method = RFConfig.method_GROEBNER,crtRootSign = [RFConfig.SIGNPOSITIVE])
# print(beta,(Nbits + R1) * beta ** 2,Nbits // 4 + R1)

# x0 = getPrime(R2)
# a0 = getPrime(2048)

# for _ in tqdm(range(1)):
#     p,q = getPrime(Nbits // 2),getPrime(Nbits // 2)
#     N = p * q
#     k = getPrime(R1)
#     kp = k * p

#     PR = PolynomialRing(Zmod(k * N),"x",1)
#     x = PR.gens()[0]
#     f = 0x8150ef932c24cabd * x + a0 - (0x8150ef932c24cabd * x0 + a0 - kp)

#     print(Integer(f(x0)) % kp)
#     Cp = Coppersmith(mode = Coppersmith.mode_MODP,beta = beta,rfconfig = rfcfg)
#     ROOTS = Cp.small_roots(fs = [f],bounds = [2 ** R2],i = 25,ROOTS = [x0])
#     print(ROOTS)

#     PR = PolynomialRing(ZZ,"x",1)
#     x = PR.gens()[0]
#     f = 0x8150ef932c24cabd * x + a0 - (0x8150ef932c24cabd * x0 + a0 - kp)
#     Cp = Coppersmith(mode = Coppersmith.mode_MODUP,modulus = N,beta = 0.5,rfconfig = rfcfg)
#     print(f,N,k)
#     # f = (x1 + x) ** 2 + a0 * (x1 + x) - (x0 ** 2 + a0 * x0 - kp)
#     print(Integer(f(x0)) % kp)
#     ROOTS = Cp.small_roots(fs = [f],bounds = [2 ** R2],i = 25,u = k,ROOTS = [x0])
#     print(ROOTS,x0)
#     if ROOTS != []:
#         win += 1
# print(win)


