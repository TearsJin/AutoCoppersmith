```
./AutoCoppersmith
    /Util
        Config.py
        Findroots.py
        Lattice.py
        Misc.py
    Coppersmith.py
```

# AutoCoppersmith

自用的Coppersmith工具，基本的框架就是 Automated Coppersmith的流程。(基本上就是JM策略，没有对线性方程做提升，线性方程可以结合Herrmann的方法做提升，后续看怎样整合进来)

目前可以解决的情况

​	mod N、mod p和mod kp的单变元(多变元)方程(方程组)



## 关于 Config

Config中有许多可以调的东西

```
RFConfig: 用来更改默认的求根方法，现在只有两个，一个gb，一个用传统的variety，CRT的方法还在调bug，后面会弄进来
		  也可以更改求根时的一些细节
ETConfig: JM扩展策略，可以设置Ts来额外的shift指定的变元，例如对于双变元方程f(x,y)，设置Ts = [0,1]会让y多shift一倍的m。
```

## 可以添加的库/环境

安装这些库/环境可以减少运行时间

flatter -> 加速LLL https://github.com/keeganryan/flatter

fgb_sage -> 加速gb https://github.com/mwageringel/fgb_sage



## 后续想要更新的方向

1. 通用有益多项式选取

2. 拆分线性化的适配

3. 扩展策略的默认策略

4. 方程组中方程的比例选择

   

# 安装

```
sage -python -m pip install .
```



# Refference

Solving the Hidden Number Problem for CSIDH and CSURF via Automated Coppersmith.

Approximate Divisor Multiples Factoring with Only a Third of the Secret CRT-Exponents
