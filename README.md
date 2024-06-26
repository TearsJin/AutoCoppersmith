```
/AutoCoppersmith
    /Util
        Config.py
        Findroots.py
        Lattice.py
        Misc.py
    Coppersmith.py
```

# AutoCoppersmith

自用的Coppersmith工具，基本的框架就是 Automated Coppersmith的流程。(基本上就是JM策略，没有对线性方程做提升，线性方程可以结合Herrmann的方法做提升，后续看怎样整合进来)


因为一开始这个package只是自己用来做Coppersmith的ctf题和做实验调试用的，如果在使用的过程中体验很差，可以找我做修改，有遇到bug也可以提一下，非常感谢！

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

2. 拆分线性化的适配 (已完成，可以通过Config中的ULConfig设置拆分线性化的等式)

3. 扩展策略的默认策略

4. 方程组中方程的比例选择

   

# 安装

```
sage -python -m pip install .
```



# Refference

Solving the Hidden Number Problem for CSIDH and CSURF via Automated Coppersmith.

Approximate Divisor Multiples Factoring with Only a Third of the Secret CRT-Exponents

# log

2024.5.31   第一次上传

2024.6.1    添加了CRT求根、拆分线性化适配（Config.ULConfig）、fix了mod kp情况下的多项式构造出错的问题和对Coppersmith.py的整体流程进行的调整，使得整个流程适配新的mod kp以及拆分线性化

2024.6.2    给大部分的类和方法写了说明