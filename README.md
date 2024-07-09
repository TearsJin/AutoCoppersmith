# AutoCoppersmith

Coppersmith工具，流程参考的是Automated Coppersmith的流程。(基本上就是JM策略，没有对线性方程做提升，线性方程可以结合Herrmann的方法做提升，后续看怎样整合进来)

这个package只是用来做调试用的Coppersmith的ctf题和做实验，如果在使用的过程中体验很差，可以找我做修改，有遇到bug也可以提一下，非常感谢！

目前可以解决的情况

​	mod N、mod p和mod kp的单变元(多变元)方程(方程组)

```
/AutoCoppersmith
    /Util
        Config.py
        Findroots.py
        Lattice.py
        Misc.py
    Coppersmith.py
```

## Coppersmith.py
里面主要是一个`Coppersmith`类


```python
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
```
* `mode`: 设置模方程的类型，现在有三种类型可以设置
    * `mode_MODN`:      设置方程是模N的。
    * `mode_MODP`：     设置方程是模N的一个未知的因子p的。
    * `mode_MODUP`：    设置方程是模N的一个未知因子的倍数UP的，且U已知。
* `modulus`: 设置方程的模数，只有将模式设置为`mode_MODUP`才需要使用。当模式设置为`mode_MODUP`时，定义的方程必须是在`ZZ`下的，例如
    ```
        PR = PolynomialRing(ZZ,"x",1)
        x = PR.gens()[0]
        f = .......
    ```
    然后在后续初始化`Coppersmith`类的时候把模数填上即可。其他的两个模式都只需要把变量定义在对应的环上面就行了，如
    ```
        PR = PolynomialRing(Zmod(N),"x",1)
        x = PR.gens()[0]
        f = .......
    ```
* `beta`: 在`mode_MODP`和`mode_MODUP`两个模式下需要给定P的大小。
* `rfconfig`(optional): 求根类`class rootsFinder`使用的配置，如果不设置则用默认配置。细节见后文"Config.py""。
* `etconfig`(optional): 扩展策略使用的配置，如果不设置则用默认配置。细节见后文"Config.py"。
* `ulconfig`(optional): 拆分线性化的配置，如果不设置则不使用拆分线性化。细节见后文"Config.py""
* `code_level`(optional): 设置代码等级。可以设置成两类
    * `code_level_ATTACK`(default): 设置为`code_level_ATTACK`时，不计算格的理论界。
    * `code_level_EXP`: 设置为`code_level_EXP`时，计算格的理论界。
* `logging_level`(optional): 设置日志信息等级。可以设置成两类，越往下输出的信息更多。
    * `logging_level_INFO`
    * `logging_level_DEBUG`(default)
---
```python
def small_roots(self, fs: list, bounds: list, i: int, u = 1, ROOTS = [], HsFilter = []) -> list:
    """
        :param fs:          The original system of equations for which the roots need to be found
        :param bounds:      Upper bounds for all the roots
        :param i:           Sets the size of the Coppersmith lattice, related to the power of the modulus, m = i * n, where n is the number of fs
        :param ROOTS:       Correct solutions, can be provided during experiments to quickly determine if the solutions exist in the lattice after LLL
        :param HsFilter:    Filters the polynomials obtained after LLL, for example, setting it to [0,1,2,3] will select the first 4 polynomials
        :param u:           Needs to be provided when mode is set to mode_MODUP
        """
```
`small_roots()`是`Coppersmith`的主方法。
* `fs`: 需要求解的原始方程组。
* `bounds`: 所有变元的上界。
* `i`: 用来设置格的大小，`i`越大格越大。
* `ROOTS`(Optional): 设置正确解 ，一般是做实验时用来检验正确性。不知道正确解直接不填即可。
* `HsFilter`(Optional): LLL以后的矩阵的过滤器。具体为选择LLL后的矩阵的哪几行，例如设置为`[0,1,2,3]`则最后用于求根的方程只会用前四行。
* `u`(Optional): 当模式设置为`mode_MODUP`时需要设置P的倍数u。

---

```python
def __exgcd(self,mons: list):
```
单项式之间的扩展欧几里得算法

---
```python
def __digCalculate(self):
```
计算格的对角线上元素的乘积，当代码等级设置为`code_level_EXP`时会使用。

---
```python
def __ConstructM(self):
```
构造Coppersmith方法使用的单项式集合.

---
```python
def __ConstructF(self):
```
构造Coppersmith方法使用的移位多项式集合。

---
```python
def __ConstructB(self):
```
构造Coppersmith方法使用的格基矩阵
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