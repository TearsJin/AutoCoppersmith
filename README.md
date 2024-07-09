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

## Findroots.py
里面主要是一个整数方程求根用的类`class rootsFinder`，求根有三种可选的求根方法，求根方法在`rfconfig`中进行配置。三种求根方法：
```
    FConfig.method_GROEBNER: 利用整数上的groebner基求根
    
    RFConfig.method_CRT:     利用Zq上的groebner基求根
     
    RFConfig.method_VARIETY: 用VARIETY方法直接求根

```
---

```python
def __init__(self,rfconfig):
    """
        Initialize rootsFinder with configuration provided by RFConfig.

        :param rfconfig: Configuration for root finding details.
    """
```
* `rfconfig`: 根据`rfconfig`里的配置进行初始化，如果没有配置则使用默认配置

---

```python
def find_roots(self,R,Hs,bounds = []):
    """
        Main method for root finding.

        :param R:       The ring where roots are located.
        :param Hs:      Polynomials obtained after LLL.
        :param bounds:  Upper bounds of roots.
    """
```
* `R`: 根所在的整数环，主要是为了从中获取变元的个数。
* `Hs`: LLL以后的矩阵构造出的多项式组，也就是用来求根的整数方程组。
* `bounds`: 根的上界。




## Config.py

Config中暂时有三个类

```
class RFConfig: 用来配置求根方法，现在有三个:
    method_Groebner:   
    method_CRT:         
    method_VARIETY: 

class ETConfig: JM扩展策略，可以设置Ts来额外的shift指定的变元，例如对于双变元方程f(x,y)，设置Ts = [0,1]会让y多shift一倍的m。

class ULConfig: 用来配置拆分线性化方法。
```

---

```python

class RFConfig:
    def __init__(self,method = method_GROEBNER,gbLimitNum = 5,fgbThreads = 2,crtRootSign = []) -> None:
            """
                Initialization of RFConfig

                :param method:      Sets the root-finding method, options are method_GROEBNER, method_CRT, method_VARIETY
                :param gbLimitNum:  Sets the minimum number of polynomials used in GROEBNER root-finding
                :param fgbThreads:  Sets the number of threads for fgb_sage root-finding
                :param crtRootSign: Sets the sign of the final root in CRT root-finding; if not set, all cases will be output
            """
```
* `method`: 设置求根方法。
    * `method_GROEBNER`: 利用整数上的groebner基求根
    * `method_CRT`:      利用Zq上的groebner基求根
    * `method_VARIETY`:  用VARIETY方法直接求根B
* `gbLimitNum`(Optional): 设置GB基求根时使用的多项式数量的最小值（因为GB基求根时用的多项式太少会非常慢）
* `fgbThreads`(Optional): 设置fgb_sage计算GB基的线程数。(装了fgb_sage才能用)
* `crtRootSign`(Optional): 设置`method_CRT`求根时求得的根的符号，如果不设置则会输出所有可能的根。(Zq上求根时永远都是整数，所以如果根是负数的时候求得的根需要转换一下)

---

```python
class ETConfig:
    """
        Class for setting details of extension strategy.
        
        Attributes:
            Ts: The additional shifts in the extension strategy, for example, when there are three variables, setting [0, 0.2, 0.7] will shift the second variable by 0.2 * m times, and the third variable by 0.7 * m times.
    """
    def __init__(self, Ts = []) -> None:
        self.Ts = Ts
```
* `Ts`: 扩展策略中的额外偏移，例如，当有三个变量时，设置[0，0.2，0.7]将使第二个变量偏移0.2 * m倍，第三个变量偏移 0.7 * m倍。

---

```python
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
```
* `qr`: 拆分线性化中的新变元之间的关系。
* `unqr`: 拆分线性化中旧变元和新变元的代换关系。
* `bounds`: 旧变元的根的上界。
使用例子
```python
    PR1.<u1,u2,u3> = PolynomialRing(ZZ,3)
    PR2.<x,y> = PolynomialRing(ZZ,2)
    qr = u2 ** 2 - u1 - u3            #  u2 ** 2 = u1 + u3        
    unqr = [x ** 2 - y, x , y]          #  u1 = x ** 2 - y, u2 = x, u3 = y 
    ulconfig = ULConfig(qr,unqr,bounds)
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