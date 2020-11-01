# The QR Factorization


复习一下线性代数

<!--more-->

## 背景

当在非结构网格上离散对流项时，利用高斯散度定理可以得到
$$
\int_{CV}\nabla \cdot (\rho \phi \mathbf{u})dV = \sum_{all\ surfaces}\int_{\Delta A_i}\mathbf{n_i}\cdot (\rho \phi \mathbf{u})dA
$$
定义一个对流通量$F_i$
$$
F_i = \int_{\Delta A_i}\mathbf{n_i}\cdot (\rho \mathbf{u})dA \simeq \mathbf{n_i}\cdot(\rho \mathbf{u})\Delta A_i
$$


当使用交错网格，此时界面上的速度矢量是已知的，上式可以直接计算出来。当使用同位网格，则需要界面上的插值技术如**Rhie-Chow插值**[^1]来避免棋盘型压力分布。

对流项可以表示为
$$
\sum_{all\ surfaces}F_i \phi_i
$$
$\phi_i$是界面上的物理量，需要满足守恒性、有界性和有向性。在非结构网格上使用一阶迎风差分是不可能的，结构网格下的一阶迎风格式在流向不与坐标轴平行时会引起假扩散。这在连正交都不能保证的非结构网格中将引入很强的假扩散，如果采用线性迎风差分或者是TVD格式，就必须要计算$\phi$的梯度（Darwish and Moukalled, 2003）。

一种可行的方法是**最小二乘梯度重建**（least-squares gradient reconstruction），参照下图：

{{< figure src="/images/The-QR-Factorization/grid.png" title="一个非结构网格">}}

$$
\phi_i - \phi_0 \simeq \left(\frac{\partial \phi}{\partial x}\right)\bigg|_0\Delta x_i + \left(\frac{\partial \phi}{\partial y}\right)\bigg|_0\Delta y_i
$$
这组成了一个超定方程，用矩阵表示为

$$
\begin{bmatrix} \Delta x_{1} & \Delta y_{1} \\\\ \Delta x_{2} & \Delta y_{2} \\\\ \Delta x_{3} & \Delta y_{3} \\\\ \vdots & \vdots \\\\ \Delta x_{N} & \Delta y_{N} \end{bmatrix}
\begin{bmatrix} \frac{\partial \phi}{\partial x}\big|_0 \\\\ \frac{\partial \phi}{\partial y}\big|_0 \end{bmatrix} \simeq
\begin{bmatrix} \phi_{1}-\phi_{0} \\\\ \phi_{2}-\phi_{0} \\\\ \phi_{3}-\phi_{0} \\\\ \vdots \\\\ \phi_{N}-\phi_{0} \end{bmatrix}
$$

那么两个方向上的梯度应该取多少使得误差最小呢，也就是求 

$$
\mathop{min} \limits_{x}||\boldsymbol{Ax-b||^2}
$$

 这其实是一个线性最小二乘问题，也就是这个方法名字的由来。一种解法是 $\boldsymbol{x = (A^{T}A)^{-1}A^TB}$，然而这可能导致精度上的重大缺陷，首先，矩阵求逆的精度较低，计算量较大，其次，也是更重要的原因：从**条件数**[^2]判断，$\boldsymbol{A^TA}$的条件数是$\boldsymbol{A}$条件数的平方。这种做法会使得方程比原来的超定方程更加病态（badly conditioned），在有限的浮点数精度下，矩阵$\boldsymbol{A^TA}$可能是奇异的，求逆就无法进行。

## QR分解

{{< figure src="/images/The-QR-Factorization/qr.png" title="QR分解两种形式">}}

对$\boldsymbol{A}$进行QR分解$\boldsymbol{A = QR}$，$\boldsymbol{Q}$是正交矩阵，满足$\boldsymbol{Q^TQ =I}$，$\boldsymbol{R}$为上三角矩阵。

那么 $\boldsymbol{x = R^{-1}Q^Tb}$，而上三角矩阵的求逆很简单。

现在问题是，如何进行分解呢？

### 1. Gram-Schmidt正交化方法

矩阵$\boldsymbol{A}$表示为$\boldsymbol{A} = \boldsymbol{[v_1, v_2,...,v_n]}$，将这些基底向量正交化：
$$
\boldsymbol{\beta_1 = v_1}\\\\ \boldsymbol{\beta_2=v_2-\langle v_2,\eta_1\rangle \eta_1}\\\\ \vdots\\\ \boldsymbol{\beta_n=v_n-\sum_{i=1}^{n-1}\langle v_n, \eta_i\rangle \eta_i}
$$

其中的$\boldsymbol{\eta_n = \frac{\beta_n}{||\beta_n||}}$即标准正交基。那么有
{{< figure src="/images/The-QR-Factorization/matrix.png">}}
显然，$\boldsymbol{\eta_i}$组成了正交矩阵，又由$\boldsymbol{||\beta_i||}\neq 0$，$\boldsymbol{R}$是一个主对角元非零的上三角矩阵。

### 2. Householder变换

> 相比与Gram-Schmidt正交化，使用Householder变换具有更好的[数值稳定性](https://zh.wikipedia.org/wiki/数值稳定性)

Householder变换是将一个向量关于某个$n-1$维子空间进行反射的操作。首先说明是如何进行反射的（暂时不考虑复矩阵）

#### 2.1 怎么旋转

Householder变换矩阵$\boldsymbol{H} = \boldsymbol{I} - \boldsymbol{2vv^T}$，这个变换具有性质$\boldsymbol{H^T = H, H^{-1} = H^T}$，即是对称正交的。其中$\boldsymbol{v}$是平行于“镜面”的单位向量，那么
$$
\boldsymbol{Hx = x-2v\langle v,x\rangle}
$$

{{< figure src="/images/The-QR-Factorization/householder.png" title="Householder反射变换">}}

很清楚，Householder就是一个反射变换。

#### 2.2 怎么应用

为了将他用在QR分解上，就得考虑如何将这个变换作用在矩阵$\boldsymbol{A}$上从而得到$\boldsymbol{Q}$，也就是说通过变换将部分基向量映射到一组维度递减的单位向量$\boldsymbol{e_1,e_2,...,e_n}$上，例如$\boldsymbol{e_1 = [1,0,0]^T, e_2 = [1,0]^T}$，得到$\boldsymbol{\Eta A = R}$，$\boldsymbol{\Eta}$将会是一系列正交变换的组合，从而也是正交的，那么$\boldsymbol{Q = \Eta ^T}$。

为了将一个向量$\boldsymbol{x}$反射到$\boldsymbol{e_1}$，可选取“镜面向量”

$$
\boldsymbol{u = x-||x||e_1,\ v = \frac{u}{||u||}}
$$

接下来对$\boldsymbol{A}$的部分基向量依次通过变换$\boldsymbol{H_i}$反射到$\boldsymbol{e_1, e_2，...，e_n}$，生成矩阵$\boldsymbol{R}$:
$$
\boldsymbol{H_n...H_2H_1A = R}
$$
注意到$\boldsymbol{H_n}$比$\boldsymbol{H_{n-1}}$维度小，可令
$$
\boldsymbol{H_n} = \begin{pmatrix}\boldsymbol{I_{k}}& 0\\\\ 0& {\mathop{H}\limits^\sim}_n \end{pmatrix}
$$
其中${\mathop{H}\limits^\sim}_n$是真正对余子式起变换作用的。



因此$\boldsymbol{Q = [H_n...H_2H_1]^T = H_1H_2...H_n}$。然而在实际的求解过程中，不需要求出$\mathbf{Q}$，因为
$$
\boldsymbol{H_n...H_2H_1A}x = \boldsymbol{H_n...H_2H_1b}\\\\ \implies \boldsymbol{Rx = z}
$$
由于$\boldsymbol{R}$是上三角矩阵，那么只需要对这个方程逐步回代就可以求解。

### 3. 其它方法

QR分解还有Givens旋转等其他方法，上述的方法还有分类如Full QR Factorization、Reduced QR Factorization等，实现算法更是有许多的技巧，误差和稳定性分析也是一门大学问。但是我不想再深究下去了。

## Practice
```python
import numpy as np

def householder(A,b):
    r = A.shape[0]
    c = A.shape[1]
    Q = np.identity(r)
    R = np.copy(A)
    z = np.copy(b)
    for cnt in range(r - 1):
        x = R[cnt:, cnt]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e 
        v = u / np.linalg.norm(u) # normalized "mirror vector"
        Q_cnt = np.identity(r)
        Q_cnt[cnt:, cnt:] -= 2.0 * np.outer(v, v)
        R = np.dot(Q_cnt, R)  # R=H(n-1)...H(2)H(1)A
        Q = np.dot(Q, Q_cnt)  # Q=H(1)H(2)...H(n)
        z = np.dot(Q_cnt, z)  # z=H(n-1)...H(2)H(1)b 
   
    # reduction
    Q = Q[:r,:c]
    R = R[:c,:c]
    z = z[:c]
    return Q,R,z

A = np.array([[1, 5, 0],[5, -1, 4],[5, 1, -4],[0, 4, 1]],dtype=float)
m = A.shape[1]
b = np.array([[1],[2],[3],[4]])
Q,R,z= householder(A,b)
print(np.linalg.solve(R,z))
```
> [[0.4644]
> 
> [0.4628]
> 
> [0.0561]]
```matlab
a = [1,5,0;5,-1,4;5,1,-4;0,4,1];b=[1;2;3;4];a\b
```

> ans =
> 
> 0.4644 
> 
> 0.4628 
> 
> 0.0561 



[^1]: 在同位网格的界面采用特殊的速度插值来避免棋盘型压力分布。特征是引入一个三阶的压力梯度项来提供阻尼（人工扩散）抑制压力震荡。
[^2]: $\kappa(\mathbf{A}) = \frac{\sigma_{max}(\mathbf{A})}{\sigma_{min}(\mathbf{A})}$，也就是矩阵的极大奇异值和极小奇异值之比，一个问题的条件数是该数量在数值计算中的容易程度的衡量，也就是该问题的适定性。一个低条件数的问题称为良置的，而高条件数的问题称为病态（或者说非良置）的。[维基百科](https://zh.wikipedia.org/wiki/%E6%9D%A1%E4%BB%B6%E6%95%B0)

## 参考

- [1] G.Strang, Linear Algebra and its applications
- [2] [QR decomposition](https://en.wikipedia.org/wiki/QR_decomposition)
- [3] [豪斯霍尔德变换](https://zh.wikipedia.org/wiki/%E8%B1%AA%E6%96%AF%E9%9C%8D%E5%B0%94%E5%BE%B7%E5%8F%98%E6%8D%A2)
- [4] [[数值计算] QR分解](https://zhuanlan.zhihu.com/p/84415000)
- [5] C.Moler, MATLAB数值计算
