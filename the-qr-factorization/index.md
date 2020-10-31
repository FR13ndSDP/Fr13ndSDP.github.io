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


当使用交错网格，此时界面上的速度矢量是已知的，上式可以直接计算出来。当使用同位网格，则需要界面上的插值技术如**Rhie-Chow**[^1]插值来避免棋盘型压力分布。

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
这组成了一个超定方程组，用矩阵表示为

$$
\begin{bmatrix} \Delta x_{1} & \Delta y_{1} \\\\ \Delta x_{2} & \Delta y_{2} \\\\ \Delta x_{3} & \Delta y_{3} \\\\ \vdots & \vdots \\\\ \Delta x_{N} & \Delta y_{N} \end{bmatrix}
\begin{bmatrix} \frac{\partial \phi}{\partial x}\big|_0 \\\\ \frac{\partial \phi}{\partial y}\big|_0 \end{bmatrix} \simeq
\begin{bmatrix} \phi_{1}-\phi_{0} \\\\ \phi_{2}-\phi_{0} \\\\ \phi_{3}-\phi_{0} \\\\ \vdots \\\\ \phi_{N}-\phi_{0} \end{bmatrix}
$$

那么两个方向上的梯度应该取多少使得误差最小呢，也就是求 

$$
\mathop{min} \limits_{x}||\mathbf{A}x-b||^2
$$

 这其实是一个线性最小二乘问题，也就是这个方法名字的由来。一种解法是 $x = \mathbf{(A^{T}A)^{-1}A^TB}$，然而这可能导致精度上的重大缺陷，首先，矩阵求逆的精度较低，计算量较大，其次，也是更重要的原因：从**条件数**[^2]判断，$\mathbf{A^TA}$的条件数是$\mathbf{A}$条件数的平方。这种做法会使得方程比原来的超定方程更加病态（badly conditioned），在有限的浮点数精度下，矩阵$\mathbf{A^TA}$可能是奇异的，求逆就无法进行。

## QR分解

{{< figure src="/images/The-QR-Factorization/qr.png" title="QR分解两种形式">}}

对$\mathbf{A}$进行QR分解$\mathbf{A = QR}$，$\mathbf{Q}$是正交矩阵，满足$\mathbf{Q^TQ =I}$，$\mathbf{R}$为上三角矩阵。

那么 $x = \mathbf{R^{-1}Q^T}b$，而上三角矩阵的求逆很简单。

现在问题是，如何进行分解呢？

### 1. Gram-Schmidt正交化方法

矩阵$\mathbf{A}$表示为$\mathbf{A} = \boldsymbol{[v_1, v_2,...,v_n]}$，将这些基底向量正交化：
$$
\boldsymbol{\beta_1 = v_1}\\\\ \boldsymbol{\beta_2=v_2-\langle v_2,\eta_1\rangle \eta_1}\\\\ \vdots\\\ \boldsymbol{\beta_n=v_n-\sum_{i=1}^{n-1}\langle v_n, \eta_i\rangle \eta_i}
$$

其中的$\boldsymbol{\eta_n = \frac{\beta_n}{||\beta_n||}}$即标准正交基。那么有
$$
\mathbf{A} = \boldsymbol{[v_1,...,v_n] = [\eta_1,...,\eta_2] \begin{bmatrix} ||\beta_1||& \frac{\langle v_2,\beta_1\rangle}{||\beta_1||}& ...&\frac{\langle v_n,\beta_1\rangle}{||\beta_1||}\\\\ 0 & ||\beta_2||& ... & \frac{\langle v_n,\beta_2\rangle}{||\beta_2||}\\\\ \vdots& \vdots& ... & \vdots\\\\ 0 & 0 & ... & ||\beta_n|| \end{bmatrix}}
$$
显然，$\boldsymbol{\eta_i}$组成了正交矩阵，又由$\boldsymbol{||\beta_i||}\neq 0$，$\mathbf{R}$是一个主对角元非零的上三角矩阵。

### 2. Householder变换

> 相比与Gram-Schmidt正交化，使用Householder变换具有更好的[数值稳定性](https://zh.wikipedia.org/wiki/数值稳定性)






### 稳定性分析

---
[^1]: 在同位网格的界面采用特殊的速度插值来避免棋盘型压力分布。特征是引入一个三阶的压力梯度项来提供阻尼（人工扩散）抑制压力震荡。
[^2]: $\kappa(\mathbf{A}) = \frac{\sigma_{max}(\mathbf{A})}{\sigma_{min}(\mathbf{A})}$，也就是矩阵的极大奇异值和极小奇异值之比，一个问题的条件数是该数量在数值计算中的容易程度的衡量，也就是该问题的适定性。一个低条件数的问题称为良置的，而高条件数的问题称为病态（或者说非良置）的。[维基百科](https://zh.wikipedia.org/wiki/%E6%9D%A1%E4%BB%B6%E6%95%B0)


