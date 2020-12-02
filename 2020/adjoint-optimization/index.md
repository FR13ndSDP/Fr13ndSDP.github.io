# Adjoint Optimization


伴随优化基本方程的推导以及OpenFOAM中的实现

<!--more-->

> 使用如下的符号记法：全导数 $d_x$ $(\nabla_x)$，偏导 $\partial_x$，微分 $d$.

## 1 伴随方法

在一个偏微分方程系统中，假设存在变量 $x\in R^{n_x}, p\in R^{n_p}$，方程 $f(x,p):R^{n_x}\times R^{n_p}\rightarrow R$以及关系 $g(x,p) = 0$，且 $g_x$ 处处非奇异，$d_pf$怎么求？

### 1.1 动机

$g(x,p)=0$ 的求解是CFD求解器的核心，给定参数 $p$，程序将计算出 $x$。例如 $p$ 可以是边界条件或初始条件的参数，也可以是物性参数，而 $x$ 是计算得到的场量，这种求解模式是*正问题*。如果引入 $f(x,p)$作为度量标准，例如用来计算 $x$ 的光滑程度，那么通常需要最小化$f$，这种求解模式是*反问题*。

梯度 $d_pf$ 很有用，可以用来计算优化问题 $min_pf$ 中 $f$ 对于参数 $p$ 变化的敏感程度，并使用梯度下降方法求解最优化问题。

### 1.2 推导

1）

考虑一个简单的函数 $f(x)$，首先由链式法则有
$$
d_pf = d_df(x(p)) = \partial_xfd_px\ (=f_xx_p)
$$
其次因为 $g(x,p) = 0$ 处处成立，则 $d_pg = 0$， 也就是
$$
g_xx_p + g_p = 0
$$
那么得到
$$
d_pf = -f_xg_x^{-1}g_p
$$
从线性代数的观点看，$f_xg_x^{-1}$ 是一个行向量乘以 $n_x\times n_x$ 的矩阵，并且是以下方程的解
$$
g_x^{T}\lambda = -f_x^{T}
$$
称为伴随方程，$\lambda$称为伴随变量。得到 $d_pf = \lambda^Tg_p$

2）

定义拉格朗日函数
$$
L(x,p,\lambda) = f(x) + \lambda^Tg(x,p)
$$
这里$\lambda$是拉格朗日乘子组成的向量，由于 $g(x,p)$ 处处为零，$\lambda$可以随意选取。
$$
d_pf(x) = d_pL = f_xx_p + d_p\lambda^Tg + \lambda^T(g_xx_p+g_p)\\\\ =(f_x+\lambda^Tg_x)x_p+\lambda^Tg_p
$$
如果选取 $g^T\lambda = -f_x^T$，第一项为零，可以避免计算 $x_p$ 并且得到 $d_pf = \lambda^Tg_p$，与第一种推导方式相同。

或者，由于$df = f_xdx +f_pdp = (f_x+\lambda^Tg_x )dx+ \lambda^Tg_pdp$，可以推至同样的结果。

## 2 PDE约束的连续伴随问题

### 2.1 一般伴随方程

给出一个特定的优化问题
$$
min\quad J = J(\mathbf{v}, p, \alpha)\\\\ s.t.\quad R(\mathbf{v},p,\alpha) = 0
$$
考虑由连续介质力学给出的约束，$\mathbf{v}$ 与 $p$ 作为变量，即等价为前文中的 $x$; $\alpha$ 作为 设计参数，即等价于前文中的$p$，给出问题
$$
min\quad J = J(\mathbf{v},p,\alpha)\\\\ s.t. \quad R^{u} = \nabla\cdot(\mathbf{vv}) + \nabla p - \nabla\cdot(2\nu D(\mathbf{v}))+\alpha \mathbf{v} = 0\\\\ R^{p} = \nabla \cdot \mathbf{v} = 0
$$
上式由达西定律引入了渗透率作为源项，假设某一区域使得目标函数增大，可以通过减小这一区域的渗透率作为惩罚。$D(\mathbf{v})=\frac{1}{2}(\nabla \mathbf{v}+\nabla \mathbf{v}^T)$ 代表应变率张量。

使用拉格朗日乘数法将其转变为无约束优化问题
$$
min \quad L = J+\int_\Omega(\mathbf{u},q)Rd\Omega\\\\ = J+ \int_{\Omega}qR^pd\Omega + \int_\Omega{\mathbf{u}\cdot R^{u}}d\Omega
$$
其中$(\mathbf{u},q)$为拉格朗日乘子，称其为伴随速度和伴随压力（后面会看到为什么），但是实际上并没有物理含义。

因为
$$
dJ = dL = L_{\alpha}d{\alpha}+ L_{\mathbf{v}}d{\mathbf{v}}+L_{p}dp
$$
注意这里没有将粘度作为微分变量，是一种近似的做法，被称为“冻结湍流”。

由于 $(\mathbf{u},q)$可以自由选取，取 适当的值使得 $L_{\mathbf{v}}d{\mathbf{v}}+L_pdp = 0$ 成立，便得到
$$
\frac{\partial J}{\partial \alpha} = L_{\alpha} = J_{\alpha} + \int_{\Omega}\mathbf{u}\cdot \mathbf{v}d\Omega
$$
当考虑网格中的目标函数关于 $\alpha$的梯度，可以得到
$$
\frac{\partial J}{\partial \alpha_i} = L_{\alpha_i} = J_{\alpha_i} + \mathbf{u_i}\cdot \mathbf{v_i}V_i
$$
其中 $V_i$ 为网格体积。

如果 $L_{\mathbf{v}}d{\mathbf{v}}+L_pdp = 0$ 成立，那么伴随方程
$$
\delta_{\mathbf{v}}J + \delta_{p}J +
\int_{\Omega}\mathbf{u}\cdot [\nabla\cdot(\mathbf{v\delta v})+\nabla\cdot(\mathbf{\delta v v})   - \nabla\cdot(2\nu D(\delta\mathbf{v}))+\alpha \delta \mathbf{v}]d\Omega \\\\ - \int_{\Omega}q\nabla \cdot \delta \mathbf{v}d\Omega + \int_{\Omega}\mathbf{u}\cdot \nabla \delta pd\Omega = 0
$$
积分使用散度公式与高斯散度定理(以及一些张量运算的推导)：
$$
\nabla\cdot (q\mathbf{v}) = \nabla q \cdot \mathbf{v}+q\nabla \cdot \mathbf{v}\\\\ \int_{\Omega}\nabla \cdot (q\mathbf{v})\ d\Omega =\int_{\Gamma}q\mathbf{ v}\cdot \mathbf{n}d\Gamma
$$
并且把 $J$ 拆分为 
$$
J = \int_{\Gamma}J_{\Gamma}\ d\Gamma + \int_{\Omega} J_{\Omega}\ d\Omega
$$
得到  (*注意*  $(\mathbf{v}\cdot \nabla)\mathbf{u}=\nabla(\mathbf{vu})$）
$$
\int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{u} \cdot \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial p}\right) \delta p+\int_{\Omega} \mathrm{d} \Omega\left(-\nabla \cdot \mathbf{u}+\frac{\partial J_{\Omega}}{\partial p}\right) \delta p\\\\ \quad+\int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{n}(\mathbf{u} \cdot \mathbf{v})+\mathbf{u}(\mathbf{v} \cdot \mathbf{n})+2 v \mathbf{n} \cdot \mathbf{D}(\mathbf{u})-q \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial \mathbf{v}}\right) \cdot \delta \mathbf{v}-\int_{\Gamma} \mathrm{d} \Gamma 2 v \mathbf{n} \cdot \mathbf{D}(\delta \mathbf{v}) \cdot \mathbf{u}\\\\ \quad+\int_{\Omega} \mathrm{d} \Omega\left(-\nabla \mathbf{u} \cdot \mathbf{v}-(\mathbf{v} \cdot \nabla) \mathbf{u}-\nabla \cdot(2 v \mathbf{D}(\mathbf{u}))+\alpha \mathbf{u}+\nabla q+\frac{\partial J_{\Omega}}{\partial \mathbf{v}}\right) \cdot \delta \mathbf{v}=0
$$

观察上面的式子，他应该对任意的 $\delta \mathrm{p}$和 $\delta \mathbf{v}$ 成立，因此各个积分分别等于0，当在 $\Omega$ 内积分时，边界积分为零，因此得到 $\Omega$ 内的伴随方程
$$
-\nabla (\mathbf{vu})-\nabla \mathbf{u}\cdot \mathbf{v} = -\nabla{q}+\nabla \cdot (2\nu D(\mathbf{u}))-\alpha \mathbf{u} - \frac{\partial J_{\Omega}}{\partial \mathbf{v}}\\\\ \nabla \cdot \mathbf{u} = \frac{\partial J_{\Omega}}{\partial p}
$$
可以看到，这组方程与 N-S 方程非常相似，因此变量 $(\mathbf{u}, q)$ 被看作是非物理的速度和压力。

### 2.2 边界条件

在边界处，可以得到边界条件为
$$
\int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{n}(\mathbf{u} \cdot \mathbf{v})+\mathbf{u}(\mathbf{v} \cdot \mathbf{n})+2 v \mathbf{n} \cdot \mathbf{D}(\mathbf{u})-q \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial \mathbf{v}}\right) \cdot \delta \mathbf{v}-\int_{\Gamma} \mathrm{d} \Gamma 2 v \mathbf{n} \cdot \mathbf{D}(\delta \mathbf{v}) \cdot \mathbf{u}=0\\\\ \int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{u} \cdot \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial p}\right) \delta p=0
$$
或者
$$
\int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{n}(\mathbf{u} \cdot \mathbf{v})+\mathbf{u}(\mathbf{v} \cdot \mathbf{n})+2 v \mathbf{n} \cdot \mathbf{D}(\mathbf{u})-q \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial \mathbf{v}}\right) \cdot \delta \mathbf{v}=0\\\\ \int_{\Gamma} \mathrm{d} \Gamma\left(\mathbf{u} \cdot \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial p}\right) \delta p + \int_{\Gamma} \mathrm{d} \Gamma 2 v \mathbf{n} \cdot \mathbf{D}(\delta \mathbf{v}) \cdot \mathbf{u}=0
$$


1） 在入口边界以及壁面处，速度值为给定值或者0，因此 $\delta \mathbf{v} = 0$，因此部分积分可以视为0，为了避免无解，显然应该选取第一种边界条件的定义。

因此
$$
\int_{\Gamma} \mathrm{d} \Gamma 2 v \mathbf{n} \cdot \mathbf{D}(\delta \mathbf{v}) \cdot \mathbf{u}=0 \\\\ u_n = -\frac{\partial J_{\Gamma}}{\partial p}
$$
其中 $u_n = \mathbf{u\cdot n}$ 为伴随速度在边界法向的分量。现在需要定义$u_t$ 在边界处的取值，这里需要用到一些近似，具体推导参看参考文献[^1]，这里直接给出入口及壁面处应该满足的边界条件：
$$
u_t = 0\\\\ u_n = -\frac{\partial J_{\Gamma}}{\partial p}\\\\ \mathbf{n}\cdot \nabla q = 0
$$
2) 出口边界常见设为压力为零和速度零梯度，因此除了自动满足的积分，边界条件剩下
$$
\mathbf{n}(\mathbf{u} \cdot \mathbf{v})+\mathbf{u}(\mathbf{v} \cdot \mathbf{n})+v(\mathbf{n} \cdot \nabla) \mathbf{u}-q \mathbf{n}+\frac{\partial J_{\Gamma}}{\partial \mathbf{v}}=0
$$

## 3 特定问题—目标函数为能量损失

选取目标函数为流动的能量损失，即最小化能量损失
$$
J = -\int_{\Gamma}(p+\frac{1}{2}v^2)\mathbf{v}\cdot \mathbf{n}\ d\Gamma
$$
负号是因为考虑到通量的方向，即最小化进出口能量之差。

$J_{\Omega} = 0, J_{\Gamma} = -(p+\frac{1}{2}v^2)\mathbf{v\cdot n}$，因此
$$
\frac{\partial J_{\Gamma}}{\partial p} = -\mathbf{v\cdot n},\quad \frac{\partial J_{\Gamma}}{\partial \mathbf{v}} = \frac{\partial(-p\mathbf{v\cdot n}-\frac{1}{2}\mathbf{(v\cdot v)v\cdot n)}}{\partial\mathbf{v}} = -p\mathbf{n}-(\mathbf{v\cdot n})\cdot \mathbf{v}
$$


代入伴随方程，得到在 $\Omega$ 内
$$
-\nabla (\mathbf{vu})-\nabla \mathbf{u}\cdot \mathbf{v} = -\nabla{q}+\nabla \cdot (2\nu D(\mathbf{u}))-\alpha \mathbf{u}\\\\ \nabla \cdot \mathbf{u} = 0
$$
入口边界条件：
$$
u_t = 0\\\\ u_n = v_n\\\\ \mathbf{n}\cdot \nabla q = 0
$$
出口边界条件：
$$
\mathbf{n}(\mathbf{u} \cdot \mathbf{v})+\mathbf{u}(\mathbf{v} \cdot \mathbf{n})+v(\mathbf{n} \cdot \nabla) \mathbf{u}-q \mathbf{n}-p\mathbf{n}-(\mathbf{v\cdot n})\cdot \mathbf{v}=0
$$

## 4. OpenFOAM中的实现

## 参考资料

[^1]:C. Othmer, A continuous adjoint formulation for the computation of topological and surface sensitivities of ducted flows
[^2]:Andrew M. Bradley, PDE-constrained optimization and the adjoint method
[^3]:Topology Optimisation of Fluids Through the Continuous Adjoint Approach in OpenFOAM
[^4]:Ulf Nilsson, Description of adjointShapeOptimizationFoam and how to implement new objective functions
