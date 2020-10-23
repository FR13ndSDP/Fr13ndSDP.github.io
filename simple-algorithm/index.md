# Simple algorithm


SIMPLE 算法在一维稳态无粘不可压流动问题中的应用

<!--more-->

#### 1.问题描述

![](/images/simple-algorithm/problem.jpg)

- 流体密度 $\rho = 1\ kg/m^3$, 不可压缩，无粘流动
- 喷管长度 $L = 0.2\ m$
- 喷管左右面积给定 $A_{A} = 0.5\ m^2, A_{E} = 0.1\ m^2$
- 左边界：恒定动压 $10\ Pa$，右边界：恒定静压$0\ Pa$

#### 2.控制方程

- 质量守恒方程
  
  $$
  \frac{d}{dx}(\rho A u) = 0
  $$
  
- 动量方程
  $$
  \rho u A \frac{du}{dx} = -A\frac{dp}{dx}
  $$

解析解可由伯努利定律得到
$$
p_0 = p_{E} + 1/2\rho M^2/(\rho A_E)^2
$$
其中$M$为质量流量，可得到精确解为$0.44271\ kg/s$

#### 3.SIMPLE算法

Semi-Implicit Method for Pressure Linked Equations.

这是一个求解速度-压力耦合的经典迭代算法

实施过程中需要注意的点：

- 边界条件的处理，在本例中入口边界使用了延迟修正（deferred correction）增强迭代稳定性
- 松弛因子的选取，一般取$\alpha_u + \alpha_p = 1$，压力松弛因子要取得相对小避免发散
- 速度的松弛放在离散动量方程的求解中，压力松弛在解出压力修正值之后

带有速度松弛的离散动量方程：
$$
\frac{a_{i,J}}{\alpha_u}u_{i,J} = \sum a_{nb}u_{nb}+(p_{I-1,J} - p_{I,J})A_{i,J}+b_{i,J}+\frac{(1-\alpha_u)a_{i,J}}{\alpha_u}u_{i,J}^{(n-1)}
$$
下面是一个代码不规范但是能说明具体实施过程的例程，对流项使用一阶迎风差分。


```python
#!/usr/bin/python3
# Copyright (c) 2020 Fr13ndSDP All rights reserved.
import numpy as np
import matplotlib.pyplot as plt

# the TDMA
def sol(aw, ap, ae, su):
    n = np.size(aw)
    p = np.zeros(n)
    q = np.zeros(n)
    phi = np.zeros(n)
    p[0] = ae[0]/ap[0]
    q[0] = su[0]/ap[0]
    for i in range(1, n):
        p[i] = ae[i]/(ap[i] - aw[i]*p[i-1])
        q[i] = (su[i] + aw[i]*q[i-1])/(ap[i]-aw[i]*p[i-1])

    phi[n-1] = q[n-1]
    for i in range(n-1, 0, -1):
        phi[i-1] = p[i-1]*phi[i]+q[i-1]

    return phi

rho = 1  # density
L = 2   # length
n = 100  # number of pressure nodals
nu = n-1    # number of velocity nodals
dx = L/n
p0 = 10  # inlet stagnation pressure
p_exit = 0  # oulet gauge pressure
mass = []   # mass flow
sp = [1]    # momentum residual (source of pressure Eq.)
m = 1  # guessed mass flow

# matrix coeff and source term coeff
apu = np.zeros(nu)
awu = np.zeros(nu)
aeu = np.zeros(nu)
suu = np.zeros(nu)
app = np.zeros(n)
awp = np.zeros(n)
aep = np.zeros(n)
sup = np.zeros(n)

# under-relaxation factor
alpha_u = 0.8
alpha_p = 0.3

d = np.zeros(nu)  # A/ap

# area of pressure node and velocity node
Ap = np.linspace(0.5, 0.1, n)
Au = np.linspace(0.5-0.25/n, 0.1+0.25/n, nu)

# initial guessed value
u = np.zeros(nu)
for i in range(nu):
    u[i] = m/(rho*Au[i])

p = np.linspace(10, 0, n)
pf = np.zeros(n)  # Pressure correction value

while(max(sp) > 1e-5):
    sp = []
    # velocity update
    for i in range(1, nu-1):
        awu[i] = rho*(u[i-1]+u[i])/2*Ap[i]  # F_w
        aeu[i] = 0  # upwind scheme
        apu[i] = rho*(u[i]+u[i+1])/2 * Ap[i+1]  # awu[i]+aeu[i]+ (F_e-F_w)
        apu[i] = apu[i]/alpha_u  # velocity relaxation
        suu[i] = Au[i]*(p[i]-p[i+1]) + (1-alpha_u)*apu[i]*u[i] 
        d[i] = Au[i]/apu[i]

    apu[0] = rho*(u[0]+u[1])/2*Ap[1] + rho * \
        (u[0]*Au[0]/Ap[0])*Ap[0]*0.5*(Au[0]/Ap[0])**2
    apu[0] = apu[0]/alpha_u  # velocity relaxation
    suu[0] = Au[0]*(p0 - p[1])+rho*(u[0]*Au[0]/Ap[0])*Ap[0] * \
        (Au[0]/Ap[0])*u[0]+(1-alpha_u)*apu[0]*u[0] #deffered correction
    d[0] = Au[0]/apu[0]

    awu[nu-1] = rho*(u[nu-2]+u[nu-1])/2*Ap[nu-1]
    apu[nu-1] = rho*u[nu-1]*Au[nu-1]
    apu[nu-1] = apu[nu-1]/alpha_u  # velocity relaxation
    suu[nu-1] = Au[nu-1]*(p[nu-1] - p[nu]) + (1-alpha_u)*apu[nu-1]*u[nu-1]
    d[nu-1] = Au[nu-1]/apu[nu-1]

    u = sol(awu, apu, aeu, suu)

    # pressure update
    for i in range(1, n-1):
        awp[i] = rho*d[i-1]*Au[i-1]
        aep[i] = rho*d[i]*Au[i]
        app[i] = awp[i]+aep[i]
        sup[i] = rho*u[i-1]*Au[i-1] - rho*u[i]*Au[i]
        sp.append(sup[i])

    # do not update boundary pressure
    # so set app to none-zero number and set sup to zero
    app[0] = 1
    sup[0] = 0
    app[n-1] = 1
    sup[n-1] = 0
    pf = sol(awp, app, aep, sup)

    # under-relaxation of pressure
    prev = p
    p = p + pf
    p[0] = p0 - 0.5*rho*u[0]*u[0]*(Au[0]/Ap[0])**2
    p = (1-alpha_p)*prev + alpha_p*p

    for i in range(nu):
        u[i] = u[i]+d[i]*(pf[i]-pf[i+1])

    # calculate mass flow
    mass.append(rho*Au[0]*u[0])

plt.figure(1)
plt.plot(mass)
plt.figure(2)
plt.plot(p)
plt.figure(3)
plt.plot(u)
print("mass flow:", mass[-1])
print("momentum residual:", max(sp))
```

> mass flow: 0.45496785377760623
>
> momentum residual: 8.283464151048747e-06



![mass flow](/images/simple-algorithm/output_1_1.png)



![pressure](/images/simple-algorithm/output_1_2.png)



![velocity](/images/simple-algorithm/output_1_3.png)

采用压力修正方程中源项的绝对值的最大值作为判断收敛依据，采用100个压力节点得到收敛解：质量流量$0.455\ kg/s$，并得到了压力、速度的分布。



#### 3.参考文献

- [1]  陶文铨, 数值传热学
- [2] H.K.Versteeg, An introduction to computational fluid dynamics : the finite volume method
- [3] J.H.Ferziger, M.Peric, Computational method for fluid dynamics
