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
其中$M$为质量流量，可得到精确解为$0.44721\ kg/s$

#### 3.SIMPLE算法

Semi-Implicit Method for Pressure Linked Equations.

这是一个求解速度-压力耦合的经典迭代算法

实施过程中需要注意的点：

- 交错网格，速度网格落后于压力网格

- 边界条件的处理，在本例中入口边界使用了延迟修正（deferred correction）增强迭代稳定性
- 松弛因子的选取，一般取$\alpha_u + \alpha_p = 1$，压力松弛因子要取得相对小避免发散
- 速度的松弛放在离散动量方程的求解中，压力松弛在解出压力修正值之后

带有速度松弛的离散动量方程：
$$
\frac{a_{i,J}}{\alpha_u}u_{i,J} = \sum a_{nb}u_{nb}+(p_{I-1,J} - p_{I,J})A_{i,J}+b_{i,J}+\frac{(1-\alpha_u)a_{i,J}}{\alpha_u}u_{i,J}^{(n-1)}
$$
下面是一个代码不规范但是能说明具体实施过程的例程，对流项使用一阶迎风差分或者TVD格式（Van Leer） 。


```python
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

def vanLeer(u,n):
    re = (u[n] - u[n-1])/(u[n+1] - u[n])
    rw = (u[n-1] - u[n-2])/(u[n] - u[n-1])
    phie = (re+abs(re))/(1+re)
    phiw = (rw+abs(rw))/(1+rw)
    return phie,phiw

def minMod(u,n):
    re = (u[n] - u[n-1])/(u[n+1] - u[n])
    rw = (u[n-1] - u[n-2])/(u[n] - u[n-1])
    if (re >0):
        phie = min(re,1)
    else:
        phie = 0
    if (rw >0):
        phiw = min(rw,1)
    else:
        phiw = 0  
    return phie,phiw

def upwind(u,n):
    phie = 0
    phiw = 0
    return phie,phiw

rho = 1  # density
L = 2   # length
n = 100  # number of pressure nodals
nu = n-1    # number of velocity nodals
dx = L/n
p0 = 10  # inlet stagnation pressure
p_exit = 0  # oulet gauge pressure
mass = []   # mass flow
residual = []
sp = [1]    # momentum residual (source of pressure Eq.)
m = 0.5  # guessed mass flow

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
    for i in range(2, nu-1):
        fw = rho*(u[i-1]+u[i])/2*Ap[i]
        fe = rho*(u[i]+u[i+1])/2*Ap[i+1]
        phie,phiw = upwind(u,i)
        awu[i] = fw*(1-0.5*phiw)
        aeu[i] = -0.5*fe*phie  # upwind scheme
        apu[i] = awu[i]+aeu[i]+fe-fw  # awu[i]+aeu[i]+ (F_e-F_w)
        apu[i] = apu[i]/alpha_u  # velocity relaxation
        suu[i] = Au[i]*(p[i]-p[i+1]) + (1-alpha_u)*apu[i]*u[i]
        d[i] = Au[i]/apu[i]

    apu[0] = rho*(u[0]+u[1])/2*Ap[1] + rho * \
        (u[0]*Au[0]/Ap[0])*Ap[0]*0.5*(Au[0]/Ap[0])**2
    apu[0] = apu[0]/alpha_u  # velocity relaxation
    suu[0] = Au[0]*(p0 - p[1])+rho*(u[0]*Au[0]/Ap[0])*Ap[0] * \
        (Au[0]/Ap[0])*u[0]+(1-alpha_u)*apu[0]*u[0]
    d[0] = Au[0]/apu[0]

    awu[1] = rho*(u[0]+u[1])/2*Ap[1]
    apu[1] = (rho*(u[1]+u[2])/2*Ap[2])/alpha_u
    suu[1] = Au[1]*(p[1] - p[2]) + (1-alpha_u)*apu[1]*u[1]
    d[1] = Au[1]/apu[1]

    awu[nu-1] = rho*(u[nu-2]+u[nu-1])/2*Ap[nu-1]
    apu[nu-1] = (rho*u[nu-1]*Au[nu-1]*Au[nu-1]/Ap[n-1])/alpha_u
    suu[nu-1] = Au[nu-1]*(p[nu-1] - p[nu]) + (1-alpha_u)*apu[nu-1]*u[nu-1]
    d[nu-1] = Au[nu-1]/apu[nu-1]

    u = sol(awu, apu, aeu, suu)

    # pressure update
    for i in range(1,n-1):
        awp[i] = rho*d[i-1]*Au[i-1]
        aep[i] = rho*d[i]*Au[i]
        app[i] = awp[i]+aep[i]
        sup[i] = rho*u[i-1]*Au[i-1] - rho*u[i]*Au[i]
        sp.append(abs(sup[i]))

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
    residual.append(max(sp))

plt.figure(1)
plt.plot(mass)
plt.title("mass flow")
plt.figure(2)
plt.plot(residual)
plt.title("max residual")

print("mass flow:", mass[-1])
print("momentum residual:", max(sp))
```

> mass flow: 0.4441588782980728 
>
> momentum residual: 8.626240521558692e-06



![mass flow](/images/simple-algorithm/output_1_1.png)



![pressure](/images/simple-algorithm/output_1_2.png)

上面的代码使用一节迎风格式，采用压力修正方程中源项的绝对值的最大值作为判断收敛依据，采用100个压力节点得到收敛解：质量流量$0.444\ kg/s$，并得到了压力、速度的分布。

不妨试一试使用Van Leer格式，如下
$$
\phi_e = \phi_P + \frac{1}{2}\Psi(r)(\phi_E - \phi_P)
$$

$$
r = \frac{\phi_P - \phi_W}{\phi_E - \phi_p}
$$

$$
\Psi(r) = \frac{r+|r|}{1+r}
$$

![TVD](/images/simple-algorithm/tvd.png)

当$\Psi(r)$取零时就转化为了一阶迎风。理论上他是二阶精度的，使用下面的方法离散动量方程（忽略我这tedious的写法），需要特别小心边界上的处理方法，采用镜象法计算出$\Psi(r)$。其中同样使用了延迟修正方法，把TVD格式中除去迎风项都放到源项里，就可以用TDMA解方程了。

```python
def vanLeer(u,n):
    re = (u[n] - u[n-1])/(u[n+1] - u[n])
    rw = (u[n-1] - u[n-2])/(u[n] - u[n-1])
    phie = (re+abs(re))/(1+re)
    phiw = (rw+abs(rw))/(1+rw)
    return phie,phiw
    
for i in range(2, nu-1):
        fw = rho*(u[i-1]+u[i])/2*Ap[i]
        fe = rho*(u[i]+u[i+1])/2*Ap[i+1]
        phie,phiw = vanLeer(u,i)
        awu[i] = fw
        apu[i] = fe/alpha_u  # awu[i]+aeu[i]+ (F_e-F_w)
        suu[i] = Au[i]*(p[i]-p[i+1]) + fw*(0.5*phiw*(u[i] - u[i-1])) -\
             fe*(0.5*phie*(u[i+1] - u[i])) + (1-alpha_u)*apu[i]*u[i]
        d[i] = Au[i]/apu[i]

    re = 2*(u[0] - u[0]*Au[0]/Ap[0])/(u[1] - u[0])
    phie0 =  (re+abs(re))/(1+re)
    apu[0] = rho*(u[0]+u[1])/2*Ap[1] + rho * u[0]*Au[0]*0.5*(Au[0]/Ap[0])**2
    apu[0] = apu[0]/alpha_u  # velocity relaxation
    suu[0] = Au[0]*(p0 - p[1]) + rho*u[0]*u[0]*Au[0]*Au[0]/Ap[0] + (1-alpha_u)*apu[0]*u[0] -\
         0.5*rho*(u[0]+u[1])/2*Ap[1]*phie0*(u[1] - u[0])
    d[0] = Au[0]/apu[0]

    fw = rho*(u[0]+u[1])/2*Ap[1]
    fe = rho*(u[1]+u[2])/2*Ap[2]
    re = (u[1] - u[0])/(u[2] - u[1])
    phie =  (re+abs(re))/(1+re)
    awu[1] = fw
    apu[1] = fe/alpha_u
    suu[1] = Au[1]*(p[1] - p[2]) + (1-alpha_u)*apu[1]*u[1] + fw*(0.5*phie0*(u[1] - u[0])) -\
             fe*(0.5*phie*(u[2] - u[1]))
    d[1] = Au[1]/apu[1]

    fw = rho*(u[nu-1]+u[nu-2])/2*Ap[n-2]
    rw = (u[nu-2] - u[nu-3])/(u[nu-1] - u[nu-2])
    phiw = (rw+abs(rw))/(1+rw)
    awu[nu-1] = rho*(u[nu-2]+u[nu-1])/2*Ap[nu-1]
    apu[nu-1] = (rho*u[nu-1]*Au[nu-1]*Au[nu-1]/Ap[n-1])/alpha_u
    suu[nu-1] = Au[nu-1]*(p[nu-1] - p_exit) + (1-alpha_u)*apu[nu-1]*u[nu-1] + 0.5*fw*phiw*(u[nu-1] - u[nu-2])
    d[nu-1] = Au[nu-1]/apu[nu-1]
```

使用50个压力节点，得到解

> mass flow: 0.4466919233987552
>
> momentum residual: 8.999975353141121e-06

嗯，果然结果更接近解析解了，只用了一半的节点数，就达到了更高的精度。

#### 3.参考文献

- [1]  陶文铨, 数值传热学
- [2] H.K.Versteeg, An introduction to computational fluid dynamics : the finite volume method
- [3] J.H.Ferziger, M.Peric, Computational method for fluid dynamics
