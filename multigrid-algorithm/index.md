# Multigrid Method


多重网格简单示例

<!--more-->

### 基本理论

考虑这样的线性方程组
$$
\mathbf{A}\vec{x} = \vec{b}
$$
当使用标准迭代求解器（Gauss-Seidel，Jacobi，SOR）时，观察到收敛速度趋于“停顿”，即在几次迭代后无法有效减少误差，细化网格后，问题更加突出。 实际上，研究表明，收敛速度是误差场频率的函数，即误差从节点到节点的梯度。 如果误差以高频模式分布，则收敛速度很快。 然而，在最初的几次迭代之后，误差场被平滑（变为低频），使得收敛速度变差。

一个显而易见的想法是：既然平滑后的误差在细网格上变为了低频，那么这个误差可以在粗网格上得到有效消除，因为细网格上的高频可能是粗网格上的高频！

迭代得到的中间值设为$\vec{y}$，则
$$
\mathbf{A}\vec{y} = \vec{b} - \vec{r}
$$
$\vec{r}$为残差向量，如果用$\vec{e} = \vec{x} - \vec{y}$表示中间值和解之间的距离，则
$$
\mathbf{A}\vec{e} = \vec{r}
$$
因此这个方程与原来要解的方程具有同样的系数矩阵，如果可以解出$\vec{e}$，那么便可修正中间值。为了有效地达成这一目的，把上面的方程放到粗网格中计算，以快速光滑残差。在这之后，将$\vec{e}$再投影回细网格来修正解，如此反复迭代。

多重网格方法常用的几个名词：

- `Agglomeration` ：所有的多重网格方法都需要在原始“精细”网格的基础上定义一系列粗网格。定义粗网格的过程涉及到所谓的聚集，即从原始网格中组合几个节点或控制量或系数，得到粗网格上的稀疏矩阵。
- `Restriction`：利用插值方法将误差向量投影到粗网格上
- `Prolongation`：Restriction的反过程

多重网格法的两种形式：

- Geometric multigrid：在几何多重网格中，生成了网格的层次结构。方程在每个级别上进行离散。 几何多重网格优于代数多重网格的优点是，对于非线性问题，前者应具有更好的性能，因为系统中的非线性通过重新离散化可以降低到粗糙的水平。
- Algebraic multigrid：该算法被称为代数多重网格方案，因为生成的粗糙层方程没有使用任何几何或在粗糙层上的重新离散。这样做的优点是不需要生成或存储粗级别的网格，也不需要在粗级别上计算通量或源项。这一特性使得AMG对于非结构化网格的使用尤为重要。

多重网格法的关键算子：

- Prolongation算子$\mathbf{P}$，将粗网格上的结果插值到细网格上

- Restriction算子 $\mathbf{R}$，$\mathbf{R} = \frac{1}{2}\mathbf{P^T}$
- 代数多重网格的粗网格矩阵：$\mathbf{A_{2h}} = \mathbf{R}\mathbf{A_h}\mathbf{P}$

### 算例

![problem](./images/multigrid-method/multigrid.png)

一个带内热源的导热问题
$$
k\frac{d^2T}{dx^2} + g = 0
$$
长度1m，截面积0.01$m^2$，$k = 5 W/m.K$，$q = 20 kW/m^3$

离散过程很简单，最后得到一个三对角矩阵，可以使用TDMA直接求解，在这里将使用多重网格进行求解。

- V型循环，三层网格，网格数依次为20、10、5
- 最后一层使用TDMA直接求解
- $\mathbf{P}$和$\mathbf{R}$都使用线性插值代替
- 粗网格矩阵可以重新离散得到或者使用插值

```python
import numpy as np
import matplotlib.pyplot as plt

k = 5
q = 20e3
a = 0.01
n = 20
n2 = 10
n3 = 5
tl = 100
tr = 500
l = 1.0


def restrict(vec):
    n = int(np.size(vec)/2)
    vec_coarse = np.zeros(n)
    for i in range(n):
        vec_coarse[i] = (vec[2*i]+vec[2*i+1])/2

    return vec_coarse


def prolong(vec):
    n = np.size(vec)
    vec_fine = np.ones(n*2)*vec[0]/2
    for i in range(1, n*2):
        if (i % 2 == 0):
            vec_fine[i] = (vec[int(i/2) - 1] + vec[int(i/2)])/2
        else:
            vec_fine[i] = vec[int(i/2)]

    return vec_fine


# iter stands for itereation times
def gsIter(aw, ap, ae, su, phi, itr):
    n = np.size(aw)
    phi0 = phi
    for i in range(itr):
        phi[0] = (ae[0]*phi0[1]+su[0])/ap[0]
        for i in range(n-1):
            phi[i] = (aw[i]*phi[i-1]+ae[i]*phi0[i+1]+su[i])/ap[i]

        phi[n-1] = (aw[n-1]*phi[n-2]+su[n-1])/ap[n-1]
        for i in range(n):
            phi0[i] = phi[i]

    return phi


def residual(su, aw, ap, ae, phi):
    n = np.size(su)
    res = np.zeros(n)
    for i in range(1, n-1):
        res[i] = su[i] - (-aw[i]*phi[i-1]+ap[i]*phi[i]-ae[i]*phi[i+1])
        res[0] = su[0] - ap[0]*phi[0] + ae[0]*phi[1]
        res[n-1] = su[n-1] - ap[n-1]*phi[n-1] + aw[n-1]*phi[n-2]
    return res


def getMatrix(num):
    aw = np.zeros(num)
    ap = np.zeros(num)
    ae = np.zeros(num)
    su = np.zeros(num)
    for i in range(1, num-1):
        aw[i] = k*a*num/l
        ae[i] = aw[i]
        su[i] = q*a*l/num
        ap[i] = aw[i]+ae[i]

    aw[0] = 0
    ae[num-1] = 0
    ae[0] = k*a*num/l
    aw[num-1] = k*a*num/l
    su[0] = q*a*l/n+2*k*a*num/l*tl
    su[num-1] = q*a*l/n+2*k*a*num/l*tr
    ap[0] = aw[0]+ae[0] + 2*k*a*num/l
    ap[num-1] = aw[num-1]+ae[num-1] + 2*k*a*num/l
    return aw, ap, ae, su


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


if __name__ == '__main__':
    # prepare matrix
    aw, ap, ae, su = getMatrix(n)
    aw2, ap2, ae2, su2 = getMatrix(int(n/2))
    aw3, ap3, ae3, su3 = getMatrix(int(n/4))
    phi = np.ones(n)*150
    phi = gsIter(aw, ap, ae, su, phi, 3)

    res = 1
    cnt = 0
    res_list = []
    while(res > 1e-6):
        # finest mesh
        phi = gsIter(aw, ap, ae, su, phi, 2)
        r = residual(su, aw, ap, ae, phi)
        res = np.mean(r)
        res_list.append(res)
        cnt += 1
        print("iter: ", cnt, " res:", res)
        # ------------- restriction -----------------
        r2 = restrict(r)
        e2 = gsIter(aw2, ap2, ae2, r2, np.zeros(int(n/2)), 10)
        r4 = residual(r2, aw2, ap2, ae2, e2)
        rr = restrict(r4)
        #e4 = gsIter(aw3,ap3,ae3,rr,np.zeros(int(n/4)),10)
        e4 = sol(aw3, ap3, ae3, rr)
        # ------------ prolongation -----------------
        eff = prolong(e4)
        e2 = e2+eff
        e2 = gsIter(aw2, ap2, ae2, r2, e2, 2)
        ef = prolong(e2)
        phi = phi + ef

    plt.yscale('log')
    plt.figure(1)
    plt.plot(res_list)
    plt.show()
```

```
iter:  1  res: 11.51138138623885
iter:  2  res: 5.799009148342661  
iter:  3  res: 3.156842353629591  
iter:  4  res: 1.7547835232202516 
iter:  5  res: 0.9784454646211088 
iter:  6  res: 0.5436277301623192 
iter:  7  res: 0.3001937730563299 
iter:  8  res: 0.1646310272677738 
iter:  9  res: 0.08963673469351363
iter:  10  res: 0.048440945656969346
iter:  11  res: 0.025991466883453995
iter:  12  res: 0.013859948257314158
iter:  13  res: 0.007354351924635694
iter:  14  res: 0.003887729670987028
iter:  15  res: 0.0020495022264498176
iter:  16  res: 0.0010782926067378184
iter:  17  res: 0.0005665134079634981
iter:  18  res: 0.0002973377378992836
iter:  19  res: 0.0001559499926912622
iter:  20  res: 8.175372336438613e-05
iter:  21  res: 4.284309166848743e-05
iter:  22  res: 2.2446576915058357e-05
iter:  23  res: 1.1758364624370188e-05
iter:  24  res: 6.158758171181944e-06
iter:  25  res: 3.2255524303081986e-06
iter:  26  res: 1.6892364982368235e-06
iter:  27  res: 8.846259248684873e-07
```

经过27次V循环，迭代收敛。

### 分析

![result](/images/multigrid-method/result.png)

可以看到，单纯的高斯赛德尔迭代需要近660步才能迭代到收敛，而采用多重网格后，只需要27个V循环，每个循环由2次细网格上的GS迭代和10次较粗网格上的GS迭代以及最粗网格的TDMA组成，大大减少了计算量。



### 参考

[1] [MIT courseware](https://ocw.mit.edu/courses/mathematics/18-086-mathematical-methods-for-engineers-ii-spring-2006/readings/am63.pdf)

[2] H.K.Versteeg, An introduction to computational fluid dynamics: the finite volume method

[3]  [CFD-online Multigrid methods](https://www.cfd-online.com/Wiki/Multigrid_methods)


