# Ferziger and Blazek on LES


Some crap about turbulence and LES

<!--more-->

> Versteeg, H. K. and W. Malalasekera (2007). An introduction to computational fluid dynamics: the finite volume method, Pearson education.

雷诺数$Re$代表了流体中惯性力（与对流效应相关）和粘性力之间的比值。实验观察到，在某一临界雷诺数 ${Re}_{crit}$ 之下，流动是平稳的，流体的层与层间是有序的滑移。如果边界条件不随时间发生改变，流动是稳态的，此时称为**层流**。

当雷诺数高于${Re}_{crit}$值时，就会发生一系列复杂的事件，最终导致流动特性的根本变化。 在最终状态下，流动行为是随机且混乱的。 即使在恒定施加的边界条件下，运动本质上也变得不稳定。 速度和所有其他流动特性以随机和混乱的方式变化，此时流动称为**湍流**。

即使在速度和压力脉动仅存一维或者二维空间的流动中，湍流带来的扰动总是具有三维的空间特征。此外，湍流中往往还存在着涡，这些涡的尺度范围非常广，最大的涡通过称为**涡伸展**（vortex stretching）的过程与主流相互作用并从中提取能量。 剪切流中平均速度梯度的存在使旋转的涡变形，从而使得适当的涡被拉伸，因为它一端被迫比另一端移动得快。

最大的涡的特征速度$\vartheta$和特征长度$\ell$具有与主流速度和长度相当的量级，可以定义一个“大涡”雷诺数${Re}_{\ell}=\vartheta\ell/\nu$，这个数值很大并将接近主流雷诺数，表明大涡中起主导作用的是惯性力，粘性力基本可以忽略。因此，这些大涡具有良好的无粘性，并且在涡伸展过程中其动量守恒，在主流做功作用下，其速度增大，截面积减小。小一些的涡虽然也受到主流影响，但较大程度上受到大涡的影响，因此能量从大涡一级级地传送到最小的涡上，这被称为**能量级联**（energy cascade）。湍流蕴含的能量也因此包含很宽的频率（$f$）或波数（$\kappa=2\pi f/U$）范围。湍流能谱如下图

{{< figure src="/images/Ferziger-on-LES/spectrum.png" title="spectrum">}}

图中，谱能量$E(\kappa)$是波数的函数，单位$m^3/s^2$，它代表单位质量，单位波数的扰动动能。图中显示了最大的涡具有最大的动能，随着波数增大$E(\kappa)$迅速减小，因此最小的涡具有最低的能量。

湍流中最小尺寸的运动大概是0.1到0.01mm的量级，频率10kHz，最小涡的雷诺数$R e_{\eta}=v \eta / \nu=1$。

> Ferziger, J. H., et al. (2020). Computational methods for fluid dynamics, Springer.

## 1. 大涡模拟方程

湍流具有很宽的空间和时间尺度；在流体中的涡的尺度如左图所示，右图展示了一个典型的在流场中一点处速度分量的历史演变，这其中扰动的尺度是很明显的。

{{< figure src="/images/Ferziger-on-LES/vortex.png" title="vortex">}}

大尺度的涡往往比小尺度的涡具有更大的能量；它们的大小和强度使得它们称为流场中守恒变量最强有力的传播者，而小的涡则更弱，并对守恒变量提供更小的输运作用。因此大涡模拟的基本思路就是对这些影响更大的大涡进行更精确的模拟。LES已经成为大气科学的主要工具，用于研究云层、降水、污染物质运输和山谷中的气流。在地球大气边界层中，LES所处理的物理过程包括浮力、旋转、夹带、凝结以及与粗糙地面和海洋表面的相互作用。LES甚至也被用于设计高压燃气透平并且在空气声学模拟中占据主导地位。

要定义需要精确计算的量很重要。我们需要一个仅包含大尺度涡的速度场，而这个速度场最好是由原速度场滤波而得来的。在这种方法中，大的或被分辨率的场，即被模拟的场，本质上是完整的场的局部平均值。我们将使用一维符号;把它推广到三维是很简单的。滤波后速度定义为
$$
\bar{u}_i(x) = \int G(x,x^{'})u_i(x^{'})dx^{'}
$$
其中$G(x,x^{'})$为滤波器核函数，这是一个局部函数。在LES中应用的核函数包括一个高斯滤波器，一个盒式滤波器（box filter）和一个截止滤波器（将波数在截止值之上的所有傅里叶系数消除掉）。每一个滤波器都有一个与之对应的长度单位$\Delta$。一般认为尺度大于$\Delta$的涡是大涡，小于$\Delta$的是小涡，对这些小涡的计算需要引入模型。

当对不可压缩流的N-S方程加以滤波，就会得到一组与URANS（unsteady RANS）方程组相似的方程：
$$
\frac{\partial \rho\bar{u}_i}{\partial t}+\frac{\partial(\rho \overline{u_i u_j})}{\partial x_j}=-\frac{\partial \bar{p}}{\partial x_i}+\frac{\partial}{\partial x_j}\left[\mu(\frac{\partial \bar{u}_i}{\partial x_j}+\frac{\partial \bar{u}_j}{\partial x_i})\right]
$$

由于连续性方程是线性的，滤波后不改变：
$$
\frac{\partial(\rho \bar{u}_i)}{\partial x_i} = 0
$$

由于很重要的一点：
$$
\overline{u_iu_j}\neq \bar{u}_i\bar{u}_j
$$

而该不等式左侧难以计算，因此在近似计算时由该不等式引起的误差表示为
$$
\tau^s_{ij} = -\rho(\overline{u_iu_j}-\bar{u}_{i}\bar{u}_{j}
$$
因此LES计算中的动量方程为
$$
\frac{\partial \rho\bar{u}_i}{\partial t}+\frac{\partial(\rho \bar{u}_i\bar{u}_j)}{\partial x_j}=\frac{\partial \tau^s_{ij}}{\partial x_j}-\frac{\partial \bar{p}}{\partial x_i}+\frac{\partial}{\partial x_j}\left[\mu(\frac{\partial \bar{u}_i}{\partial x_j}+\frac{\partial \bar{u}_j}{\partial x_i})\right]
$$
在LES中，$\tau^s_{ij}$被称为亚网格尺度雷诺应力（subgrid-scale Reynolds stress）。“应力”源自于对于这一项的处理方式而非物理本质，它实际上是由小尺度或未解析尺度的作用所引起的大尺度动量通量。“亚网格尺度”这个名字也具有误导性，因为滤波窗口宽度，$\Delta$，除了明显需要满足的条件$\Delta \geq h$外不必与网格大小$h$有关。当今，用来描述亚网格尺度雷诺应力的模型分为subgrid-scale(SGS)和subfilter-scale(SFS)模型。

亚网格尺度雷诺应力包含小尺度场的局部平均值，因此它的模型应该基于局部速度场，或者可能基于局部流体运动的历史。后者可以通过求解偏微分方程模型来实现，从而获得确定亚网格尺度雷诺应力所需的参数。

## 2. 滤波[^3] 

滤波，一种平均化操作。在DNS中，速度场$\boldsymbol U(\vec{x},t)$需要在非常细密的网格上（直到Kolmogorov 尺度$\eta$）被解析，在LES中，则需要对速度场进行低通滤波，使得滤波后的速度场$\overline{\boldsymbol U}{(\vec{x},t)}$在较粗的网格上得以被解析。

滤波操作表示为
$$
\overline{\boldsymbol U}{(\vec{x},t)} = \int_{\Omega}G(\vec{r}, t)\boldsymbol{U}(\vec{x}-\vec{r},t)d\vec{r}
$$
其中$\Omega$表示整个流域，滤波函数满足归一化条件
$$
\int_{\Omega}G(\vec{r},\vec{x})d\vec{r} = 1
$$
在最简单的情况下，滤波函数是各向均匀的，即和$\vec{x}$无关。

假设$U(x)$为定义在$-\infty < x < \infty$ 上的标量函数，暂且称之为速度场，很容易可以将其扩展到三维。使用各项均匀滤波函数$G(r)$，滤波后的速度场为
$$
\overline{U}(x) = \int_{-\infty}^{\infty}G(r)U(x-r)dr
$$

以盒式滤波器（Box Filter）为例，$\overline{U}(x)$就是在$x-1/2\Delta < x^{\prime} < x+1/2\Delta$上对于$U(x^\prime)$进行平均得到的。高斯滤波函数选取为
$$
G(r) = \sqrt{\frac{6}{\pi \Delta^2}} exp({-\frac{6r^2}{\Delta^2}})
$$
即均值为0，方差为$\sigma^2 = 1/12\Delta^2$的高斯分布。$\sigma^2$的选取是根据其二阶矩$\int_{-\infty}^{\infty}r^2G(r)dr = 1/12\Delta^2$，保持与盒式滤波器一致而确定的。

```matlab
% 功能：对一维信号的高斯滤波，头尾r/2的信号不进行滤波
% r     :高斯滤波模板数（关联点数）
% sigma :标准差
% y     :需要进行高斯滤波的序列
% grid_scale: 网格分辨率
% delta :滤波尺度
grid_scale = 0.5;
x = 1:grid_scale:50;
y = x + rand(1,length(x))*10;

delta = 3; 
r        = 5; 
sigma    = sqrt(delta^2/12);
y_filted = Gaussianfilter(r, sigma, y);
fluc = y - y_filted;
% 作图对比
plot(x, y, x, y_filted, x, fluc);
title('高斯滤波');
legend('滤波前','滤波后','扰动')

function y_filted = Gaussianfilter(r, sigma, y)

% 生成一维高斯滤波模板
GaussTemp = ones(1,r*2-1);
for i=1 : r*2-1
    GaussTemp(i) = exp(-(i-r)^2/(2*sigma^2))/(sigma*sqrt(2*pi));
end

% 高斯滤波
y_filted = y;
for i = r : length(y)-r+1
    y_filted(i) = y(i-r+1 : i+r-1)*GaussTemp';
end
end
```

{{< figure src="/images/Ferziger-on-LES/gaussian.png" title="gaussian filter on scalar function">}}

## 3.  Smagorinsky与相关模型

最早和最常用的次网格尺度模型是由Smagorinsky[^1]提出的。这是一个涡-粘模型。所有这些模型都是基于这样一种假设，即亚网格雷诺应力的主要影响是输运作用和耗散的增加。由于这些现象是由于层流中的粘性造成的，通过合理地假设，一个可行的模型可能是
$$
\tau^s_{ij}-\frac{1}{3}\tau^s_{kk}\delta_{ij} = \mu_t\left(\frac{\partial \bar{u}_i}{\partial x_j}+\frac{\partial \bar{u}_j}{\partial x_i}\right)=2\mu_t\bar{S}_{ij}
$$
其中$\mu_t$表示涡粘度，$\bar{S}_{ij}$为大尺度场的应变率张量。Wyngaard[^2]展示了Lilly是如何在1967年使用$\tau^s_{ij}$演化方程推导出这个模型的。类似的模型也经常用于RANS方程。

亚网格尺度涡粘度可以由量纲理论推导出来，其形式为：
$$
\mu_t = {C_s}^2\rho\Delta^2|\bar{S}|
$$
其中$C_s$为无量纲模型参数，$\Delta$为滤波尺度，而 $|\bar{S}| =({\bar{S}_{ij}}^2)^{1/2}$

涡流粘度的这种形式可以以多种方式得到。 理论提供了参数的估计。 这些方法中的大多数仅适用于各向同性湍流，它们都得出$C_s\approx 0.2$。 不幸的是，$C_s$不是恒定的。 它可能是雷诺数和/或其他无量纲参数的函数，并可能在不同的工况中采用不同的值。

Smagorinsky模型尽管较为成功，但是仍然存在问题，它的使用也在下降，取而代之的是下面描述的更复杂的模型。虽然不推荐使用，但它仍然被使用。如果将该模型用于管道流，需要做多处修改。流体中$C_s$参数的值必须从0.2降低到大约0.065，这将涡粘性降低了几乎一个数量级。所有剪切流都需要这种幅度的变化。在靠近通道表面的区域，该值必须进一步减小（在不存在湍流的壁面处的涡粘度应该为0）。一种成功的方法是借用长期以来在RANS模型中用于降低近壁面涡粘性的van Driest阻尼:
$$
C_s = C_{s0}(1-e^{-n^+/A^+})^2
$$
其中$n^+$为无量纲壁面距离（$n^+=nu_{\tau}/\nu$，其中$u_{\tau}$为剪切速度，$u_{\tau} = \sqrt{\tau_{wall}/\rho}$，$\tau_{wall}$为壁面处切应力）。$A^+$为常数并通常取为25左右。尽管这个修改产生了期望的结果，但是很难在LES的视角下证明它是正确的。

更进一步的问题是，在靠近壁面的地方，流动的结构是各向异性的。这产生了低速和高速流体区域(条纹);在展向和法线方向上，它们都是大约1000个粘性单元长，30-50个粘性单元宽。为了解析这些条纹，需要高度各向异性的网格，并且SGS模型中滤波尺度$\Delta$的选择并不明确。通常选择为$(\Delta_1\Delta_2\Delta_3)^{1/3}$或者$(\Delta_1^2+\Delta_2^2+\Delta_3^2)^{1/2}$，其中$\Delta_i$为$i$方向的网格尺度。

因此，Smagorinsky模型存在许多困难。 如果我们希望模拟更复杂和/或更高的雷诺数流，那么就需要拥有更准确的模型。事实上，基于DNS数据得出的结果的详细测试表明，Smagorinsky模型在表示亚网格尺度应力的细节方面相当差。特别地，涡粘模型强制将$\tau^s_{ij}$与$\bar{S}_{ij}$关联起来，而*这种关联在现实上是不存在的*。

## 4. 动态(动力)模型

小尺度的未解析的涡在很多方面表现出和解析尺度下的小涡的相似性，基于这个想法提出另一种亚格子尺度模型，即尺度相似模型。其主要想法是强调解析尺度中的较大涡（尺度小于$\Delta$）和解析尺度中的较小涡（尺度大于$\Delta$）之间的相互作用。基于这个想法有模型
$$
\tau_{i j}^{\mathrm{s}}=-\rho\left(\overline{\bar{u}_{i} \bar{u}_{j}}-\overline{\bar{u}}_{i} \overline{\bar{u}}_{j}\right)
$$
其中双横线表示一个量被滤波两次，这个模型和实际的SGS雷诺应力非常相关，但是却几乎不能耗散任何能量因此也就无法作为一个独立的SGS模型。这一模型中从大尺度到小尺度的能量传递（正级联 forward scatter）和从最小的解析尺度到最大尺度的能量传递（反级联 back scatter）几乎相等，这一点很有用。为了纠正耗散过小的问题，可以结合Smagorinsky模型来构建一个混合模型。

在动态模型中，涡粘度表达式中的模型常数被替换为一个在时间和空间中演化的参数$C_d$
$$
\mu_t = C_d(\vec{r},t)\rho \Delta^2|\bar{S}|
$$
Germano提出使用“试滤波（test filter）” $\hat{\Delta}$ 来进行二次滤波，试滤波的尺度需要大于$\Delta$，通常$\hat{\Delta} = 2\Delta$，二次滤波后引入*亚试尺度应力(subtest-scale stress)* $\tau_{ij}^{ST}$
$$
\tau_{i j}^{S T}=-\rho(\widehat{\overline{u_{i} u_{j}}}-\hat{\bar{u}}_{i} \hat{\bar{u}}_{j})
$$
并且有
$$
\tau_{ij}^{ST} - \frac{1}{3}\tau_{ij}^{ST}\delta_{ij} = 2C_s^2\rho \hat{\Delta}|\hat{\bar{S}}|\hat{\bar{S}}_{ij}
$$
而在原滤波尺度下
$$
\tau_{i j}^{\mathrm{s}}=-\rho\left(\overline{u_{i} u_{j}}-\bar{u}_{i} \bar{u}_{j}\right)\\\\ \tau_{i j}^{\mathrm{s}}-\frac{1}{3} \tau_{k k}^{\mathrm{s}} \delta_{i j}=2 C_{S}^{2} \rho \Delta^{2}|\bar{S}| \bar{S}_{i j}
$$
$\tau_{ij}^{ST}$与亚格子尺度应力$\tau_{ij}^s$ 的关系表示为 Germano Identidy:
$$
\mathcal{L}_{i j}=\tau_{ij}^{ST} - \tau_{ij}^{s} = \rho\left(\widehat{\bar{u}_{i} \bar{u}_{j}}-\hat{\bar{u}}_{i} \hat{\bar{u}}_{j}\right)=2 C_{S}^{2} \rho\left(\hat{\Delta}|\hat{\bar{S}}| \hat{\bar{S}}_{i j}-\Delta^{2} \widehat{|\bar{S}| \bar{S}_{i j}}\right)+\frac{1}{3} \delta_{i j} \mathcal{L}_{k k}
$$
其中$\mathcal{L}_{kk}$包含各向同性项并且只有$C_s^2$是未知的。类比SGS应力，得到
$$
\mathcal{L}_{ij}-\frac{1}{3}\delta_{ij}\mathcal{L}_{ij} = 2C_d M_{ij}
$$
其中
$$
M_{ij} =\left(\hat{\Delta}|\hat{\bar{S}}| \hat{\bar{S}}_{i j}-\Delta^{2} \widehat{|\bar{S}| \bar{S}_{i j}}\right)
$$
注意到Germano Iddentity实际上对于一个未知数$C_s^2$给出了六个方程（对于每一个$S_{ij}$分量），因此方程是超定的。为了解决这一问题，Lilly使用最小二乘法给出了$C_d(\vec{r},t)$的表达式：
$$
C_d(\vec{r},t) = \frac{1}{2}\rho\frac{\mathcal{L}_{ij}M_{ij}}{M_{mn}M_{mn}}
$$
在实际使用中，常对分子分母取系综平均，即
$$
C_d(\vec{r},t) = \frac{1}{2}\rho\frac{\langle \mathcal{L}_{ij}M_{ij}\rangle}{\langle M_{mn}M_{mn}\rangle}
$$

动态模型的优点：

- 在剪切流中Smagorinsky系数要远小于在均匀各向同性湍流中的取值，动态模型可以自动做到这一点
- 模型参数在接近壁面处需要被进一步减小，动态模型自动在近壁区以合适的方式减小参数值
- 动态模型适合于各向异性网格，滤波尺度选取不确定的问题在一定程度上被动态过程修正。

存在的问题：

- 动态过程计算得到的模型参数是随时间和空间坐标快速变化的量，因此参数可能会取得非常大的正/负值。尽管负的涡粘度可以被视为是能量反级联的表现，但是过大的负值或者持续过长的时间有可能会引起数值不稳定。

## 5. LES算例

### Surface mounted cube

{{< figure src="/images/Ferziger-on-LES/surface mounted cube.png" title="surface mounted cube">}}

入口为充分发展湍流，出口边界条件为
$$
\frac{\partial \phi}{\partial t}+U\frac{\partial \phi}{\partial n}=0
$$
壁面为无滑移边界条件，侧面为周期性边界。使用$240\times 128\times128$的网格，二阶精度。

Ferziger的计算结果：

{{< figure src="/images/Ferziger-on-LES/result.png" title="result">}}

### Flow over a sphere in sub and supercritical Re

毕业论文中，在亚临界雷诺数$Re = 10^4$和超临界雷诺数$Re =10^6$下的球体绕流大涡模拟，表现出明显的湍流边界层再附着和尾迹区宽度的明显改变以及阻力系数的显著改变。

{{< figure src="/images/Ferziger-on-LES/sub_critic.png" title="subcritical">}}

{{< figure src="/images/Ferziger-on-LES/super_critic.png" title="supercritical">}}

{{< figure src="/images/Ferziger-on-LES/wallstress-1.png" title="subcritical wall stress">}}

{{< figure src="/images/Ferziger-on-LES/wallstress-2.png" title="supercritical wall stress">}}

[^1]:Smagorinsky, J. (1963). "General circulation experiments with the primitive equations: I. The basic experiment." Monthly weather review 91(3): 99-164.
[^2]:Wyngaard, J. C. (2010). Turbulence in the Atmosphere, Cambridge University Press.
[^3]:Pope. Turbulent Flows
