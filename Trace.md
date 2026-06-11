# 有关离散$\text{Gauss}$分布闭合轨迹

我们期望讨论在一条简单参数曲线上，一个质点（或人）试图沿参数曲线运动产生的真实轨迹是类似什么样的，以及如何理想地模拟这一系列离散点的分布. 所以提出一种构想，在每个时刻向期望方向按照$\text{Gauss}$分布前进，下面试图研究这种模拟方法，并分析结果.

## 1.瞬时期望方向

首先引入瞬时期望方向的概念.假设在一条简单参数曲线$\Gamma (t),(t\in D)$外(上)有一点$P$，对于我们的理想模型，需要确定一个单位方向，用于刻画在此刻试图回到$\Gamma (t)$并沿$t$增大方向前进的趋势.这一向量用$\boldsymbol{r}_\text{Exp}(P)$表示.  
首先描述距离$P$最近点的曲线上点，直接采用最短欧式距离
$$
D_\text{min}=\min \{r>0|\exists t\in D,r=\|\boldsymbol{r}_P-\boldsymbol{r}(t)\|\}
$$
此时即有
$$
t_\text{min}=\argmin_{t\in D} \|\boldsymbol{r}_P-\boldsymbol{r}(t)\|,Q=\Gamma(t_\text{min})
$$
所以回归方向向量
$$
\boldsymbol{n}=\boldsymbol{r}_\text{back}=\frac{\boldsymbol{r}_Q-\boldsymbol{r}_P}{\|\boldsymbol{r}_Q-\boldsymbol{r}_P\|}
$$
前进方向向量直接取切向量
$$
\boldsymbol{\tau}=\boldsymbol{r}_\text{forward}=\frac{\Gamma'(t_\text{min})}{\|\Gamma'(t_\text{Min})\|}
$$
综合来看，我们给出线性优化
$$
\boldsymbol{r}_\text{Exp}(\lambda,P)=\frac{\boldsymbol{\tau}+\lambda\boldsymbol{n}}{\|\boldsymbol{\tau}+\lambda\boldsymbol{n}\|}
$$
这样得出的$\boldsymbol{r}_\text{Exp}$是线性组合，$\lambda$表示“回归”这一运动的权重.类似于此可以给出另一种更为直观的$\boldsymbol{r}_\text{Exp}$表示.对于离散的情况，简单地选取切向量方向其实并不合适，更合适的考虑是直接选择割线，所以对$\boldsymbol{\tau}$的选择需要改进.从$Q$点考虑前进长度$L$到达$Q_+$，那么在离散情况下，向$Q_+$前进显然更适合描述.所以定义
$$
\exist T,s.t. \int_{t_\text{min}}^{t_\text{min}+T}\|\Gamma'(t)\|\mathrm{d} t=L,\boldsymbol{r}_{Q_+}:=\boldsymbol{r}(t_\text{min}+T)
$$
$$
\boldsymbol{r}_\text{Exp}(P,L)=\frac{\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P}{\|\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P\|}      \tag {*}
$$
实际上局部来看，$d_\text{min}=\|\boldsymbol{r}_Q-\boldsymbol{r}_P\|$，则该式近似为
$$
\boldsymbol{r}_{\text{Exp}}\approx\frac{L\boldsymbol{\tau}+d_\text{min}\boldsymbol{n}}{\sqrt{L^2+d_{\text{min}}^2}}$$
也即它相当于上一个版本取
$$
\lambda\approx\frac{d}{L}
$$
的特殊情况.综合来看，我们使用$(*)$标记的形式.  

## 2.$\text{Gauss}$分布路径跟随
接下来考虑在$P$点处的随机游走，用角度按$\text{Gauss}$分布作近似处理，期望值按照选定的
$$
\theta_{\text{Exp}}=\arg \boldsymbol{r}_\text{Exp}(P,L)=\arg \frac{\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P}{\|\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P\|}
$$
所以在$P$处选取前进方向
$$
\theta=\theta_\text{Exp}+\varepsilon,\varepsilon\sim N(0,\sigma^2),\boldsymbol{r}_{P_{k+1}}=\boldsymbol{r}_{P_k}+v\Delta t \cdot\boldsymbol{r}_k
$$
其中$\boldsymbol{r}_k=(\cos \theta,\sin \theta)$，这就是带有随机误差的运动.综合来看，先考虑直线运动的例子.设直线$L:x=0$，初始从$P_0=O$开始运动，给定一组$\text{Gauss}$分布的随机变量列
$$
\{\varepsilon_{n,\sigma}\},\forall i,\varepsilon_{n,\sigma}\sim N(0,\sigma^2)
$$
$L$参数选取为$y$坐标，即$y=t$.给定$P(x,y)$点时，有$t_\text{min}=y$，显然$T=L$，$\boldsymbol{Q_+}=(0,y+L)$，所以
$$
\boldsymbol{r}_\text{Exp}(P,L)=\frac{\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P}{\|\boldsymbol{r}_{Q_+}-\boldsymbol{r}_P\|}=\frac{1}{\sqrt{x^2+L^2}}(-x,L)
$$
所以
$$
\theta_\text{Exp}=\operatorname{atan2}(L,-x)=\frac{\pi}{2}+\arctan{\frac{x}{L}},\theta=\frac{\pi}{2}+\arctan{\frac{x}{L}}+\varepsilon
$$
所以有递推
$$
\Delta(x,y)=v\Delta t\left(\cos{\left(\frac{\pi}{2}+\arctan \frac{x}{L}+\varepsilon\right)},\sin \left(\frac{\pi}{2}+\arctan \frac{x}{L}+\varepsilon\right)\right)
$$
将$\sin$和$\cos$展开，得到
$$
\Delta (x,y)=-\frac{v\Delta t}{\sqrt{x^2+L^2}}\left(x\cos \varepsilon+L\sin \varepsilon,x\sin \varepsilon-L\cos\varepsilon\right)
$$
可以得到递推
$$
x_{k+1}-x_k=-\frac{v\Delta t}{\sqrt{x_k^2+L^2}}(x_k\cos \varepsilon_k+L\sin \varepsilon_k)
$$
$$
y_{k+1}-y_k=-\frac{v\Delta t}{\sqrt{x_k^2+L^2}}(x_k\sin \varepsilon_k-L\cos \varepsilon_k)
$$
再譬如取标准椭圆
$$
\frac{x^2}{a^2}+\frac{y^2}{b^2}=1
$$
参数式写为
$$
(x,y)=(a\cos t,b\sin t),t\in [0,2\pi)
$$
先解决$L$的问题，考虑弧长
$$
\int _{t_\text{min}}^{t_\text{min}+T}\sqrt{b^2\cos ^2 t+a^2 \sin^2 t}\space \mathrm{d}t=L
$$
记$t_\text{min}=t_0$，$a^2-b^2=\Delta$，则$f(t)=\sqrt{b^2+\Delta \sin^2 t}$，作$\text{Taylor}$展开
$$
f(t_0+u)=f(t_0)+f'(t_0)u+\frac{1}{2}f''(t_0)u^2+O(u^3)
$$
$$
f'(t)=\frac{\Delta \sin 2t}{2f(t)},f''(t)=\frac{\Delta \cos 2t}{f(t)}-\frac{\Delta ^2\sin^2 2t}{4f^3(t)}
$$
所以积分
$$
\int_{t_0}^{t_0+T}f(t)\mathrm{d}t=\int_0^T f(t_0+u)\mathrm{d}u
$$
代入$\text{Taylor}$展开得
$$
\int_0^Tf(t_0)+\frac{\Delta \sin 2t_0}{2f(t_0)}u+\frac{1}{2}\left(\frac{\Delta \cos 2t_0}{f(t_0)}-\frac{\Delta ^2\sin^2 2t_0}{4f^3(t_0)}\right)u^2+O(u^3)\space \mathrm{d}u
$$
$$
=f(t_0)T+\frac{\Delta \sin 2t_0}{4f(t_0)}T^2+\frac{1}{6}\left(\frac{\Delta \cos 2t_0}{f(t_0)}-\frac{\Delta ^2\sin^2 2t_0}{4f^3(t_0)}\right)T^3+O(T^4)
$$
一阶近似下，可以取$L\approx f(t_0)T$，即
$$
T=\frac{L}{f(t_0)}=\frac{L}{\sqrt{b^2\cos ^2t_\text{min}+a^2\sin^2t_\text{min}}}
$$
二阶近似下，可以取$\displaystyle L\approx f(t_0)T+\frac{\Delta \sin 2t_0}{4f(t_0)}T^2$，此时直接给出近似反解结果
$$
T\approx \frac{L}{f(t_0)}-\frac{\Delta \sin 2t_0}{4f^4(t_0)}L^2
$$
下面是最近点$t_\text{min}$的选择问题，显然有
$$
t_\text{min}=\argmin_{t\in [0,2\pi)}D^2(t)=\argmin_{t\in [0,2\pi)}\left((a\cos t-p_x)^2+(b\sin t-p_y)^2\right)

$$
对$t$显式求导
$$
\frac{\mathrm{d}}{\mathrm{d}t}\left((a\cos t-p_x)^2+(b\sin t-p_y)^2\right)=0
$$
整理得到
$$
ap_x\sin t-bp_y\cos t+(b^2-a^2)\sin t\cos t=0
$$
一般来说也可以用$\text{Lagrange}$乘数处理，代入$\lambda$后得到方程（默认乘子使用$\mathcal{L}=f+\lambda \varphi$）
$$
\mathcal{L}=(x-p_x)^2+(y-p_y)^2+\lambda(\frac{x^2}{a^2}+\frac{y^2}{b^2}-1),\nabla_{x,y,\lambda} \mathcal{L}=0
$$
通过消元可以得到稳定的一维方程
$$
\frac{a^2 p_x^2}{(a^2+\lambda)^2}+\frac{b^2p_y^2}{(b^2+\lambda)^2}=1
$$
解得$\lambda$后最近点坐标就是
$$
Q\left(\frac{a^2 p_x}{a^2+\lambda},\frac{b^2 p_y}{b^2+\lambda}\right),t_\text{min}=\operatorname{atan2}\left(\frac{b p_y}{b^2+\lambda},\frac{a p_x}{a^2+\lambda}\right)
$$
于是根据定义，可以简单地给出$\boldsymbol{r}_\text{Exp}$的数值.后面给出三组模拟测试的输出结果（分别为直线、椭圆一阶近似、椭圆二阶近似）.

## 3.进一步的修正

### 3.1 关于$L-T$计算的$\text{Taylor}$简化
我们考虑之前给出的$T-L$定义
$$
\exist T,s.t. \int_{t_\text{min}}^{t_\text{min}+T}\|\Gamma'(t)\|\mathrm{d} t=L,\boldsymbol{r}_{Q_+}:=\boldsymbol{r}(t_\text{min}+T)
$$
一般若$\|\Gamma'(t)\|$无初等原函数，不妨使用$\text{Taylor}$级数解决，记$\|\Gamma'(t)\|=\varphi(t)$将$\varphi(t)$在$t_0$处展开
$$
    \varphi(t_0+u)=\sum_{n=0}^\infty \frac{\varphi^{(n)}(t_0)}{n!}u^n=\sum_{n=0}^N\frac{\varphi^{(n)}(t_0)}{n!}u^n+O(u^{N+1})
$$
所以
$$
\int_{t_\text{min}}^{t_\text{min}+T}\|\Gamma'(t)\|\mathrm{d} t=\int_0^T \varphi(t_0+u)\mathrm{d}u=\sum_{n=0}^N\frac{\varphi^{(n)}(t_0)}{(n+1)!}T^{n+1}+O(T^{N+2})
$$
取切线近似时，则$\|\Gamma'(t_0)\|T=L$，$\displaystyle T=\frac{L}{\|\Gamma'(t_0)\|}$，这就退化到了第一种$\boldsymbol{r}_\text{Exp}$的形式. 所以可以考虑近似到二阶或更高阶，这里需要幂级数反演操作，根据幂级数展开唯一性可以证明反演的合法性.关于幂级数反演，这里的主要思想在于，现在的形式为
$$
L=L(T)=\sum C_i T^i 
$$
级数，而我们需要
$$
T=T(L)=\sum C_i'L^i
$$
级数形式，只要先假设下式成立，代入上式反解系数即可. 对一般的二阶或三阶，不加验证地直接给出前四阶结果
$$
T=\frac{L}{A_0}-\frac{A_1}{2A_0^3}L^2+\frac{3A_1^2-A_0A_2}{6A_0^5}L^3+\frac{-15A_1^3+10A_0A_1A_2-A_0^2A_3}{24A_0^7}L^4+O(L^5)
$$
其中$A_k=\varphi^{(k)}(t_0)$.  可见普通的幂级数反演是相当复杂的. 所以给出$\text{Lagrange}$反演公式：$L=f(u)$，$f(0)=0$，$f'(0)=a_1\ne 0$，则
$$
u=f^{-1}(L)=\sum_{n=1}^\infty  \frac{1}{n}[z^{n-1}]\left(\frac{z}{f(z)}\right)^n L^n
$$
其中$[z^{n-1}]F(z)$表示取形式幂级数$F(z)$中$z^{n-1}$的系数，也即这个系数实际上是可以由$\text{Taylor}$公式给出的，是一个机械的流程，这里就不给出了. 但总而言之，对于不大的$L$，我们还是根据前四阶（或更少）的反演级数计算就足够了，给出$\text{Lagrange}$反演实际上是在面对$\varphi(t)$不初等可积时，一些$\kappa(t_0)$大或$\|\Gamma'(t)\|$变化复杂的地方时至少提供了提高阶数（当然还有减小$L$）以确定$\boldsymbol{r}_\text{Exp}$的可能性.

### 3.2 $d_\text{min}$计算的迭代
在实际计算中，最近点$Q_k$的求解不宜在每一步都作为独立的全局优化问题处理. 一方面，若每次更新 $P_k$后都重新应用 $\text{Lagrange}$乘数法或全局搜索计算最近点，当步数$n$很大时计算成本较高；另一方面，对于空间中不同曲线段相互接近的情形，全局最近点可能落在与当前运动分支无关的参数区间内，从而导致路径跟踪意义上的错误. 因此，在模拟中更合理的做法是将最近点参数 $t_k$作为状态量进行递推跟踪. 具体而言，可由上一时刻的 $t_k$或其切向预测值给出$t_{k+1}$的初始估计，再在该初值附近的局部参数区间内进行$\text{Newton}$迭代或局部搜索，从而得到与上一时刻连续一致的最近点$Q_{k+1}=\Gamma(t_{k+1})$. 上面的椭圆例子已经体现出了反复计算对总计算量的贡献是极大的，但这还不是严重的错误，为了解决这一计算问题，采用$\text{Newton}$迭代，显然有初始值
$$
t_{k+1}^{(0)}=t_k+\frac{(\boldsymbol{r}_{P_{k+1}}-\boldsymbol{r}_{P_k})\cdot\Gamma'(t_k)}{\|\Gamma'(t_k)\|^2}
$$
假设寻找的最近点满足驻点方程
$$
(\boldsymbol{r}_{P_{k+1}}-\Gamma(t))\cdot \Gamma'(t)=0
$$
代入$\text{Newton}$迭代公式
$$
t^{(m+1)}=t^{(m)}-\frac{(\boldsymbol{r}_{P_{k+1}}-\Gamma(t^{(m)}))\cdot \Gamma'(t^{(m)})}{-\|\Gamma'\left(t^{(m)}\right)\|^2+(\boldsymbol{r}_{P_{k+1}}-\Gamma(t^{(m)}))\cdot \Gamma''(t^{(m)})}
$$
一般而言，对不那么病态的曲线，迭代次数$M$不会很大. 这在数值计算上也是方便的. 对于第二个曲线“跳跃”的问题，使用局部投影的方法解决，定义
$$
I_k:=[t_{k+1}^{(0)}-\delta,t_{k+1}^{(0)}+\delta]
$$
$$
t_{k+1}=\argmin_{t\in I_k}\|\boldsymbol{r}_{P_{k+1}}-\Gamma(t)\|^2
$$
这与$\text{Newton}$迭代的流程存在一定的冲突：即寻找的驻点不一定在$I_k$之内. 因此，我们选取保守一点的策略，即：  
实际计算时，在预测参数附近取局部区间$I_k$先用$\text{Newton}$迭代求$h'(t)=0$的内部驻点；若驻点存在且位于$I_k$内，则将其与区间端点共同作为候选点，比较$h(t)=\|\boldsymbol{r}_{P_{k+1}}-\Gamma(t)  \|^2$
的大小并取最小者；若$\text{Newton}$迭代未收敛或驻点不在区间内，则直接取两个端点中距离较小者. 这样既保留了局部最近点的意义，又避免了全局最近点跳转到其他曲线段的问题. 
### 3.3 转向惯性系数$\alpha$
对于真实的跑步前进模型，我们不能期望人的前视会立即反应在调整角度上，也即对于真实的转向，需要一个惯性的量度，称之为$\alpha$，$\alpha \in (0,1]$，用于刻画转向的惯性. 之前的模型，我们认为
$$
\theta_{k}=\theta_{\text{Exp},k}+\varepsilon_k\space,\space \theta_{\text{Exp},k}=\arg \boldsymbol{r}_{\text{Exp,k}}\space,\space \varepsilon_k\sim N(0,\sigma^2)
$$
现在我们维护$\theta_k$作为一个递推的变量，即保留惯性
$$
\theta_{k+1}=\theta_k+\alpha \operatorname{wrap}(\theta_{\text{Exp},k+1}-\theta_k)+\varepsilon_{k+1}
$$
其中$\operatorname{wrap}(x)$表示将$x$调整到$(-\pi ,\pi]$间. 这样$\alpha \in (0,1]$，当$\alpha$接近$0$时，转向惯性极大，即基本不会转向；相反地，$\alpha$接近$1$时，转向惯性极小，基本上立即转向新的$\boldsymbol{r}_\text{Exp}$方向.

### 3.4 $v$的递推处理&间隔采样

既然3.3提到“惯性”，我们也不难想到实际物理运动中在标量$v$的性质上，也有惯性的重要作用，所以我们对$v$也采用递推处理的操作，并赋予随机噪声，表述如下
$$
v_{k+1}=v_k+\beta(v_{\text{Exp}}-v_k)+\varepsilon_{v,k+1}\space,\space \varepsilon_{v,k+1}\sim N(0,\sigma_v^2)
$$
其中$v_\text{Exp}$表示期望的平均速度，$\beta \in(0,1]$表示速度调整系数，$\beta$越大，则$v_{k+1}$越倾向于快速调整到$v_\text{Exp}$，当然一般而言，根据我们的真实感受，对于一个对速度不太敏感的人而言，$\beta$的值应该是比较小的. 完整来看，我们的综合递推应该写作：（为了不造成歧义，对之前的$\text{Gauss}$随机变量加上了下标）
$$
\theta_{k+1}=\theta_k+\alpha \operatorname{wrap}(\theta_{\text{Exp},k+1}-\theta_k)+\varepsilon_{\theta,k+1}\space,\space \varepsilon_{\theta,k+1}\sim N(0,\sigma_\theta^2)
$$
$$
v_{k+1}=v_k+\beta(v_{\text{Exp}}-v_k)+\varepsilon_{v,k+1}\space,\space \varepsilon_{v,k+1}\sim N(0,\sigma_v^2)
$$
$$
\boldsymbol{r}_{P_{k+1}}=\boldsymbol{r}_{P_k}+v_{k+1}\Delta t(\cos \theta_{k+1},\sin \theta_{k+1})
$$
不过实际计算时也可以考虑加上正值限制
$$
v_{k+1}=\max(0,v_k+\beta(v_{\text{Exp}}-v_k)+\varepsilon_{v,k+1})
$$
当然，在$\text{Gauss}$效果够好的情况下，这也可以忽略. 另外，对于真实的“设备测量”，我们有必要区别运动模拟间隔和采样间隔. 最简单也最合理的策略，就是间隔一个固定的$N\Delta t$采样. 
$$
 \forall k,P_{Nk}\space\text{exists}\space,\space {P_{Nk}}\in\{P_{\text{scan},n}\}
$$

### 3.5 测量随机误差
这是最简单的修正，只需要考虑每一项测量值的随机误差即可，假定$\Delta x,\Delta y$都是$\text{Gauss}$分布的随机变量
$$
\Delta x,\Delta y \sim N(0,\sigma_{\text{measure}}^2)
$$
那么极坐标分布函数
$$
f(\rho,\theta)=\frac{1}{2\pi \sigma^2_\text{measure}}\exp\left(-\frac{x^2}{2\sigma_\text{measure}^2}\right)\cdot\exp\left(-\frac{y^2}{2\sigma_\text{measure}^2}\right)$$
$$
=\frac{1}{2\pi \sigma^2_\text{measure}}\exp\left(-\frac{\rho^2}{2\sigma_\text{measure}^2}\right)
$$
因为$\displaystyle \frac{\partial f(\rho,\theta)}{\partial \theta}=0$，所以分布保证圆对称性.  所以直接
$$
P_{\text{scan},n}\leftarrow P_{\text{scan},n}+(\Delta x,\Delta y)
$$
即可. 

