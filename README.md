# CZ段路明方案Note #

## 基本模型 ##

### 模型 ###

### ![CZ1a](/CZ1a.png) ###

考虑单个光子与腔内原子的相互作用。哈密顿算符如下：
$$
H=(\Delta-i\gamma_e/2)|e><e|+g(a|e><g|+a^{\dagger}|g><e|)+i\sqrt{\kappa/2\pi}\int_{-\omega_b}^{\omega_b}d\omega(a^{\dagger}b(\omega)-ab^{\dagger}(\omega))+\int_{-\omega_b}^{\omega_b}d\omega(\omega b^{\dagger}(\omega)b(\omega)) +(\Delta-i\gamma/2)a^{\dagger}a
$$
其中$\Delta$表示红失谐，$|s>$态就是图中的$|0>$态，$|g>$态就是图中的$|1>$态。

系数g是原子与腔的耦合系数，$\gamma_e$是原子在e态自发辐射导致的损耗因子，$\gamma$表示腔模中的光子逃逸到系统之外导致的损失，$\kappa$表示腔和波导的耦合系数。

$a$表示腔模，$b$表示波导上的模式。

$\omega_b$表示频域上计入考虑的范围

系统的态可以表示为：
$$
|\Psi>=(c_g |g>|1>_a+c_e|e>|0>_a)|0>_b+c_s|s>|1>_a|0>_b+|s>|0>_a \int_{-\omega_b}^{\omega_b}d\omega c_{\omega,s}b^{\dagger}(\omega)|0>+|g>|0>_a \int_{-\omega_b}^{\omega_b}d\omega c_{\omega,g}b^{\dagger}(\omega)|0>
$$
$c_g$表示光子在腔模上且原子在g态的机率幅；

$c_e$表示光子被吸收原子处于e态的机率幅；

$c_s$表示光子在腔模上且原子在s态的机率幅；

$c_{\omega,s}$表示原子在s态上且光子在波导$\omega$​模上的机率幅；

$c_{\omega,g}$表示原子在g态上且光子在波导$\omega$模上的机率幅。

### 解析计算 ###

$H|\Psi>=i\hbar\frac{d}{dt}|\Psi>$

解得：
$$
\dot{c_e}=(-i\Delta-\gamma_e/2)c_e-igc_g\\
\dot{c_g}=(-i\Delta-\gamma/2)c_g-igc_e+\sqrt{\kappa/2\pi}\int_{-\omega_b}^{\omega_b}d\omega c_{\omega,g}\\
\dot{c_s}=(-i\Delta-\gamma/2)c_s+\sqrt{\kappa/2\pi}\int_{-\omega_b}^{\omega_b}d\omega c_{\omega,s}\\
\dot{c_{\omega,g}}=-i\omega c_{\omega,g}-\sqrt{\kappa/2\pi}c_g\\
\dot{c_{\omega,s}}=-i\omega c_{\omega,s}-\sqrt{\kappa/2\pi}c_s\\
$$


## 数值模拟 ##

### 计算思路 ###

##### 原方法： #####

将模型离散化：
$$
H=(\Delta-i\gamma_e/2)|e><e|+g(a|e><g|+a^{\dagger}|g><e|)+i\sqrt{\kappa/2\pi}\Sigma_{j=1}^{N}\delta\omega(a^{\dagger}b(\omega_j)-ab^{\dagger}(\omega_j))+\Sigma_{j=1}^{N}\delta\omega(\omega_j b^{\dagger}(\omega_j)b(\omega_j)) +(-i\Delta-i\gamma/2)a^{\dagger}a
$$

$$
|\Psi>=(c_g |g>|1>_a+c_e|e>|0>_a)|0>_b+c_s|s>|1>_a|0>_b+|s>|0>_a \Sigma_{j=1}^{N}\delta\omega c_{\omega_j,s}b^{\dagger}(\omega_j)|0>+|g>|0>_a \Sigma_{j=1}^{N}\delta\omega c_{\omega_j,g}b^{\dagger}(\omega_j)|0>
$$

令$\kappa'=\sqrt{\delta\omega\kappa/2\pi}, c_{j,s}=c_{\omega_j,s}\sqrt{\delta\omega}, c_{j,g}=c_{\omega_j,g}\sqrt{\delta\omega}, \omega_j=\delta\omega*(j-N/2)$

薛定谔方程计算结果：
$$
\dot{c_e}=(-i\Delta-\gamma_e/2)c_e-igc_g\\
 \dot{c_g}=(-i\Delta-\gamma/2)c_g-igc_e+\kappa'\Sigma_{j=1}^{N} c_{j,g}\\
 \dot{c_s}=(-i\Delta-\gamma/2)c_s+\kappa'\Sigma_{j=1}^{N} c_{j,s}\\
 \dot{c_{j,g}}=-i\omega_j c_{j,g}-\kappa'c_g\\
 \dot{c_{j,s}}=-i\omega_j c_{j,s}-\kappa'c_s\\
$$
利用向后差分计算时间演化即可。

##### 改进： #####

原方法的传播相位变化过大会影响精度，故设
$$
c_{j,g}'=c_{j,g}e^{i\omega_j t}\\
c_{j,s}'=c_{j,s}e^{i\omega_j t}\\
\Rightarrow\\
\dot{c_e}=(-i\Delta-\gamma_e/2)c_e-igc_g\\
 \dot{c_g}=(-i\Delta-\gamma/2)c_g-igc_e+\kappa'\Sigma_{j=1}^{N} c_{j,g}'e^{-i\omega_j t}\\
 \dot{c_s}=(-i\Delta-\gamma/2)c_s+\kappa'\Sigma_{j=1}^{N} c_{j,s}'e^{-i\omega_j t}\\
\dot{c_{j,g}'}=-\kappa'c_ge^{i\omega_j t}\\
 \dot{c_{j,s}'}=-\kappa'c_se^{i\omega_j t}\\
$$

##### 进一步改进： #####

欧拉方法在演化时间较长时会产生较大误差，故采取改良欧拉方法，采取“预估矫正”，详见链接：

[改进欧拉法(预估校正法)]: https://blog.csdn.net/honyniu/article/details/110911640

采用此种方法，可以较为精确地求解该演化方程

### 输入输出相位变化 ###

下图横坐标均为频率，纵坐标为输出脉冲的相位

取$g=5,\kappa=10,\gamma_e=0, \gamma=0, \Delta=0$，$\omega_b=1$

原子处在初态$|g>$态。

![g_g5](/g_g5.png)

原子处在初态$|s>$态。

![s_g5](/s_g5.png)

可见二者在$\omega=0$处确实有$\pi$翻转

### 损耗分析（未完成） ###

##### 光子损耗计算方式： #####

出射时：$P_{loss}=1-\sum_j(|c_{j,g}|^2+|c_{j,s}|^2)$

##### 保真度计算方式： #####

理想状态下，$c_{j,g,out}=c_{j,g,in}, c_{j,s,out}=e^{i\pi}c_{j,s,in}$（除去传播相位）

但实际上，由于脉冲存在线宽，实际上不同频率上的相位改变不同。

故为了表示保真度，我们采用

$WaveDot=\sum_j (c_{j,g,out}*c_{j,g,in}+ c_{j,s,out}*e^{i\pi}c_{j,s,in})$

$Fidelity=|WaveDot|^2$

这实际上是考虑了损耗之后的综合保真度。

#### 改变g ####

取$\kappa=10,\gamma_e=1, \gamma=0, \Delta=0$，$\omega_b=10$，入射脉冲在频域上的线宽为1，原子处在初态$\frac{1}{\sqrt{2}}(|g>+|s>)$态。

![F_gs_g(1-20)_kappa1_gammae1](/F_gs_g(1-20)_kappa1_gammae1.bmp)

横坐标为$g$ ，纵坐标为综合保真度Fidelity。



#### 改变$\kappa$ ####

取$g=5,\gamma_e=1, \gamma=0, \Delta=0$，$\omega_b=10$，入射脉冲在频域上的线宽为1，原子处在初态$\frac{1}{\sqrt{2}}(|g>+|s>)$态。

![F_gs_g10_kappa(1-10)_gammae1](/F_gs_g5_kappa(1-10)_gammae1.bmp)

横坐标为$\kappa$，纵坐标为综合保真度Fidelity。

取$g=10,\gamma_e=1, \gamma=0, \Delta=0$，$\omega_b=10$，入射脉冲在频域上的线宽为1，原子处在初态$\frac{1}{\sqrt{2}}(|g>+|s>)$态。

![F_gs_g10_kappa(1-10)_gammae1](/F_gs_g10_kappa(1-10)_gammae1.bmp)

横坐标为$\kappa$，纵坐标为综合保真度Fidelity。

#### 考虑协同系数C ####

将上述三组数据的$g,\kappa,\gamma_e$全部用于计算$C=\frac{2g^2}{\gamma\kappa}$，考虑协同系数对保真度和损耗的影响。

![F_C_gs](F_C_gs.bmp)

![L_C_gs](L_C_gs.bmp)

![Fp_C_gs](Fp_C_gs.bmp)

$Fidelityp=Fidelity/(1-Loss)$

可见损耗和保真度均为协同系数的函数（因为采用的数据中$g,\kappa$都会变化，但是只要协同系数相同损耗和保真度都相同）。

从图中还可得知，综合保真度会随着协同系数单调递增，而损耗会先增后减，峰值约在C=6.25时取到，Loss=0.426。



## 与理论结果比较 ##

##### 无展宽 #####

这里采取的理论结果参考了

[A quantum phase switch between a single solid-state spin and a photon]: https://www.nature.com/articles/nnano.2015.334

的补充材料。
$$
|Spin>=\frac{1}{\sqrt{2}}(|\uparrow>+|\downarrow>)\\
|photon_{input}>=|y>\\
|\Psi_{origin}>=\frac{1}{\sqrt{2}}(|\uparrow>+|\downarrow>)|y>\\
After\ the \ switch:\\ 
\Rightarrow 
|\Psi_{final}>=\frac{1}{\sqrt{2}}(f_{\uparrow}|\uparrow>+f_{\downarrow}|\downarrow>)|y>\\
The\ ideal \ result \ is:
|\Psi_{final,ideal}>=\frac{1}{\sqrt{2}}(|\uparrow>-|\downarrow>)|y>\\
Fidelity=|<\Psi_{final,ideal}|\Psi_{final}>|^2\\
Since\ f_{\uparrow}=\sqrt{N}(1-\frac{2\alpha}{1+C}),f_{\downarrow}=-\sqrt{N}(2\alpha-1)\\
Fidelity=N\alpha^2(\frac{C}{1+C})^2.\\
=(\frac{\kappa_{ex}}{\kappa})^2*\eta_{in}\eta_{out}*(\frac{C}{1+C})^2\\
P.S.\\
\alpha=\frac{\kappa_{ex}}{\kappa}*\frac{\sqrt{\eta_{in}\eta_{out}}}{\sqrt{\eta_{in}\eta_{out}}+\sqrt{(1-\eta_{in})\eta_{filter}}}\\
N=(\sqrt{\eta_{in}\eta_{out}}+\sqrt{(1-\eta_{in})\eta_{filter}})^2
$$
各参数含义见链接补充材料第5部分，式 (19~26)。

注：

这篇文章中的方法核心思想是将入射和出射的单光子脉冲中与腔模不符合的分量和符合的分量分开处理。腔模记作$\hat{a}$，与腔模不符合的分量（偏振不符，频率展宽失谐等）整体等价为另一个模式$\hat{v}$。因此，该方法并没有将频域上的展宽做详细讨论，这一点是否合理有待商榷。

在我的模拟中，假设腔理想，不会泄露光子，光子的损失仅来自激发态的自发辐射。故$\frac{\kappa_{ex}}{\kappa}=1$，$Fidelity=
=\eta_{in}\eta_{out}*(\frac{C}{1+C})^2\\$

且偏振条件均理想，故$\eta_{in}\eta_{out}$仅来自频域展宽。

经计算比较，模拟结果和理论结果并不十分符合。

##### 有展宽 #####

这里参考了文章：Optimal cavity design for minimizing errors in cavity-QED-based atom-photon entangling gates with finite temporal duration
$$
L_0=\frac{-\kappa_{ex}+\kappa_{in}+i\Delta}{\kappa_{ex}+\kappa_{in}+i\Delta}\\
L_1=\frac{-\kappa_{ex}+\kappa_{in}+i\Delta+\frac{g^2}{\gamma+i\Delta}}{\kappa_{ex}+\kappa_{in}+i\Delta+\frac{g^2}{\gamma+i\Delta}}\\
Fidelity=|<\Psi_{final,ideal}|\Psi_{final}>|^2=|\frac{1}{2}L_0-\frac{1}{2}L_1|^2\\
Let\ \kappa_{in}=0\\
L_0=\frac{-\kappa+i\Delta}{\kappa+i\Delta}\\
L_1=\frac{-\kappa+i\Delta+\frac{g^2}{\gamma+i\Delta}}{\kappa+i\Delta+\frac{g^2}{\gamma+i\Delta}}\\
C_{0,out}(\omega)=\frac{-\kappa+i\omega}{\kappa+i\omega}C_{0,in}(\omega)\\
C_{1,out}(\omega)=\frac{-\kappa+i\omega+\frac{g^2}{\gamma+i\omega}}{\kappa+i\omega+\frac{g^2}{\gamma+i\omega}}C_{1,in}(\omega)\\
$$






## 附录 ##

### 输出文件名称说明： ###

每一个image文件夹表示一个数据点，image中的Condition.txt保存的是这个数据点的各个参数取法，LossAnalysis.txt存的是损失和保真度的值。

abs(f_in(omega)).png 为输入频谱图；abs(fg_out(omega))与推导中的$c_{j,g}$的最终值对应；abs(fs_out(omega))与推导中的$c_{j,s}$的最终值对应；arg(fg_out(omega)).png为最终$c_{j,g}$ 的相位图，arg(fs_out(omega)).png 类似。

C_e(t).png为演化过程中的$C_e$变化曲线，C_s(t).png , C_g(t).png 类似。

n_external(t).png 为演化过程中波导上光子数变化图，n_internal()t).png为腔中的，n_total(t).png为系统中总光子数。

fomegaseries文件夹下是演化过程中每隔一段时间输出一次$c_{j,g},c_{j,s}$相关的信息，有幅值频谱图，也有相位频谱图。

![image-20221031182201911](C:\Users\李达宇\AppData\Roaming\Typora\typora-user-images\image-20221031182201911.png)
