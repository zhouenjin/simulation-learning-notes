# XPBD 学习笔记

返回：[[00_Index]]

前置：[[01_PBD]]

## XPBD 想解决什么

PBD 里有一个典型问题：

- 同样的约束
- 只改迭代次数
- 看起来软硬就不一样

XPBD 想做的是：

让“软硬”更多由约束参数本身决定，而不是主要由求解器迭代次数决定。

## XPBD 比 PBD 多了什么

### 1. `compliance`

表示约束柔度。

- `compliance = 0`：很硬
- `compliance` 越大：越软

### 2. `lambda_total`

表示某个约束在当前时间步里累计施加了多少修正。

它让 XPBD 不再只是“看当前误差猛修”，而是带一点“历史累计”的感觉。

## 为什么 `xpbd_soft` 里 error 会变大

这里的 `total length error` 更像“总伸长量”，而不只是“数值没算好”。

如果：

- 重力一直往下拉
- `compliance` 又比较大

那么 XPBD 会允许绳子逐渐伸长，最后靠近一个新的平衡状态。

所以前几帧里 `error` 逐步增大，通常是在反映：

**软约束在外力下逐渐发生形变。**

## XPBD 的核心更新

```math
\Delta \lambda =
\frac{-C(x) - \tilde{\alpha}\lambda}
{\sum_i w_i ||\nabla C_i||^2 + \tilde{\alpha}}
```

其中：

```math
\tilde{\alpha} = \frac{\alpha}{dt^2}
```

再用：

```math
\Delta x_i = w_i \nabla C_i \Delta \lambda
```

更新位置。

## 你要抓住的直觉

- PBD：当前错多少，就直接修多少
- XPBD：当前错多少 + 之前已经修了多少 + 这个约束本来有多软

## 对应实验

### 查看对比图

![[xpbd_compare.png]]

### 常用命令

```powershell
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode list
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode compare
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode plot --preset xpbd_rigid
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode plot --preset xpbd_soft
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode text --preset xpbd_soft --frames 5
```

## 重点对比

最值得对着看的是：

- `pbd_soft` vs `pbd_stiff`
- `xpbd_soft` vs `xpbd_soft_more_iters`

## 一句话过渡到 FEM

PBD / XPBD 都是在“位置投影”视角下做约束。  
而 [[03_FEM_Intro]] 会带你切换到“从能量出发”的视角。
