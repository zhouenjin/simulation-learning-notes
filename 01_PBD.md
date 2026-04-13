# PBD 学习笔记

返回：[[00_Index]]

## PBD 是什么

`PBD = Position-Based Dynamics`

核心思想：

1. 先按速度和外力预测位置
2. 再根据约束直接修正位置
3. 最后用修正后的位置反推速度

## 一帧流程

```python
for each particle:
    v += dt * acceleration
    x_old = x
    x = x + dt * v

repeat K times:
    solve constraints

for each particle:
    v = (x - x_old) / dt
```

## 距离约束

相邻两个点希望满足：

```math
\|x_1 - x_2\| = L
```

写成约束函数：

```math
C(x_1, x_2) = \|x_1 - x_2\| - L = 0
```

## 你应该记住的现象

- 迭代次数少：绳子更软
- 迭代次数多：绳子更硬
- 固定点本质上是 `inv_mass = 0`
- 逆质量越大，越容易在投影时被改位置

## 为什么 `soft` 和 `stiff` 不一样

因为 PBD 不是一次精确解所有约束，而是反复局部修正。

所以：

- 迭代少，约束满足得不够彻底
- 迭代多，约束满足得更充分

这就是为什么 PBD 里的“软硬”会依赖迭代次数。

## 对应实验

### 查看对比图

![[pbd_compare.png]]

### 常用命令

```powershell
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode list
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode compare
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode plot --preset soft
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode plot --preset stiff
```

## 代码里最关键的地方

你之后复习时，重点看：

- `step()`
- `solve_distance_constraint()`

## 一句话过渡到 XPBD

PBD 好用，但它的“刚度”会依赖迭代次数和时间步长，所以我们需要 [[02_XPBD]]。
