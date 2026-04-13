# 具身智能物理学习索引

> 这是一套面向“数学系大二，快速上手 PBD / XPBD / FEM”的实验型笔记。

## 学习路线

推荐顺序：

1. [[01_PBD]]
2. [[02_XPBD]]
3. [[03_FEM_Intro]]
4. [[04_FEM_Corotational]]

为什么这样学：

- `PBD` 最容易做出能跑的东西
- `XPBD` 是在 PBD 上加入更稳定的“软硬”参数
- `FEM` 更接近连续介质力学和现代仿真

## 实验文件

- [[01_PBD]] 对应脚本：`pbd_rope_demo.py`
- [[02_XPBD]] 对应脚本：`xpbd_rope_demo.py`
- [[03_FEM_Intro]]、[[04_FEM_Corotational]] 对应脚本：`fem_triangle_demo.py`

## 图片总览

### PBD 参数对比

![[pbd_compare.png]]

### XPBD 参数对比

![[xpbd_compare.png]]

## 快速复习

### PBD 一句话

先预测位置，再把位置投影回满足约束的状态。

### XPBD 一句话

在 PBD 的基础上加入 `compliance` 和 `lambda_total`，让“软硬”更像材料参数。

### FEM 一句话

从单元形变能量出发，通过能量导数得到力，再更新系统。

## 我现在应该先看什么

如果你是回来复习：

- 先看 [[01_PBD]] 里的“核心骨架”
- 再看 [[02_XPBD]] 里的“为什么 error 会增长”
- 再进入 [[03_FEM_Intro]]
- 最后看 [[04_FEM_Corotational]] 里的“为什么纯旋转不该有能量”
