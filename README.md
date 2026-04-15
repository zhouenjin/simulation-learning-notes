# Rope Demos

这个目录是 `PBD / XPBD / FEM` 学习实验区。

如果你想放进 Obsidian，建议直接把整个 `rope_demos` 文件夹作为一个知识库子目录使用。

## 推荐打开顺序

1. `[[00_Index]]`
2. `[[01_PBD]]`
3. `[[02_XPBD]]`
4. `[[03_FEM_Intro]]`
5. `[[04_FEM_Corotational]]`
6. 对照 `pbd_rope_demo.py`、`xpbd_rope_demo.py`、`fem_triangle_demo.py`

## 主要文件

- `pbd_rope_demo.py`
- `xpbd_rope_demo.py`
- `pbd_compare.png`
- `xpbd_compare.png`
- `00_Index.md`
- `01_PBD.md`
- `02_XPBD.md`
- `03_FEM_Intro.md`
- `04_FEM_Corotational.md`

## 常用命令

```powershell
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode compare
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode compare
python D:\学习\编程学习\具身智能\rope_demos\pbd_rope_demo.py --mode text --preset stiff --frames 5
python D:\学习\编程学习\具身智能\rope_demos\xpbd_rope_demo.py --mode text --preset xpbd_soft --frames 5
python D:\学习\编程学习\具身智能\具身智能\rope_demos\fem_triangle_demo.py --mode all
python D:\学习\编程学习\具身智能\具身智能\rope_demos\fem_triangle_demo.py --case rotate
```

## License

本仓库内容采用 [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) 许可协议：

- 允许转载、分享、改编
- 需要保留署名
- 不允许商业使用
