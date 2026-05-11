# Abaqus 重构版代码解析

本文档对 [Project_GSRS_abaqus.py](Project_GSRS_abaqus.py) 进行逐层讲解，结合具体代码片段说明**为什么这样组织**、**每一步在做什么**、**与原始 [Project_GSRS.py](Project_GSRS.py) 相比改进在哪**。

## 整体结构

重构后脚本分为四个职责清晰的层次：

```
┌─────────────────────────────────────────┐
│  vec helpers (顶部函数)                 │  向量运算原语
├─────────────────────────────────────────┤
│  ShellGeometry (纯几何类)               │  节点坐标 + 杆件拓扑，不依赖 Abaqus
├─────────────────────────────────────────┤
│  search_best_adjustment (搜索函数)      │  复用 ShellGeometry 搜最优系数
├─────────────────────────────────────────┤
│  GeodesicShellPart (Abaqus 建模类)      │  把几何结果落地为 Part
└─────────────────────────────────────────┘
```

**核心思想**：把"算几何"和"建 Part"彻底分开。原始脚本在 `search_best_adc` 和 `Project_dcx.part` 两处重复了 100 多行几何计算——任何几何改动都要同步两处，极易遗漏。重构后几何只算一遍，两个调用方共用。

## 第一层：向量工具与数学辅助

[Project_GSRS_abaqus.py:24-57](Project_GSRS_abaqus.py#L24-L57)

```python
def vec_add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def vec_sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def vec_scale(a, k): return (a[0] * k, a[1] * k, a[2] * k)
```

纯函数，三维元组进三维元组出。对应原脚本的 `plus` / `cut` / `double`，语义一致，名字改得更通用。

### `project_on_sphere` — 核心简化

```python
def project_on_sphere(p, radius):
    length = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
    if length == 0:
        return (0.0, 0.0, 0.0)
    k = radius / length
    return (p[0] * k, p[1] * k, p[2] * k)
```

这是对原 [xyz_coordinates](Project_GSRS.py#L26-L52) 的替代。原函数用 27 行、`atan` 分象限、`math.cos` / `math.sin` 多次调用实现同一件事：**把点缩放到半径 r 的球面上**。

**等价性推导**：设 `|p|=L`、`l_xy=√(x²+y²)`。原式输出 `(r·cos(z)·|cos(xy)|·sign(x), r·cos(z)·|sin(xy)|·sign(y), r·sin(z))`。代入 `cos(z)=l_xy/L`、`sin(z)=z/L`、`|cos(xy)|=|x|/l_xy`、`|sin(xy)|=|y|/l_xy` 化简即得 `(r·x/L, r·y/L, r·z/L)`，完全一致。

**顺带修掉的 bug**：原函数在 `x=y=0` 时 `xy_angle` 未定义，又在 [Project_GSRS.py:38](Project_GSRS.py#L38) 用了 `pi` 但只 `import math`，纯属靠"那条分支从不被触达"在苟活。新版直接没有这些分支。

### `_interpolate` — 消除重复的等分循环

```python
def _interpolate(start, end, steps):
    step = vec_scale(vec_sub(end, start), 1.0 / steps)
    return [vec_add(start, vec_scale(step, j)) for j in range(steps)]
```

原脚本里形如

```python
a = double(cut(c, (0,0,1)), 1.0/f1)
coordinates = [plus((0,0,1), double(a, m)) for m in range(f1+1)]
```

这种模式在上下部重复出现 **7 次**，每次写法略有不同（端点含不含、步长取 `1/f1` 还是 `1/(f1-i)`）。统一封装后消歧义。

注意函数返回的是**左闭右开**区间：`[start, start+step, ..., end-step]`，不含 `end`。这对下部构建的"相邻 seed 列之间填充"正好合适——`end` 是下一列的起点，会在下一轮自然补上。

### `_variance`

```python
def _variance(values):
    n = len(values)
    if n == 0: return 0.0
    mean = sum(values) / n
    return sum([(v - mean) ** 2 for v in values]) / n
```

替换原脚本的 `np.var`。一来去掉 numpy 依赖（Abaqus 内嵌 Python 加载慢），二来**必须用列表推导而非生成器表达式**——Abaqus 2020+ 的内嵌 `sum` 对生成器类型做了校验，传生成器会抛 `TypeError: arg1; found 'generator', expecting a recognized type`。这是个踩过的坑。

## 第二层：`ShellGeometry` — 纯几何

[Project_GSRS_abaqus.py:62-235](Project_GSRS_abaqus.py#L62-L235)

### 构造函数：几何参数全部前置

```python
def __init__(self, f1, f2, span, rise, adjustment):
    self.f1 = f1
    self.f2 = f2
    self.radius = rise / 2.0 + (span * span / 8.0) / rise
    self.z_angle = math.acos((self.radius - rise) / self.radius)
    base_split = f1 / float(f1 + f2)
    self.z_angle_upper = (base_split + adjustment) * self.z_angle
    self.upper_layers = self._build_upper()
    self.lower_layers = self._build_lower()
```

`radius` 是由跨度 s 和矢高 h 反算的球面半径：

$$R = \frac{h}{2} + \frac{s^2}{8h}$$

`z_angle` 是球冠的半张角（从球心看，顶点到边缘的夹角）。`base_split` 把这个张角按 `f1:f2` 分给上下两部分，`adjustment` 是外部搜出来的微调量，用于均匀杆长。

**关键设计**：一个 `ShellGeometry` 实例 = 一组确定的几何。`adjustment` 通过构造参数注入，**不**作为状态被反复重设。这让"搜索最优系数"简化为"构造一系列实例取最优"。

### `_build_upper` — 上部球冠

[Project_GSRS_abaqus.py:75-105](Project_GSRS_abaqus.py#L75-L105)

**第 1 步：定义顶部辐射五束方向**

```python
apex = (0.0, 0.0, 1.0)
base_ring = [
    (math.sin(self.z_angle_upper) * math.cos(TWO_PI / PENTAGON * k),
     math.sin(self.z_angle_upper) * math.sin(TWO_PI / PENTAGON * k),
     math.cos(self.z_angle_upper))
    for k in range(PENTAGON)
]
```

从单位球顶点 `(0,0,1)` 出发，往 `z_angle_upper` 的纬度上均分 5 个方向——这就是正二十面体上顶点周围的五重对称。

**第 2 步：沿每条边从顶点线性等分到底环**

```python
edges = [
    [vec_add(apex, vec_scale(vec_sub(base, apex), m / float(f1)))
     for m in range(f1 + 1)]
    for base in base_ring
]
```

`edges[k][m]` 是第 k 条辐射边上的第 m 个分点。`edges[k][0] == apex`、`edges[k][f1] == base_ring[k]`。

**第 3 步：在相邻两条边之间再等分，得到三角形面的内部节点**

```python
layers = [[(0.0, 0.0, radius)]]  # ring 0 是顶点
for m in range(1, f1 + 1):
    ring = []
    for k in range(PENTAGON):
        a = edges[k][m]
        b = edges[(k + 1) % PENTAGON][m]
        ring.extend(_interpolate(a, b, m))
    layers.append(ring)
```

对第 m 层（`m=1..f1`），把相邻两条辐射边上的第 m 个分点 `edges[k][m]` 和 `edges[(k+1)%5][m]` 连起来再等分成 m 份——恰好得到 `PENTAGON * m` 个节点。这正是 Class I 弦分法在五重对称扇区里的标准做法。

**注意**：原脚本用 `if k+1 < PENTAGON ... else ... 0` 手写回绕，这里用 `(k + 1) % PENTAGON` 一行搞定。

**第 4 步：投影到球面**

```python
for m in range(1, f1 + 1):
    layers[m] = [project_on_sphere(p, radius) for p in layers[m]]
```

前三步得到的是正二十面体**平面三角形**上的节点，第四步才把它们径向拉到球面上。这是 Class I 方法的标准流程。

### `_build_lower` — 下部过渡带

[Project_GSRS_abaqus.py:107-172](Project_GSRS_abaqus.py#L107-L172)

下部比上部复杂得多，因为它是一个带锯齿的环带——上边界是 5 个顶点（与上部底环共享），下边界是 10 个顶点（5 主 + 5 辅交替）。

**第 1 步：三套参考环**

```python
ring_top = [
    vec_scale(self.upper_layers[f1][i], 1.0 / radius)
    for i in range(PENTAGON * f1) if i % f1 == 0
]
```

从上部底环中**每隔 f1 个取一个**——正好是 5 个原正二十面体顶点。除以 `radius` 把它们缩回单位球。

```python
ring_bot = [... 5 个下主顶点 ...]
ring_aux = [... 5 个下辅顶点 ...]
```

`ring_bot` 是下部 10 个底边顶点中的 5 个"主"顶点，`ring_aux` 是穿插其间的 5 个"辅"顶点。两者角度相差 `aux_offset`。

**第 2 步：10 条 seed 列**

```python
seeds = []
m = n = 0
while n < PENTAGON:
    bottom = ring_bot[n]
    top = ring_top[m]
    step = vec_scale(vec_sub(bottom, top), 1.0 / f2)
    seeds.append([vec_add(top, vec_scale(step, j)) for j in range(f2 + 1)])
    if abs(m - n) > 2: break
    if m > n:
        n += 1; continue
    if m == n:
        m += 1
        ring_bot[n] = ring_aux[n]   # 下一条 seed 改用辅助顶点
    if m == PENTAGON:
        m = 0
```

这段是从原脚本 [Project_GSRS.py:370-388](Project_GSRS.py#L370-L388) 照搬过来的状态机——交替切换主顶点和辅顶点，生成 10 条从上到下的"脊柱"。我保留了原有的 m/n 交替逻辑和原有行为，没有重写它（重写会冒险引入拓扑偏差）。

**第 3 步：层内填充**

```python
layers = []
total = 2 * PENTAGON
for i in range(f2 + 1):
    ring = []
    for m in range(total):
        next_m = (m + 1) % total
        if m % 2 == 0 and i != f1:
            ring.extend(_interpolate(seeds[m][i], seeds[next_m][i], f1 - i))
        elif m % 2 != 0 and i != 0:
            ring.extend(_interpolate(seeds[m][i], seeds[next_m][i], i))
    layers.append(ring)
```

第 i 层（`i=0..f2`）的环在相邻 seed 列之间填充。**偶数段**填 `f1-i` 个点（越往下越少，因为锯齿逐层收敛），**奇数段**填 `i` 个点（越往下越多，对应辅助顶点侧的扩张）。这两种计数正好合起来是每层 `PENTAGON * f1` 个节点，与上部底环对齐。

**第 4 步：投影 + 共享节点**

```python
for i in range(f2 + 1):
    for q in range(PENTAGON * f1):
        layers[i][q] = project_on_sphere(layers[i][q], radius)
layers[0] = list(self.upper_layers[f1])
```

投影到球面后，把下部的**顶环**直接替换为上部的底环——保证两部分共享节点，连续无缝。

### 杆件发射：`_collect_upper_bars` / `_collect_lower_bars`

[Project_GSRS_abaqus.py:181-235](Project_GSRS_abaqus.py#L181-L235)

这是重构的**最大技巧**：一套代码既做"计算杆长"又做"在 Abaqus 里连线"。

```python
def _collect_upper_bars(self, out, connect_only, connect=None):
    ...

@staticmethod
def _emit(out, connect_only, connect, a, b):
    if connect_only:
        connect(a, b)
    else:
        out.append(distance(a, b))
```

- 搜索系数时调用：`out=lengths, connect_only=False` → 追加杆长
- Abaqus 建模时调用：`connect=WirePolyLine, connect_only=True` → 生成线单元

**好处**：两种调用走的是同一份拓扑逻辑，从根源上保证搜索时优化的杆长和最终生成的几何**绝对一致**。原脚本搜索和建模用两份几乎相同但略有差异的代码——这种差异正是 bug 的温床。

#### 上部杆件规则

```python
for m in range(1, f1 + 1):
    ring_len = PENTAGON * m
    for q in range(ring_len):
        self._emit(..., layers[m][q], layers[m][(q + 1) % ring_len])
```

**环向杆**：每层内节点依次首尾相接。

```python
for m in range(f1, 0, -1):
    ...
    for q in range(ring_len):
        if q % m == 0:
            prev_idx = (q // m) * (m - 1) if prev_len > 0 else 0
            self._emit(..., layers[m][q], layers[m - 1][prev_idx])
        else:
            non_corner.append(q)
    for idx, q in enumerate(non_corner):
        p = layers[m][q]
        a = layers[m - 1][idx] if prev_len > 0 else layers[m - 1][0]
        b = layers[m - 1][(idx + 1) % prev_len] if prev_len > 0 else layers[m - 1][0]
        self._emit(..., p, a); self._emit(..., p, b)
```

**径向/斜向杆**：第 m 层的节点分两类。

1. **角点**（`q % m == 0`）：第 m 层每个五重扇区的角，是正二十面体棱边上的点。它只连到 m-1 层的**同位角点**。
2. **非角点**：位于扇区内部，连到 m-1 层的**两个相邻节点**，形成三角网格。

`prev_len == 0` 是 m=1 的退化情况——上一层只有顶点一个，所有连线都指向它。

#### 下部杆件规则

```python
for m in range(1, f2 + 1):
    for q in range(ring_len):
        self._emit(..., layers[m][q], layers[m][(q + 1) % ring_len])
```

**环向杆**：跳过第 0 层（那是上部底环，已经在上部代码里连好了）。

```python
for m in range(f2):
    for q in range(ring_len):
        q_prev = q - 1 if q > 0 else ring_len - 1
        self._emit(..., layers[m][q], layers[m + 1][q])
        self._emit(..., layers[m][q], layers[m + 1][q_prev])
```

**斜向杆**：每个节点连到下一层的"正下方" 和"左下方"两个邻居，形成密铺三角形。

## 第三层：`search_best_adjustment`

[Project_GSRS_abaqus.py:240-253](Project_GSRS_abaqus.py#L240-L253)

```python
def search_best_adjustment(f1, f2, span, rise, candidates=None):
    if candidates is None:
        candidates = [i / 100.0 for i in range(-2, 13)]
    best_adj = 0.0
    best_var = float("inf")
    for adj in candidates:
        geom = ShellGeometry(f1, f2, span, rise, adj)
        var = _variance(geom.bar_lengths())
        if var < best_var:
            best_var = var
            best_adj = adj
    return best_adj
```

**对比原脚本**：[search_best_adc](Project_GSRS.py#L57-L231) 有 175 行，全部是复制 `part()` 的坐标计算。重构后 **14 行**，核心就是"尝试 15 个候选系数，取方差最小那个"。

`candidates` 参数化后也便于将来改网格或做三分搜索。

## 第四层：`GeodesicShellPart` — Abaqus 适配层

[Project_GSRS_abaqus.py:258-270](Project_GSRS_abaqus.py#L258-L270)

```python
class GeodesicShellPart(object):
    def __init__(self, part_name, f1, f2, span, rise, model_name='Model-1'):
        adjustment = search_best_adjustment(f1, f2, span, rise)
        self.geometry = ShellGeometry(f1, f2, span, rise, adjustment)
        mdb.models[model_name].Part(
            name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.part = mdb.models[model_name].parts[part_name]

    def build(self):
        geom = self.geometry
        connect = lambda a, b: self.part.WirePolyLine(points=(a, b), meshable=ON)
        geom._collect_upper_bars(None, connect_only=True, connect=connect)
        geom._collect_lower_bars(None, connect_only=True, connect=connect)
```

**惟一和 Abaqus 耦合的地方**。`build()` 用一个 lambda 把 `WirePolyLine` 注入到几何类的发射循环里，几何类本身不知道 Abaqus 的存在。

这就是为什么 SAP2000 版本能轻松派生——把 `ShellGeometry` 复制过去，换一个 `AddByCoord` 的 lambda 就行。

## 调用方式

[Project_GSRS_abaqus.py:273-279](Project_GSRS_abaqus.py#L273-L279)

```python
def main(part_name, f1, f2, span, rise):
    GeodesicShellPart(part_name, f1, f2, span, rise).build()
    print("Completed!")

if __name__ == "__main__":
    main(part_name='duanchengxian_project01', f1=5, f2=2, span=20, rise=9)
```

入口保持与原脚本一致，迁移无痛。

## 总结：重构前后对比

| 维度 | 原脚本 | 重构版 |
| --- | --- | --- |
| 总行数 | 450 | 280 |
| 几何代码重复 | 搜索 + 建模各一份 | 单一来源 |
| `xyz_coordinates` | 27 行（含潜在 bug） | 5 行 |
| `search_best_adc` | 175 行（复制 part） | 14 行 |
| 依赖 | numpy | 纯 math |
| 命名 | c1/c2/c3/c4/c00/c21 | ring_top/ring_bot/seeds/layers |
| Python 2/3 | 整型除法隐患 | `from __future__ import division` |
| 扩展性 | 耦合 Abaqus | 几何可移植（已派生 SAP2000 版） |

## 设计启示

1. **分离"算什么"和"用什么算"**——纯几何类不应该知道 Abaqus、SAP2000 或任何下游消费者的存在
2. **用回调注入差异**——`_emit` + `connect` 让一份拓扑代码服务两种消费模式
3. **构造参数而非可变状态**——`adjustment` 通过 `__init__` 注入，让每次搜索都是一个干净的新实例
4. **语言习惯胜过裸循环**——`% N` 优于 `if ... else ... [0]`，列表推导优于累加循环
5. **修复和重构同步做**—在大重构里顺手把 `xyz_coordinates` 的零向量 bug、`pi` 未定义、生成器不兼容等问题都处理掉

基于python的短程线网壳建模脚本
