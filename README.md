# PCR 引物互补性筛选台

在网页上粘贴引物序列 JSON、设置盐浓度与退火温度，后端计算反向互补的局部
配对，分别检测 **发卡结构（Hairpin）** 与 **引物二聚体（含自二聚体）**，
按简化 Tm 阈值筛出高风险位点，返回 **可定位的配对区间**（1-based，5′→3′）
与 **建议淘汰的引物/组合**。任务历史持久化到 SQLite，结果可一键导出 JSON。

任务结果页还提供 **条件敏感性分析**：输入盐浓度与 Tm 阈值列表（各不超过 8
项、组合不超过 64 格），后端 **复用同一批引物** 逐格重算发卡/二聚体风险，
按事件类型、引物及配对区间对齐，以热力表标出 **首次达到 critical 的条件**，
并提供每格分析 JSON 下载。

任务结果页还提供 **多重 PCR 分池（Multiplex Pooling）**：从历史任务中勾选
**2～12 条引物**、指定互不重叠的 **必须同池引物对** 与 **池数上限（2～6）**，
后端按任务保存的盐浓度 / Tm 阈值 / min_bp **全量重算所选引物间的异源二聚体**
（不读取可能已截断的筛选明细），任一 high/critical 位点即形成不可同池边；
先将同池引物对折叠为分池单元，再搜索不超过上限的 **最少池数**，给出各池成员
与冲突位点区间；单元内冲突或上限内无法分池时说明原因，绝不输出无效方案。

技术栈：**FastAPI + SQLite + 原生 HTML/CSS/JS（零构建）**，Docker / Compose
一键启动。

---

## 一、快速启动

### 方式 A：Docker Compose（推荐）

```bash
docker compose up -d --build
```

启动后打开：

```
http://localhost:8000
```

停止：

```bash
docker compose down          # 停止并移除容器（数据保留在 volume）
docker compose down -v       # 同时删除任务历史数据
```

### 方式 B：直接用 Docker

```bash
docker build -t pcr-primer-screen:latest .
docker run -d --name pcr-primer-screen \
  -p 8000:8000 \
  -v pcr-primer-screen-data:/data \
  pcr-primer-screen:latest
```

### 方式 C：本地 Python 运行（开发用）

需要 Python 3.11+：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 端口与配置

| 项目 | 默认值 | 说明 |
| --- | --- | --- |
| 服务端口 | `8000` | Compose 中通过 `ports: "8000:8000"` 映射，改宿主机端口写成 `"18000:8000"` |
| `PCR_DB_PATH` | `/data/jobs.db`（容器内） | SQLite 数据库文件路径，挂载到 named volume `pcr-primer-screen-data` 持久化 |
| 盐浓度 `[Na⁺]` | `0.05` M | 页面/API 可调，范围 0.001–1.5 M |
| 退火温度 / Tm 阈值 | `25.0` °C | 页面/API 可调，范围 0–95 °C |
| 最小互补 bp 数 | `4` | 粗筛门槛，范围 3–30 |

本地直接运行时可用环境变量覆盖数据库位置：

```bash
PCR_DB_PATH=./data/jobs.db uvicorn app.main:app --reload
```

健康检查：`GET http://localhost:8000/healthz` → `{"status":"ok"}`

---

## 二、使用方法

1. 打开页面，点击 **「载入示例」**（或自行粘贴引物 JSON）；
2. 设置 **盐浓度** 与 **退火温度（Tm 阈值）**；
3. 点击 **「开始筛选」**；
4. 查看概览卡片、建议淘汰清单、发卡/二聚体位点（含坐标与三行配对展示）；
5. 在结果页底部 **「条件敏感性分析」** 输入盐浓度与 Tm 阈值列表，点击
   **「分析」** 查看热力表（★ 标出首次达到 critical 的条件），点击格子查看
   淘汰引物与相对原始条件的风险变化，并可 **「下载分析 JSON」**；
6. 在结果页底部 **「🧪 多重 PCR 分池」** 勾选 2～12 条引物、添加必须同池
   引物对（互不重叠）并设置池数上限（2～6），点击 **「计算分池」** 查看各池
   成员与冲突位点区间；调整配对或池数上限后可随时重算，并可下载分池 JSON；
7. 点击 **「导出 JSON」** 下载完整结构化结果；历史任务可在页面底部重新载入或删除。

### 输入 JSON 格式

数组形式（推荐）：

```json
[
  {"id": "Primer_F", "sequence": "ATGCGTACGTTAGCGGATCCGGATCCGGAT"},
  {"id": "Primer_R", "sequence": "TACGATCGATTACGGGATCCGGATCCGGAT"}
]
```

或扁平对象形式（页面会自动归一化）：

```json
{ "Primer_F": "ATGCGTACGTTAGCGGATCCGGATCCGGAT",
  "Primer_R": "TACGATCGATTACGGGATCCGGATCCGGAT" }
```

也可以直接提交完整请求体（含参数与备注，见 `/api/example` 的返回结构）：

```json
{
  "name": "IL-6 qPCR 初筛",
  "salt_conc": 0.05,
  "temperature": 25.0,
  "min_bp": 4,
  "primers": [ {"id": "F", "sequence": "ACGTACGTACGTACGT"} ]
}
```

序列只允许 A/C/G/T（大小写、空白均可，会自动规范化），长度 4–120 nt，
引物 id 需唯一（字母/数字/`._-`，1–64 字符）。

### 示例数据

点击「载入示例」或调用 `GET /api/example` 可取得内置示例；仓库中也放了一份
可直接提交的 [`examples/sample_request.json`](examples/sample_request.json)：

```bash
curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' \
  -d @examples/sample_request.json | python3 -m json.tool
```

示例包含：

| 引物 | 设计意图 |
| --- | --- |
| `Hairpin_F1` | 14 bp 茎 + `TTT` 环、3′ 臂贴 3′ 端 → critical 发卡 |
| `Primer_F2` / `Primer_R2` | 共享 16 nt 含回文（GGATCC）的 3′ 尾 → critical 的 3′ 锚定异源二聚体 + 自二聚体 |
| `Normal_F3` | 常规干净引物，作为对照 |

---

## 三、结果如何解读

每个位点都给出 **可定位区间**（1-based、闭区间，方向为引物自身 5′→3′）
和三行配对展示（顶行 5′→3′、中间互补标记、底行 3′→5′）。

### 发卡结构（Hairpin）

- `arm_5prime` / `arm_3prime`：形成茎区的两臂坐标；
- `loop_length`：两臂之间未配对的环长度（模型强制 ≥ 3 nt）；
- `at_edge`：茎区是否贴在 `3-prime` / `5-prime` 端，或位于链内部。
  **3′ 端发卡最危险**——聚合酶无法从被锁住的 3′ 端起延。

### 引物二聚体（Dimer）

- `primer_a` / `primer_b`：参与配对的两条引物（相同则为自二聚体）；
- `interval_a` / `interval_b`：两条引物各自的配对区间；
- `anchored_three_prime`：**双方配对区间都终止于各自 3′ 端**，即 3′ 锚定
  二聚体——聚合酶可直接沿自由 3′ 端延伸，是最差的组合；
- 自二聚体天然对称，结果已做去重，不会重复计数。

### 风险分级（简化 Tm 阈值模型）

对每个局部互补区间：

```
Tm_base = 4·(G+C) + 2·(A+T)                       # Wallace 经验式（基准 ~50 mM）
Tm      = Tm_base + 16.6·log10( f([Na⁺]) / f(0.05) )
          f(x) = x / (1 + 0.7x)                   # 相对 50 mM 基准的 SantaLucia 盐校正
```

| 级别 | 判定 |
| --- | --- |
| `critical` | 3′ 锚定（二聚体）或 3′ 端发卡，且 Tm ≥ 阈值 |
| `high` | Tm ≥ 阈值 |
| `medium` | Tm 在阈值以下 5 °C 以内且配对 ≥ 8 bp |
| `low` | 其余随机短互补，不进入报告 |

- 任引物关联到 `high` / `critical` 事件 → 进入 **建议淘汰** 清单并附原因；
- 每类事件最多保留 200 条（按风险、Tm 排序），被截断时会在响应 `capped`
  字段和页面上注明。

> 注：这是教学/初筛用途的 **简化 Tm 模型**，未使用完整最近邻热力学参数，
> 绝对 Tm 会偏保守；请以相对风险分级为主，关键引物请用专业工具复核。

### JSON 导出

页面「导出 JSON」下载的是完整响应（与 `POST /api/screen` 返回一致）：

```jsonc
{
  "job_id": "a02d91240a1b",
  "result": {
    "salt_conc": 0.05, "temperature": 25.0, "min_bp": 4,
    "primer_count": 4,
    "hairpins": [
      { "type": "hairpin", "primer_id": "Hairpin_F1",
        "arm_5prime": [1, 14], "arm_3prime": [18, 31],
        "loop_length": 3, "bp": 14, "gc_pairs": 7,
        "tm_base": 42.0, "tm": 42.0,
        "at_edge": "3-prime", "risk": "critical",
        "recommendation": "...",
        "alignment": { "top_5to3": "...", "pairs": "||||...", "bottom_3to5": "..." } }
    ],
    "dimers": [
      { "type": "dimer", "primer_a": "Primer_F2", "primer_b": "Primer_R2",
        "interval_a": [17, 30], "interval_b": [17, 30],
        "anchored_three_prime": true, "bp": 14, "gc_pairs": 8,
        "tm_base": 44.0, "tm": 44.0, "risk": "critical",
        "recommendation": "...", "alignment": { "..." : "..." } }
    ],
    "discard": ["Hairpin_F1", "Primer_F2", "Primer_R2"],
    "summary": { "...": "每引物聚合统计与淘汰原因" },
    "events_total": 12,
    "capped": { "hairpins": 0, "dimers": 0 }
  }
}
```

---

## 三·补、条件敏感性分析

针对一个 **已保存的历史任务**，提供盐浓度与 Tm 阈值两个列表，后端复用该
任务的 **同一批引物与 min_bp**，对每个（盐浓度 × Tm 阈值）格子重新计算
发卡 / 二聚体风险。

- **为什么可以复用**：哪些位置互补、配对区间多长只由序列决定；盐浓度只
  影响 Tm，Tm 阈值只影响分级。配对几何只扫描一次，逐格重算 Tm 与风险；
- **对齐方式**：每个事件以「类型（hairpin/dimer）+ 引物 + 配对区间」为
  唯一键，全部格子共用同一组事件，可直接比较同一发卡/二聚体的风险变化；
- **每格返回**：high / critical 数量（发卡、二聚体分列）、淘汰引物清单，
  以及相对任务 **原始条件（baseline）** 的变化——新进入 / 跌出高风险的
  事件、新增 / 解除淘汰的引物；
- **首次 critical（★）**：从 **最低盐浓度** 行开始，行内自 **最高 Tm
  阈值** 向左扫描，第一个出现 critical 事件的格子（即在最温和条件下仍
  critical 的边界）。每条对齐事件也给出各自的首次 critical 坐标；
- 前端以 **热力表** 展示（行=盐浓度升序、列=Tm 阈值升序），加粗边框为
  baseline 格，点击任意格子查看详情。

### 请求

`POST /api/jobs/{job_id}/sensitivity`

```json
{
  "salt_concs": [0.01, 0.05, 0.15, 0.5],
  "temperatures": [20.0, 25.0, 30.0, 35.0, 40.0]
}
```

两个列表均为 1–8 项、不可重复；盐浓度 0.001–1.5 M，Tm 阈值 0–95 °C；
组合不超过 **64 格**（8×8）。仓库提供了可直接提交的
[`examples/sensitivity_request.json`](examples/sensitivity_request.json)。

```bash
# 1) 先创建（或取一个已有）任务
JOB=$(curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' \
  -d @examples/sample_request.json | python3 -c 'import sys,json;print(json.load(sys.stdin)["job_id"])')

# 2) 对该任务执行条件敏感性分析
curl -s -X POST "http://localhost:8000/api/jobs/${JOB}/sensitivity" \
  -H 'Content-Type: application/json' \
  -d @examples/sensitivity_request.json | python3 -m json.tool
```

### 错误响应（4xx）

| 状态码 | 触发条件 |
| --- | --- |
| `404` | `job_id` 不存在或已被删除 |
| `400` | 盐浓度 × Tm 阈值组合超过 64 格 |
| `422` | 列表为空 / 超过 8 项 / 含重复值或非数值 / 超出取值范围 / 请求体不是 JSON |

### 响应结构（节选）

```jsonc
{
  "job_id": "a02d91240a1b",
  "sensitivity": {
    "salt_concs": [0.01, 0.05, 0.15, 0.5],
    "temperatures": [20.0, 25.0, 30.0, 35.0, 40.0],
    "baseline": { "salt_conc": 0.05, "temperature": 25.0, "discard": ["..."] },
    "cell_count": 20,
    "first_critical": {
      "salt_index": 0, "temp_index": 2,
      "salt_conc": 0.01, "temperature": 30.0,
      "critical_total": 4, "high_risk_total": 7
    },
    "cells": [
      {
        "salt_index": 0, "temp_index": 0, "salt_conc": 0.01, "temperature": 20.0,
        "is_baseline": false,
        "high_risk_hairpins": 0, "high_risk_dimers": 4,
        "critical_hairpins": 0, "critical_dimers": 4,
        "high_risk_total": 8, "critical_total": 4,
        "medium_hairpins": 0, "medium_dimers": 1, "events_reported": 9,
        "discard": ["Hairpin_F1", "Primer_F2", "Primer_R2"], "discard_count": 3,
        "changes_vs_baseline": {
          "new_high": [], "new_high_count": 0, "new_high_capped": false,
          "resolved": ["dimer|Primer_F2|Primer_R2|20-30|1-11"],
          "resolved_count": 1, "resolved_capped": false,
          "discard_added": [], "discard_removed": []
        }
      }
      // ... 共 20 格，顺序为 salt_concs 外层 × temperatures 内层（均升序）
    ],
    "events": [
      {
        "key": "dimer|Primer_F2|Primer_R2|17-30|17-30",
        "type": "dimer", "primer_a": "Primer_F2", "primer_b": "Primer_R2",
        "interval_a": [17, 30], "interval_b": [17, 30],
        "anchored_three_prime": true, "bp": 14, "gc_pairs": 8, "tm_base": 44.0,
        "tms": [36.5, 44.0, 48.2, 53.4],
        "risk_grid": [["critical", "critical", "..."], ["..."]],
        "first_critical": { "salt_index": 0, "temp_index": 2,
                            "salt_conc": 0.01, "temperature": 30.0 }
      }
      // 按全网格最高风险排序；超过 400 条时仅保留最高风险的 400 条，
      // 总数与是否截断见 events_total / events_shown / events_capped
    ],
    "events_total": 68, "events_shown": 68, "events_capped": false
  }
}
```

> 说明：`tms` 每个盐浓度一个值（Tm 不随 Tm 阈值变化）；`risk_grid` 为
> `salt_concs × temperatures` 的风险等级矩阵。事件键格式：
> `hairpin|引物|5′臂区间|3′臂区间`、
> `dimer|引物A|引物B|区间A|区间B`（区间均为 1-based 闭区间）。

---

## 三·补二、多重 PCR 分池

针对一个 **已保存的历史任务**，选择其中 **2～12 条引物**、指定若干
**互不重叠的必须同池引物对**（每条引物至多属于一对，最多 6 对）与
**池数上限**（2～6），后端：

1. 按任务保存的 **盐浓度 / Tm 阈值 / min_bp**，对所选引物 **全量重算**
   所有两两 **异源二聚体**（自二聚体与发卡不影响同池约束）；
2. 任一异源二聚体位点达到 `high` / `critical`，即在两条引物之间建立
   **不可同池边**。配对几何与 Tm 全部实时重算，**不从筛选任务中可能已被
   截断（每类事件上限 200 条）的事件明细推断冲突**；
3. 用并查集把必须同池引物对 **折叠为分池单元**（未配对引物自成单元）；
4. 单元内部存在冲突边（被强制同池的引物本身 high/critical）→ 直接判定
   **无解并给出冲突位点**，不输出方案；
5. 以单元为顶点建冲突图，从 1 池起做 **精确 k-着色**（DSATUR 回溯 + 颜色
   对称剪枝），返回不超过上限的 **最少池数**；
6. 上限内不可着色 → 返回精确色数下界、**最大两两冲突组（最大团）** 与全部
   冲突明细，说明需要至少多少个池，**不输出任何无效方案**。

> 分池方案为稳定输出：池按「成员中最早的引物在任务中的顺序」编号，同一输入
> 每次结果一致。页面上调整勾选、同池对或池数上限后点击「计算分池」即可重算。

### 请求

`POST /api/jobs/{job_id}/multiplex`

```json
{
  "primer_ids": ["Hairpin_F1", "Primer_F2", "Primer_R2", "Normal_F3"],
  "force_pairs": [["Hairpin_F1", "Normal_F3"]],
  "max_pools": 4
}
```

- `primer_ids`：2～12 条，必须全部属于该任务，不可重复；
- `force_pairs`：可省略或为空；每对两个不同 id；各对互不重叠；
- `max_pools`：2～6 的整数。

仓库提供可直接提交的
[`examples/multiplex_request.json`](examples/multiplex_request.json)。

```bash
JOB=$(curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' \
  -d @examples/sample_request.json | python3 -c 'import sys,json;print(json.load(sys.stdin)["job_id"])')

curl -s -X POST "http://localhost:8000/api/jobs/${JOB}/multiplex" \
  -H 'Content-Type: application/json' \
  -d @examples/multiplex_request.json | python3 -m json.tool
```

### 错误响应

| 状态码 | 触发条件 |
| --- | --- |
| `404` | `job_id` 不存在或已被删除 |
| `422` | 选中引物数量不在 2–12 / 含重复或不属于该任务；同池对结构非法、自配对、重复或互相重叠；池数上限不在 2–6（响应体 `detail` 给出具体原因） |

> 注意：`feasible: false`（单元内冲突 / 池数上限不足）属于 **正常业务结果**，
> 返回 **200** 并在响应中给出 `reason`、人类可读的 `message` 与冲突/最大团
> 诊断，而不是错误码。

### 响应结构（节选）

可行时：

```jsonc
{
  "job_id": "80868ce2f218",
  "multiplex": {
    "salt_conc": 0.05, "temperature": 25.0, "min_bp": 4,
    "primer_ids": ["Hairpin_F1", "Normal_F3", "Primer_F2", "Primer_R2"],
    "primer_count": 4,
    "units": [
      { "unit_id": "u0", "members": ["Hairpin_F1", "Normal_F3"], "forced_pair": true, "pool": 1 },
      { "unit_id": "u1", "members": ["Primer_F2"], "forced_pair": false, "pool": 1 },
      { "unit_id": "u2", "members": ["Primer_R2"], "forced_pair": false, "pool": 2 }
    ],
    "force_pair_count": 1,
    "max_pools": 4,
    "pairs_scanned": 6,
    "conflicting_primer_pairs": 1,
    "feasible": true,
    "reason": null,
    "message": "在池数上限 4 内找到最少 2 池方案。",
    "pool_count": 2,
    "pools": [
      { "pool": 1, "unit_ids": ["u0", "u1"], "members": ["Hairpin_F1", "Normal_F3", "Primer_F2"], "size": 3 },
      { "pool": 2, "unit_ids": ["u2"], "members": ["Primer_R2"], "size": 1 }
    ],
    "conflicts": [
      {
        "primer_a": "Primer_F2", "primer_b": "Primer_R2",
        "unit_a": "u1", "unit_b": "u2",
        "pool_a": 1, "pool_b": 2, "same_pool": false,
        "site_count": 2, "worst_risk": "critical",
        "sites": [ /* 与筛选接口 dimer 事件同构：区间、bp、Tm、三行配对展示… */ ]
      }
    ]
  }
}
```

不可行时（`feasible: false`）不返回 `pools`，并附带：

- `reason: "forced_pair_conflict"`：`forced_conflicts` 列出每个单元内部冲突
  的同池对及其 high/critical 位点；
- `reason: "pool_limit_exceeded"`：`chromatic_number`（冲突图精确色数）与
  `max_clique`（最大两两冲突组的单元与成员，色数下界）。

两种情况下 `conflicts` 都包含全量重算得到的 **全部** 不可同池引物对及位点
区间（含 1-based 闭区间、bp/GC/Tm 与三行配对展示），不受 200 条事件上限
影响。

---

## 四、HTTP API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/api/screen` | 提交引物与参数，执行筛选并落库 |
| `POST` | `/api/jobs/{job_id}/sensitivity` | 对历史任务执行盐浓度 × Tm 阈值条件敏感性分析 |
| `POST` | `/api/jobs/{job_id}/multiplex` | 对历史任务中 2–12 条选中引物执行多重 PCR 分池 |
| `GET` | `/api/example` | 获取内置示例请求体 |
| `GET` | `/api/jobs` | 历史任务列表（最近 50 条） |
| `GET` | `/api/jobs/{job_id}` | 某任务的请求 + 完整结果 |
| `DELETE` | `/api/jobs/{job_id}` | 删除历史任务 |
| `GET` | `/healthz` | 健康检查 |
| `GET` | `/docs` | FastAPI 自动生成的 Swagger 文档 |

命令行示例：

```bash
curl -s http://localhost:8000/api/example | curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' -d @- | python3 -m json.tool
```

---

## 五、算法说明（简版）

两条核酸链 Watson-Crick 配对时方向相反。固定链 A 按 5′→3′、链 B 在复合
物中按 3′→5′ 摆放，则 A 上位置 `i` 与 B 上位置 `k` 的互补碱基满足
`i + k = const`，即落在同一 **反对角线** 上：

1. 枚举所有反对角线（允许两链带悬挂端相对滑动）；
2. 沿对角线逐格判断互补，切出所有极大连续互补段（无 gap，结果可解释）；
3. 发卡检测：链与自身配对，丢弃两臂重叠或环长 < 3 nt 的平凡自配对；
4. 二聚体检测：自身（自二聚体，不要求环长）与两两异源配对，对称结果去重；
5. 对每个互补段计算 Wallace + 盐校正 Tm，按上表分级并聚合成淘汰建议；
6. 敏感性分析复用第 1–4 步的几何结果（配对几何只由序列与 `min_bp` 决定），
   只对盐浓度/Tm 阈值逐格重算 Tm 与分级，因此事件按类型+引物+区间天然对齐；
7. 多重 PCR 分池按任务保存的条件 **重新执行** 第 1–5 步扫描所选引物的全部
   两两异源二聚体（不读已截断明细），以 high/critical 位点为不可同池边；
   必须同池对先经并查集折叠为单元，再对单元冲突图做精确 DSATUR k-着色，
   从 1 池向上取不超过上限的最小可行池数，无解时以最大团给出色数下界。

纯标准库实现，无 numpy/生物信息学重依赖；30 条约 20–28 nt 的引物全量
（含所有两两组合）筛选约 0.4 秒。

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 路由 + 静态托管
│   ├── screen.py        # 核心：反对角线扫描、发卡/二聚体、Tm、分级
│   ├── models.py        # Pydantic 请求模型与输入校验
│   ├── db.py            # SQLite 任务历史
│   └── static/          # 原生前端 index.html / styles.css / app.js
├── data/                # 本地运行时的 SQLite 目录（已 gitignore）
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```
