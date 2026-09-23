"""PCR 引物互补性筛选核心算法。

只用 Python 标准库实现，便于在最小镜像中运行。

算法概览
========
两条核酸链之间的 Watson-Crick 配对发生在 **反平行** 方向上。固定第一条链
按 5'->3' 排列，第二条链（在实际复合物中）按 3'->5' 排列，则第一条链上
位置 i 与第二条链上位置 k 的碱基互为反向互补时，二者满足 ``i + k == const``，
即落在同一条 *反对角线 (anti-diagonal)* 上。

因此，局部互补区间的检测可以转化为：

1. 枚举两条长度分别为 n / m 的链之间的所有反对角线；
2. 沿每条对角线逐格判断碱基是否互补，得到一条 0/1 信号序列；
3. 从信号序列中切出所有「最长或较长」的连续（允许少量错配/中断）配对段；
4. 对发卡结构，同一条链与自身配对，过滤掉两臂位置重叠或环 (loop) 过短的段；
5. 用简化 Tm 公式（Wallace 经验式 + SantaLucia 盐校正）估算配对区间的
   Tm，与用户给定的退火温度阈值比较，给出风险等级。

位置约定
--------
所有对外返回的下标均为 **1-based、闭区间**，方向均为引物自身的 5'->3'。
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable

# ---------------------------------------------------------------------------
# 常量与限制
# ---------------------------------------------------------------------------

VALID_BASES = set("ACGT")

# 允许一条配对段中出现的连续非互补中断数。0 表示只接受完全连续的互补段。
# 简化模型下保持无 gap、无错配，结果最容易解释；同时我们在切段时会自动丢弃
# 被中断污染的段，不会把两次独立的互补事件错误合并。
MAX_INTERNAL_BREAKS = 0

# 默认最小互补 bp 数（比 Tm 阈值更早的一道粗筛）。
DEFAULT_MIN_BP = 4

# 发卡环（两臂之间未配对碱基）的最小长度。
MIN_HAIRPIN_LOOP = 3

# 每类事件最多保留的数量，防止序列很多时响应无限膨胀。
EVENT_CAP = 200

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------


@dataclass
class Segment:
    """一条反对角线上的局部互补段（算法内部使用，0-based）。"""

    length: int
    matches: int
    # 链 A 上的区间 [a_start, a_end)
    a_start: int
    a_end: int
    # 链 B 上的区间（按 B 自身 5'->3'）[b_start, b_end)
    # 注意：在反平行排列中 b_start 对应 a_end-1。
    b_start: int
    b_end: int


@dataclass
class HairpinEvent:
    primer_id: str
    arm5: list  # [start, end] 1-based 闭区间
    arm3: list
    loop_len: int
    bp: int
    matches: int
    gc: int
    tm_base: float
    tm: float
    at_edge: str  # "3-prime" | "5-prime" | "internal"
    risk: str
    recommendation: str
    top: str  # 5'->3'
    middle: str  # 配对标记
    bottom: str  # 3'->5' 排列展示


@dataclass
class DimerEvent:
    primer_a: str
    primer_b: str
    same_primer: bool
    interval_a: list
    interval_b: list
    a_is_three_prime: bool
    b_is_three_prime: bool
    anchored_three_prime: bool
    bp: int
    matches: int
    gc: int
    tm_base: float
    tm: float
    risk: str
    recommendation: str
    top: str
    middle: str
    bottom: str


@dataclass
class PrimerSummary:
    primer_id: str
    length: int
    hairpins: int = 0
    dimers: int = 0
    worst_risk: str = "low"
    discard: bool = False
    reasons: list = field(default_factory=list)


@dataclass
class ScreenResult:
    salt_conc: float
    temperature: float
    min_bp: int
    primer_count: int
    hairpins: list
    dimers: list
    discard: list
    summary: dict
    events_total: int
    capped: dict


# ---------------------------------------------------------------------------
# 基础工具
# ---------------------------------------------------------------------------


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def clean_sequence(seq: str) -> str:
    """规范化引物序列：去空白、大写。非法碱基由调用方负责报错。"""

    return "".join(seq.split()).upper()


def salt_correction(seq_len: int, salt_conc: float) -> float:
    """相对于 50 mM 基准的盐校正项（摄氏度），基于 SantaLucia 1998。

    Wallace 经验式（4·GC+2·AT）本就是在约 50 mM 一价盐条件下总结出来的，
    因此这里把 50 mM 时的校正归零：

        ΔTm = 16.6 · log10( f([Na]) / f(0.05) ),  f(x) = x/(1+0.7x)

    盐浓度高于 50 mM 时 Tm 上升，低于 50 mM 时下降。结果夹在 [-20, 25]。
    """

    if seq_len <= 0:
        return 0.0
    na = max(min(salt_conc, 1.5), 1e-4)

    def f(x: float) -> float:
        return x / (1.0 + 0.7 * x)

    delta = 16.6 * math.log10(f(na) / f(0.05))
    return max(-20.0, min(delta, 25.0))


def simplified_tm(seq: str, salt_conc: float) -> tuple[float, float]:
    """简化 Tm 估算。

    - ``tm_base``：经典 Wallace 规则 ``Tm = 4·GC + 2·AT``（°C），适用于
      约 14-40 nt 的短链，是引物设计中最常用的「简化」近似（基准盐浓度
      约 50 mM）；
    - ``tm``：在 Wallace 基础上加「相对于 50 mM 基准」的 SantaLucia 盐
      校正，作为给定盐浓度下的配对稳定性估计。

    返回 ``(tm_base, tm)``。
    """

    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    tm_base = 4.0 * gc + 2.0 * at
    tm = tm_base + salt_correction(len(seq), salt_conc)
    return tm_base, tm


# ---------------------------------------------------------------------------
# 反对角线扫描
# ---------------------------------------------------------------------------


def _extract_segments(
    signals: list[tuple[int, int, int]],
    max_breaks: int,
) -> list[tuple[int, int, int]]:
    """从一条对角线上的 (i, k, match) 序列中切出互补段。

    切分规则（max_breaks = 0 时）：
      - 遇到非互补格即断开；
      - 每个互补段取「极大连续互补区间」。

    ``max_breaks > 0`` 时允许段内最多出现 ``max_breaks`` 个连续非互补格，
    但为了 bp 数 / Tm 不被错配虚高，段两端会被修剪到互补位置。返回每个段的
    ``(a_start, a_end, matches)``，其中 a_end 为开区间；B 侧区间由调用方
    根据反对角线关系推导。
    """

    segments: list[tuple[int, int, int]] = []
    n = len(signals)
    p = 0
    while p < n:
        if signals[p][2] == 0:
            p += 1
            continue
        q = p
        matches = 0
        pending_break = 0
        seg_end = p  # 最后一个互补位置 +1（修剪后的右端）
        while q < n:
            if signals[q][2] == 1:
                matches += 1
                pending_break = 0
                seg_end = q + 1
            else:
                if pending_break >= max_breaks:
                    break
                pending_break += 1
            q += 1
        a_start = signals[p][0]
        a_end = signals[seg_end - 1][0] + 1
        if a_end - a_start >= 1:
            segments.append((a_start, a_end, matches))
        # 从被修剪掉的尾巴之后继续，避免重复产出同一个互补段
        p = seg_end
    return segments


def scan_antiparallel(
    seq_a: str,
    seq_b: str,
    *,
    min_bp: int = DEFAULT_MIN_BP,
    min_loop: int | None = None,
) -> list[Segment]:
    """扫描两条链之间所有反平行局部互补段。

    配对关系：链 A 的位置 i（A 按 5'->3'）与链 B 的位置 k（B 按自身
    5'->3'）配对时满足 ``i + k = c``（因为 B 在复合物中反向摆放）。

    参数
    ----
    min_loop:
        若不为 None，表示发卡模式：只保留两臂（B 侧即同链自身）之间至少
        隔着 ``min_loop`` 个未配对碱基的段。对于 A 侧索引递增的一格
        (i, k)，两臂分别是 i 和 k，需要 ``k - i - 1 >= min_loop``，即
        只保留 ``i`` 较小（5' 臂）、``k`` 较大（3' 臂）的格子。
    """

    n, m = len(seq_a), len(seq_b)
    found: list[Segment] = []

    # 反平行配对：A 的 i 与 B 的 k 配对时满足 i + k = c（B 在复合物中反向
    # 摆放）。允许两条链带悬挂端相对滑动，故 c ∈ [0, n+m-2]。
    for c in range(0, n + m - 1):
        # i 的合法范围：0 <= i < n 且 0 <= c-i < m
        i_lo = max(0, c - m + 1)
        i_hi = min(n - 1, c)
        if i_lo > i_hi:
            continue

        signals: list[tuple[int, int, int]] = []
        for i in range(i_lo, i_hi + 1):
            k = c - i
            if min_loop is not None:
                # 发卡：i 必须是 5' 臂，k 必须是 3' 臂，且环足够长。
                if k - i - 1 < min_loop:
                    signals.append((i, k, 0))
                    continue
            signals.append((i, k, 1 if seq_a[i] == _rc_lut(seq_b[k]) else 0))

        for a_start, a_end, matches in _extract_segments(signals, MAX_INTERNAL_BREAKS):
            if matches < min_bp:
                continue
            # 反对角线上 i+k=c：A 末端 a_end-1 对应 B 的较小索引
            b_start = c - (a_end - 1)
            b_end = c - a_start + 1
            found.append(
                Segment(
                    length=a_end - a_start,
                    matches=matches,
                    a_start=a_start,
                    a_end=a_end,
                    b_start=b_start,
                    b_end=b_end,
                )
            )
    return found


def _rc_lut(base: str) -> str:
    """单碱基互补（不反向）。"""

    return base.translate(_COMPLEMENT)


# ---------------------------------------------------------------------------
# 展示字符串
# ---------------------------------------------------------------------------


def alignment_strings(
    seq_a: str,
    seq_b: str,
    seg: Segment,
) -> tuple[str, str, str]:
    """生成三行配对展示。

    顶行：链 A 区间，5'->3'
    中行：``|`` 表示 Watson-Crick 互补，``.`` 表示错配（当前模型不会产生）
    底行：链 B 区间按复合物中的 3'->5' 方向摆放（即 B 自身 5'->3' 序列的反序）
    """

    top = seq_a[seg.a_start : seg.a_end]
    # B 自身 5'->3' 区间是 [b_start, b_end)，复合物中从右向左配对，
    # 所以展示时把它反转，使其按 3'->5' 排列在底行。
    bottom = seq_b[seg.b_start : seg.b_end][::-1]
    marks = "".join(
        "|" if top[i] == _rc_lut(bottom[i]) else "." for i in range(len(top))
    )
    return top, marks, bottom


# ---------------------------------------------------------------------------
# 风险分级
# ---------------------------------------------------------------------------

_RISK_ORDER = {"low": 0, "medium": 1, "high": 2, "critical": 3}


# 「中等风险」的 Tm 容差：Tm 落在 [阈值 - 该值, 阈值) 视为中等。
MEDIUM_TM_MARGIN = 5.0


def _hairpin_risk(seg: Segment, tm: float, temp: float, edge: str) -> str:
    if edge == "3-prime" and tm >= temp:
        return "critical"
    if tm >= temp:
        return "high"
    if tm >= temp - MEDIUM_TM_MARGIN:
        return "medium" if edge != "internal" or seg.matches >= 8 else "low"
    return "low"


def _dimer_risk(seg: Segment, tm: float, temp: float, anchored: bool) -> str:
    if anchored and tm >= temp:
        return "critical"
    if tm >= temp:
        return "high"
    if tm >= temp - MEDIUM_TM_MARGIN:
        # 接近阈值但未锚定 3' 端的短互补段，降为背景低风险。
        return "medium" if seg.matches >= 8 else "low"
    return "low"


def _hairpin_advice(risk: str, edge: str) -> str:
    if risk == "critical":
        return (
            f"3' 端可折回形成稳定发卡（Tm 已达阈值），会直接阻碍延伸，"
            "建议替换该引物 3' 端序列或整体淘汰。"
        )
    if risk == "high":
        return "发卡结构较稳定，建议缩短茎区、打破连续 GC 或重新设计该引物。"
    if risk == "medium":
        return "存在一定发卡倾向，建议结合实际退火条件复核，必要时调整序列。"
    return "发卡风险较低。"


def _dimer_advice(risk: str, anchored: bool, same: bool) -> str:
    pair_word = "自身互补（自二聚体）" if same else "异源二聚体"
    if risk == "critical":
        return (
            f"两条引物 3' 端互补形成{pair_word}，聚合酶可直接延伸该结构，"
            "强烈建议淘汰这对组合（至少重设计其中一条引物的 3' 端）。"
        )
    if risk == "high":
        return f"形成较稳定的{pair_word}，建议避免同时使用或重设计 3' 端。"
    if risk == "medium":
        return f"存在中等强度{pair_word}倾向，建议提高退火温度或重新设计。"
    return "二聚体风险较低。"


# ---------------------------------------------------------------------------
# 事件转换
# ---------------------------------------------------------------------------


def _gc_count(seq: str) -> int:
    return seq.count("G") + seq.count("C")


def _hairpin_edge(seq_len: int, seg: Segment) -> str:
    """发卡两臂中只要有一侧贴着链端点，即视为边缘发卡。"""

    at_three = seg.b_end == seq_len  # 3' 臂贴 3' 端
    at_five = seg.a_start == 0  # 5' 臂贴 5' 端
    if at_three:
        return "3-prime"
    if at_five:
        return "5-prime"
    return "internal"


def _segments_to_hairpins(
    primer_id: str,
    seq: str,
    segments: Iterable[Segment],
    *,
    salt_conc: float,
    temperature: float,
) -> list[HairpinEvent]:
    events: list[HairpinEvent] = []
    seen: set[tuple] = set()
    length = len(seq)
    for seg in segments:
        key = (seg.a_start, seg.a_end, seg.b_start, seg.b_end)
        if key in seen:
            continue
        seen.add(key)
        top, middle, bottom = alignment_strings(seq, seq, seg)
        tm_base, tm = simplified_tm(top, salt_conc)
        edge = _hairpin_edge(length, seg)
        risk = _hairpin_risk(seg, tm, temperature, edge)
        events.append(
            HairpinEvent(
                primer_id=primer_id,
                arm5=[seg.a_start + 1, seg.a_end],
                arm3=[seg.b_start + 1, seg.b_end],
                loop_len=seg.b_start - seg.a_end,
                bp=seg.length,
                matches=seg.matches,
                gc=_gc_count(top),
                tm_base=round(tm_base, 2),
                tm=round(tm, 2),
                at_edge=edge,
                risk=risk,
                recommendation=_hairpin_advice(risk, edge),
                top=top,
                middle=middle,
                bottom=bottom,
            )
        )
    return events


def _segments_to_dimers(
    id_a: str,
    seq_a: str,
    id_b: str,
    seq_b: str,
    segments: Iterable[Segment],
    *,
    salt_conc: float,
    temperature: float,
) -> list[DimerEvent]:
    events: list[DimerEvent] = []
    same = id_a == id_b
    seen: set[tuple] = set()
    len_a, len_b = len(seq_a), len(seq_b)
    for seg in segments:
        # 自配对时反对角线会把同一个互补段按 (i,k) / (k,i) 扫两遍，去重。
        key = tuple(sorted(((seg.a_start, seg.a_end), (seg.b_start, seg.b_end))))
        if same and key in seen:
            continue
        seen.add(key)
        top, middle, bottom = alignment_strings(seq_a, seq_b, seg)
        tm_base, tm = simplified_tm(top, salt_conc)
        a_3p = seg.a_end == len_a
        b_3p = seg.b_end == len_b
        anchored = a_3p and b_3p
        risk = _dimer_risk(seg, tm, temperature, anchored)
        events.append(
            DimerEvent(
                primer_a=id_a,
                primer_b=id_b,
                same_primer=same,
                interval_a=[seg.a_start + 1, seg.a_end],
                interval_b=[seg.b_start + 1, seg.b_end],
                a_is_three_prime=a_3p,
                b_is_three_prime=b_3p,
                anchored_three_prime=anchored,
                bp=seg.length,
                matches=seg.matches,
                gc=_gc_count(top),
                tm_base=round(tm_base, 2),
                tm=round(tm, 2),
                risk=risk,
                recommendation=_dimer_advice(risk, anchored, same),
                top=top,
                middle=middle,
                bottom=bottom,
            )
        )
    return events


# ---------------------------------------------------------------------------
# 顶层筛选流程
# ---------------------------------------------------------------------------


def _sort_events(events: list) -> list:
    return sorted(
        events,
        key=lambda e: (
            -_RISK_ORDER[e.risk],
            -(e.tm if hasattr(e, "tm") else 0),
            -e.bp,
        ),
    )


def screen_primers(
    primers: list[tuple[str, str]],
    *,
    salt_conc: float,
    temperature: float,
    min_bp: int = DEFAULT_MIN_BP,
) -> ScreenResult:
    """对一批引物执行发卡 + 二聚体筛选。

    ``primers`` 为 ``[(id, 5'->3' 序列), ...]``。
    """

    hairpins: list[HairpinEvent] = []
    dimers: list[DimerEvent] = []

    # 1) 发卡：每条链与自身做反平行配对，强制最小环长。
    for pid, seq in primers:
        segs = scan_antiparallel(
            seq,
            seq,
            min_bp=min_bp,
            min_loop=MIN_HAIRPIN_LOOP,
        )
        hairpins.extend(
            _segments_to_hairpins(
                pid,
                seq,
                segs,
                salt_conc=salt_conc,
                temperature=temperature,
            )
        )

    # 2) 二聚体：包含自二聚体（无 min_loop 限制）与两两异源配对。
    for idx, (pid_a, seq_a) in enumerate(primers):
        # 自二聚体：自身与自身，不做环长限制（末端连续互补即可）。
        self_segs = scan_antiparallel(seq_a, seq_a, min_bp=min_bp)
        dimers.extend(
            _segments_to_dimers(
                pid_a, seq_a, pid_a, seq_a, self_segs,
                salt_conc=salt_conc, temperature=temperature,
            )
        )
        for pid_b, seq_b in primers[idx + 1 :]:
            segs = scan_antiparallel(seq_a, seq_b, min_bp=min_bp)
            dimers.extend(
                _segments_to_dimers(
                    pid_a, seq_a, pid_b, seq_b, segs,
                    salt_conc=salt_conc, temperature=temperature,
                )
            )

    # 只报告通过阈值粗筛的位点（风险 medium 及以上，即 Tm 不低于阈值 10°C
    # 以上）。大量随机的 2-4 bp 低风险互补属于背景噪声，不进入报告；
    # 它们不影响引物是否被淘汰（淘汰只取决于 high / critical）。
    hairpins = [e for e in _sort_events(hairpins) if e.risk != "low"]
    dimers = [e for e in _sort_events(dimers) if e.risk != "low"]

    capped = {
        "hairpins": max(0, len(hairpins) - EVENT_CAP),
        "dimers": max(0, len(dimers) - EVENT_CAP),
    }
    hairpins = hairpins[:EVENT_CAP]
    dimers = dimers[:EVENT_CAP]

    # 3) 按引物聚合成淘汰建议。
    summaries = {
        pid: PrimerSummary(primer_id=pid, length=len(seq))
        for pid, seq in primers
    }

    def mark(pid: str, risk: str, reason: str) -> None:
        s = summaries[pid]
        if _RISK_ORDER[risk] > _RISK_ORDER[s.worst_risk]:
            s.worst_risk = risk
        if reason not in s.reasons:
            s.reasons.append(reason)

    for ev in hairpins:
        s = summaries[ev.primer_id]
        s.hairpins += 1
        mark(ev.primer_id, ev.risk, f"发卡结构（{ev.risk}，Tm={ev.tm}°C）")

    for ev in dimers:
        mark(ev.primer_a, ev.risk, f"引物二聚体 {ev.primer_a}↔{ev.primer_b}（{ev.risk}）")
        mark(ev.primer_b, ev.risk, f"引物二聚体 {ev.primer_a}↔{ev.primer_b}（{ev.risk}）")
        summaries[ev.primer_a].dimers += 1
        summaries[ev.primer_b].dimers += 1

    # high / critical 才建议淘汰。
    discard_list = []
    for pid, s in summaries.items():
        s.discard = s.worst_risk in ("high", "critical")
        if s.discard:
            discard_list.append(pid)

    high_hp = sum(1 for e in hairpins if e.risk in ("high", "critical"))
    high_dm = sum(1 for e in dimers if e.risk in ("high", "critical"))

    summary = {
        "total_primers": len(primers),
        "hairpin_events": len(hairpins),
        "dimer_events": len(dimers),
        "high_risk_hairpins": high_hp,
        "high_risk_dimers": high_dm,
        "primers_to_discard": sorted(discard_list),
        "primers": [
            {
                "primer_id": s.primer_id,
                "length": s.length,
                "hairpins": s.hairpins,
                "dimers": s.dimers,
                "worst_risk": s.worst_risk,
                "discard": s.discard,
                "reasons": s.reasons,
            }
            for s in sorted(summaries.values(), key=lambda x: x.primer_id)
        ],
    }

    return ScreenResult(
        salt_conc=salt_conc,
        temperature=temperature,
        min_bp=min_bp,
        primer_count=len(primers),
        hairpins=[_hairpin_to_dict(e) for e in hairpins],
        dimers=[_dimer_to_dict(e) for e in dimers],
        discard=sorted(discard_list),
        summary=summary,
        events_total=len(hairpins) + len(dimers),
        capped=capped,
    )


def _hairpin_to_dict(e: HairpinEvent) -> dict:
    return {
        "type": "hairpin",
        "primer_id": e.primer_id,
        "arm_5prime": e.arm5,
        "arm_3prime": e.arm3,
        "loop_length": e.loop_len,
        "bp": e.bp,
        "matches": e.matches,
        "gc_pairs": e.gc,
        "tm_base": e.tm_base,
        "tm": e.tm,
        "at_edge": e.at_edge,
        "risk": e.risk,
        "recommendation": e.recommendation,
        "alignment": {"top_5to3": e.top, "pairs": e.middle, "bottom_3to5": e.bottom},
    }


def _dimer_to_dict(e: DimerEvent) -> dict:
    return {
        "type": "dimer",
        "primer_a": e.primer_a,
        "primer_b": e.primer_b,
        "same_primer": e.same_primer,
        "interval_a": e.interval_a,
        "interval_b": e.interval_b,
        "a_is_three_prime": e.a_is_three_prime,
        "b_is_three_prime": e.b_is_three_prime,
        "anchored_three_prime": e.anchored_three_prime,
        "bp": e.bp,
        "matches": e.matches,
        "gc_pairs": e.gc,
        "tm_base": e.tm_base,
        "tm": e.tm,
        "risk": e.risk,
        "recommendation": e.recommendation,
        "alignment": {"top_5to3": e.top, "pairs": e.middle, "bottom_3to5": e.bottom},
    }


# ---------------------------------------------------------------------------
# 条件敏感性分析
# ---------------------------------------------------------------------------
#
# 配对几何（哪些位置互补、区间多长）只由序列与 min_bp 决定，**与盐浓度、
# Tm 阈值无关**；变化的只是每个配对段的 Tm（盐校正）与风险分级。因此敏感性
# 分析只做一次反对角线扫描，随后对每个（盐浓度 × Tm 阈值）格子复用同一批
# 几何事件逐格重算 Tm / 风险。事件以「类型 + 引物 + 配对区间」为键，在所有
# 格子之间天然对齐，可直接比较同一发卡/二聚体在不同条件下的风险变化。

# 对齐事件明细最多返回的条数（聚合计数始终基于全部事件，不受此上限影响）。
SENSITIVITY_EVENT_CAP = 400

# 每格风险变化清单中，新增/消退事件键名最多列出的条数（另给准确计数）。
SENSITIVITY_CHANGE_CAP = 50


@dataclass
class _GeoEvent:
    """一个与条件无关的几何配对事件（算法内部使用）。"""

    kind: str  # "hairpin" | "dimer"
    key: str
    involved: tuple  # hairpin: (pid,)；dimer: (id_a, id_b)
    seg: Segment
    tm_base: float
    gc: int
    edge: str | None  # 仅发卡
    anchored: bool | None  # 仅二聚体
    same: bool | None  # 仅二聚体


def _hairpin_geo_key(pid: str, seg: Segment) -> str:
    return f"hairpin|{pid}|{seg.a_start + 1}-{seg.a_end}|{seg.b_start + 1}-{seg.b_end}"


def _dimer_geo_key(id_a: str, id_b: str, seg: Segment) -> str:
    return (
        f"dimer|{id_a}|{id_b}|"
        f"{seg.a_start + 1}-{seg.a_end}|{seg.b_start + 1}-{seg.b_end}"
    )


def _collect_geometries(
    primers: list[tuple[str, str]],
    min_bp: int,
) -> list[_GeoEvent]:
    """扫描一次全部发卡 / 二聚体几何配对段（扫描规则与 screen_primers 一致）。"""

    geos: list[_GeoEvent] = []

    # 发卡：自身配对，强制最小环长，按区间去重。
    for pid, seq in primers:
        segs = scan_antiparallel(
            seq, seq, min_bp=min_bp, min_loop=MIN_HAIRPIN_LOOP
        )
        seen: set[tuple] = set()
        length = len(seq)
        for seg in segs:
            k = (seg.a_start, seg.a_end, seg.b_start, seg.b_end)
            if k in seen:
                continue
            seen.add(k)
            top = seq[seg.a_start : seg.a_end]
            geos.append(
                _GeoEvent(
                    kind="hairpin",
                    key=_hairpin_geo_key(pid, seg),
                    involved=(pid,),
                    seg=seg,
                    tm_base=float(4 * _gc_count(top) + 2 * (len(top) - _gc_count(top))),
                    gc=_gc_count(top),
                    edge=_hairpin_edge(length, seg),
                    anchored=None,
                    same=None,
                )
            )

    # 二聚体：自二聚体（无环长限制）+ 两两异源配对。
    for idx, (id_a, seq_a) in enumerate(primers):
        self_segs = scan_antiparallel(seq_a, seq_a, min_bp=min_bp)
        seen_self: set[tuple] = set()
        len_a = len(seq_a)
        for seg in self_segs:
            k = tuple(sorted(((seg.a_start, seg.a_end), (seg.b_start, seg.b_end))))
            if k in seen_self:
                continue
            seen_self.add(k)
            top = seq_a[seg.a_start : seg.a_end]
            geos.append(
                _GeoEvent(
                    kind="dimer",
                    key=_dimer_geo_key(id_a, id_a, seg),
                    involved=(id_a, id_a),
                    seg=seg,
                    tm_base=float(4 * _gc_count(top) + 2 * (len(top) - _gc_count(top))),
                    gc=_gc_count(top),
                    edge=None,
                    anchored=(seg.a_end == len_a and seg.b_end == len_a),
                    same=True,
                )
            )
        for id_b, seq_b in primers[idx + 1 :]:
            segs = scan_antiparallel(seq_a, seq_b, min_bp=min_bp)
            len_b = len(seq_b)
            for seg in segs:
                top = seq_a[seg.a_start : seg.a_end]
                geos.append(
                    _GeoEvent(
                        kind="dimer",
                        key=_dimer_geo_key(id_a, id_b, seg),
                        involved=(id_a, id_b),
                        seg=seg,
                        tm_base=float(
                            4 * _gc_count(top) + 2 * (len(top) - _gc_count(top))
                        ),
                        gc=_gc_count(top),
                        edge=None,
                        anchored=(seg.a_end == len_a and seg.b_end == len_b),
                        same=False,
                    )
                )
    return geos


def _geo_risk(geo: _GeoEvent, tm: float, temperature: float) -> str:
    if geo.kind == "hairpin":
        return _hairpin_risk(geo.seg, tm, temperature, geo.edge)
    return _dimer_risk(geo.seg, tm, temperature, geo.anchored)


def _sensitivity_event_descriptor(
    item: dict,
    salts: list[float],
) -> dict:
    geo: _GeoEvent = item["geo"]
    seg = geo.seg
    first = item["first_critical"]
    descriptor = {
        "key": geo.key,
        "type": geo.kind,
        "bp": seg.length,
        "matches": seg.matches,
        "gc_pairs": geo.gc,
        "tm_base": round(geo.tm_base, 2),
        # 每个盐浓度一格的 Tm（Tm 只随盐浓度变化，不随 Tm 阈值变化）。
        "tms": item["tms"],
        # 风险矩阵：行对应 salt_concs（升序），列对应 temperatures（升序）。
        "risk_grid": item["rows"],
        "first_critical": (
            {
                "salt_index": first[0],
                "temp_index": -first[1],
                "salt_conc": salts[first[0]],
                "temperature": None,  # 由调用方按索引补齐
            }
            if first is not None
            else None
        ),
    }
    if geo.kind == "hairpin":
        descriptor.update(
            {
                "primer_id": geo.involved[0],
                "arm_5prime": [seg.a_start + 1, seg.a_end],
                "arm_3prime": [seg.b_start + 1, seg.b_end],
                "at_edge": geo.edge,
            }
        )
    else:
        descriptor.update(
            {
                "primer_a": geo.involved[0],
                "primer_b": geo.involved[1],
                "same_primer": geo.same,
                "interval_a": [seg.a_start + 1, seg.a_end],
                "interval_b": [seg.b_start + 1, seg.b_end],
                "anchored_three_prime": geo.anchored,
            }
        )
    return descriptor


def sensitivity_scan(
    primers: list[tuple[str, str]],
    *,
    salt_concs: list[float],
    temperatures: list[float],
    min_bp: int = DEFAULT_MIN_BP,
    baseline_salt: float,
    baseline_temp: float,
) -> dict:
    """对同一批引物做（盐浓度 × Tm 阈值）条件敏感性分析。

    输出包含：每个格子的高风险 / critical 数量、淘汰引物、相对任务原始
    条件的风险变化，以及按「类型 + 引物 + 配对区间」对齐的事件风险矩阵。
    """

    salts = sorted(set(salt_concs))
    temps = sorted(set(temperatures))
    primer_ids = [pid for pid, _ in primers]

    geos = _collect_geometries(primers, min_bp)

    # 逐事件：几何只扫一次，Tm 只随盐浓度变，风险逐格判定。
    items: list[dict] = []
    for geo in geos:
        tms = [
            round(geo.tm_base + salt_correction(geo.seg.length, s), 2)
            for s in salts
        ]
        rows: list[list[str]] = []
        first = None  # (salt_index, -temp_index)：最小者即首次 critical
        any_reported = False
        for si, tm in enumerate(tms):
            row: list[str] = []
            for ti, t in enumerate(temps):
                risk = _geo_risk(geo, tm, t)
                row.append(risk)
                if risk != "low":
                    any_reported = True
                if risk == "critical":
                    cand = (si, -ti)
                    if first is None or cand < first:
                        first = cand
            rows.append(row)

        # 基线（任务原始条件）下非 low 的事件也保留，便于呈现「风险消退」。
        base_tm = geo.tm_base + salt_correction(geo.seg.length, baseline_salt)
        base_risk = _geo_risk(geo, base_tm, baseline_temp)
        if any_reported or base_risk != "low":
            items.append(
                {
                    "geo": geo,
                    "tms": tms,
                    "rows": rows,
                    "first_critical": first,
                    "base_risk": base_risk,
                }
            )

    base_levels = {item["geo"].key: item["base_risk"] for item in items}
    base_high = {k for k, r in base_levels.items() if r in ("high", "critical")}
    base_discard = {
        pid
        for item in items
        if item["base_risk"] in ("high", "critical")
        for pid in item["geo"].involved
    }

    cells: list[dict] = []
    first_critical_cell = None
    for si, s in enumerate(salts):
        for ti, t in enumerate(temps):
            med_h = med_d = 0
            high_h = high_d = 0
            crit_h = crit_d = 0
            reported = 0
            primer_worst = {pid: 0 for pid in primer_ids}
            high_now: set[str] = set()

            for item in items:
                risk = item["rows"][si][ti]
                if risk == "low":
                    continue
                reported += 1
                geo = item["geo"]
                level = _RISK_ORDER[risk]
                for pid in geo.involved:
                    if level > primer_worst[pid]:
                        primer_worst[pid] = level
                if risk in ("high", "critical"):
                    high_now.add(geo.key)
                if geo.kind == "hairpin":
                    if risk == "medium":
                        med_h += 1
                    elif risk == "high":
                        high_h += 1
                    else:
                        crit_h += 1
                else:
                    if risk == "medium":
                        med_d += 1
                    elif risk == "high":
                        high_d += 1
                    else:
                        crit_d += 1

            discard = sorted(
                pid for pid, lv in primer_worst.items() if lv >= _RISK_ORDER["high"]
            )
            discard_now = set(discard)
            new_high = sorted(high_now - base_high)
            resolved = sorted(base_high - high_now)

            cell = {
                "salt_index": si,
                "temp_index": ti,
                "salt_conc": s,
                "temperature": t,
                "is_baseline": s == baseline_salt and t == baseline_temp,
                "events_reported": reported,
                "medium_hairpins": med_h,
                "medium_dimers": med_d,
                "high_risk_hairpins": high_h,
                "high_risk_dimers": high_d,
                "critical_hairpins": crit_h,
                "critical_dimers": crit_d,
                "high_risk_total": high_h + high_d + crit_h + crit_d,
                "critical_total": crit_h + crit_d,
                "discard": discard,
                "discard_count": len(discard),
                "changes_vs_baseline": {
                    "new_high": new_high[:SENSITIVITY_CHANGE_CAP],
                    "new_high_count": len(new_high),
                    "new_high_capped": len(new_high) > SENSITIVITY_CHANGE_CAP,
                    "resolved": resolved[:SENSITIVITY_CHANGE_CAP],
                    "resolved_count": len(resolved),
                    "resolved_capped": len(resolved) > SENSITIVITY_CHANGE_CAP,
                    "discard_added": sorted(discard_now - base_discard),
                    "discard_removed": sorted(base_discard - discard_now),
                },
            }
            cells.append(cell)

            if crit_h + crit_d > 0:
                cand = (si, -ti)
                if first_critical_cell is None or cand < first_critical_cell[0]:
                    first_critical_cell = (cand, cell)

    # 对齐事件明细：按全网格最高风险、最高 Tm 排序后截断。
    items_sorted = sorted(
        items,
        key=lambda it: (
            -max(_RISK_ORDER[r] for row in it["rows"] for r in row),
            -max(it["tms"]),
            it["geo"].key,
        ),
    )
    events_out: list[dict] = []
    for item in items_sorted[:SENSITIVITY_EVENT_CAP]:
        descriptor = _sensitivity_event_descriptor(item, salts)
        if descriptor["first_critical"] is not None:
            descriptor["first_critical"]["temperature"] = temps[
                descriptor["first_critical"]["temp_index"]
            ]
        events_out.append(descriptor)

    first_critical = None
    if first_critical_cell is not None:
        c = first_critical_cell[1]
        first_critical = {
            "salt_index": c["salt_index"],
            "temp_index": c["temp_index"],
            "salt_conc": c["salt_conc"],
            "temperature": c["temperature"],
            "critical_total": c["critical_total"],
            "high_risk_total": c["high_risk_total"],
        }

    return {
        "min_bp": min_bp,
        "salt_concs": salts,
        "temperatures": temps,
        "baseline": {
            "salt_conc": baseline_salt,
            "temperature": baseline_temp,
            "discard": sorted(base_discard),
        },
        "cell_count": len(cells),
        "first_critical": first_critical,
        "cells": cells,
        "events": events_out,
        "events_total": len(items),
        "events_shown": len(events_out),
        "events_capped": len(items) > SENSITIVITY_EVENT_CAP,
    }


# ---------------------------------------------------------------------------
# 多重 PCR 分池
# ---------------------------------------------------------------------------
#
# 从一个历史任务中选取 2-12 条引物做多重 PCR 分池：
#
# 1. 用任务保存的盐浓度 / Tm 阈值 / min_bp **重新全量扫描**所选引物两两之间的
#    异源二聚体——筛选结果里的 dimers 列表只含 medium 及以上且最多 200 条，
#    已被截断或被判为 low 的位点不能用来推断冲突；
# 2. 任一 high / critical 异源二聚体位点构成一条「不可同池」边，同一条引物对
#    可以挂多个位点（页面展示各自的冲突区间）；
# 3. 先按用户指定、互不重叠的「必须同池」引物对把引物折叠为单元（每个单元
#    至多含 2 条引物），单元内部若存在冲突边则直接判定不可行，不产出方案；
# 4. 对单元冲突图求精确 chromatic number（顶点 ≤ 12，回溯足够快），返回不超过
#    池数上限的最少池数着色；上限内无解时说明原因。

POOL_SELECTION_MIN = 2
POOL_SELECTION_MAX = 12
POOL_LIMIT_MIN = 2
POOL_LIMIT_MAX = 6


class _DSU:
    """极简并查集，用于把必须同池的引物折叠成单元。"""

    def __init__(self, items: Iterable[str]):
        self.parent = {x: x for x in items}

    def find(self, x: str) -> str:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def _heterodimer_site(
    id_a: str,
    seq_a: str,
    id_b: str,
    seq_b: str,
    seg: Segment,
    *,
    salt_conc: float,
    temperature: float,
) -> dict:
    """单个异源二聚体位点的对外描述（1-based 区间 + 三行配对展示）。"""

    top, middle, bottom = alignment_strings(seq_a, seq_b, seg)
    tm_base, tm = simplified_tm(top, salt_conc)
    anchored = seg.a_end == len(seq_a) and seg.b_end == len(seq_b)
    risk = _dimer_risk(seg, tm, temperature, anchored)
    return {
        "interval_a": [seg.a_start + 1, seg.a_end],
        "interval_b": [seg.b_start + 1, seg.b_end],
        "bp": seg.length,
        "matches": seg.matches,
        "gc_pairs": _gc_count(top),
        "tm_base": round(tm_base, 2),
        "tm": round(tm, 2),
        "anchored_three_prime": anchored,
        "risk": risk,
        "alignment": {
            "top_5to3": top,
            "pairs": middle,
            "bottom_3to5": bottom,
        },
    }


def _max_clique(adj: dict[int, set[int]]) -> list[int]:
    """Bron-Kerbosch（带 pivot）求最大团；顶点 ≤ 12，直接枚举即可。"""

    n = len(adj)
    best: set[int] = set()

    def bk(r: set[int], p: set[int], x: set[int]) -> None:
        nonlocal best
        if len(r) + len(p) <= len(best):
            return
        if not p and not x:
            best = set(r)
            return
        u = max(p | x, key=lambda v: len(adj[v] & p))
        for v in list(p - adj[u]):
            bk(r | {v}, p & adj[v], x & adj[v])
            p = p - {v}
            x = x | {v}

    bk(set(), set(range(n)), set())
    return sorted(best)


def _exact_k_coloring(adj: dict[int, set[int]], k: int) -> list[int] | None:
    """判断图是否可 k 着色，可行时返回每个顶点的颜色编号（0-based）。

    按度数降序赋值；候选颜色只取「已用过的颜色 + 一个新颜色」，利用颜色
    对称性剪枝。顶点最多 12 个，穷尽搜索很快。
    """

    n = len(adj)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda v: (-len(adj[v]), v))
    colors = [-1] * n

    def assign(pos: int, used_colors: int) -> bool:
        if pos == n:
            return True
        v = order[pos]
        blocked = {colors[u] for u in adj[v] if colors[u] >= 0}
        # 颜色对称：只尝试已出现的颜色和恰好一个新颜色，不尝试跳过新颜色。
        for c in range(min(k, used_colors + 1)):
            if c in blocked:
                continue
            colors[v] = c
            if assign(pos + 1, max(used_colors, c + 1)):
                return True
        colors[v] = -1
        return False

    return colors if assign(0, 0) else None


def _sorted_sites(sites: list[dict]) -> list[dict]:
    return sorted(sites, key=lambda s: (-_RISK_ORDER[s["risk"]], -s["tm"], s["bp"]))


def _intra_unit_reason(intra: list[dict]) -> str:
    parts: list[str] = []
    for edge in intra:
        top = edge["sites"][0]
        parts.append(
            f"必须同池的引物对 {edge['primer_a']} 与 {edge['primer_b']} 之间存在 "
            f"{edge['site_count']} 个 high/critical 异源二聚体位点"
            f"（最严重为 {edge['worst_risk']}：{edge['primer_a']} "
            f"{top['interval_a'][0]}-{top['interval_a'][1]} ↔ {edge['primer_b']} "
            f"{top['interval_b'][0]}-{top['interval_b'][1]}，Tm {top['tm']}°C），"
            "二者无法放入同一池"
        )
    return "；".join(parts) + "。"


def multiplex_pooling(
    selected: list[tuple[str, str]],
    *,
    forced_pairs: list[tuple[str, str]],
    max_pools: int,
    salt_conc: float,
    temperature: float,
    min_bp: int = DEFAULT_MIN_BP,
) -> dict:
    """对选中的引物执行多重 PCR 分池。

    ``selected`` 为按用户选择顺序排列的 ``[(id, 5'->3' 序列), ...]``，
    长度 2-12；``forced_pairs`` 为互不重叠的必须同池引物对；``max_pools``
    为池数上限（2-6）。冲突位点完全按给定的盐浓度 / Tm 阈值 / min_bp
    重新扫描，不读取历史筛选结果中可能已被截断的事件明细。
    """

    ids = [pid for pid, _ in selected]

    # 1) 先把必须同池引物对折叠为单元（迭代顺序即选择顺序，单元编号稳定）。
    dsu = _DSU(ids)
    for a, b in forced_pairs:
        dsu.union(a, b)
    root_to_idx: dict[str, int] = {}
    units: list[dict] = []
    for pid in ids:
        root = dsu.find(pid)
        if root not in root_to_idx:
            root_to_idx[root] = len(units)
            units.append({"unit_id": f"U{len(units) + 1}", "members": [], "forced": False})
        units[root_to_idx[root]]["members"].append(pid)
    forced_pair_set = {frozenset(pair) for pair in forced_pairs}
    for unit in units:
        unit["forced"] = frozenset(unit["members"]) in forced_pair_set
    unit_idx_of = {
        pid: idx for idx, unit in enumerate(units) for pid in unit["members"]
    }

    # 2) 全量重算所选引物两两异源二聚体，high/critical 位点才构成冲突边。
    primer_edges: dict[tuple[int, int], list[dict]] = {}
    for i, (id_a, seq_a) in enumerate(selected):
        for j in range(i + 1, len(selected)):
            id_b, seq_b = selected[j]
            segs = scan_antiparallel(seq_a, seq_b, min_bp=min_bp)
            sites = []
            for seg in segs:
                site = _heterodimer_site(
                    id_a, seq_a, id_b, seq_b, seg,
                    salt_conc=salt_conc, temperature=temperature,
                )
                if site["risk"] in ("high", "critical"):
                    site["primer_a"] = id_a
                    site["primer_b"] = id_b
                    sites.append(site)
            if sites:
                primer_edges[(i, j)] = _sorted_sites(sites)

    # 3) 折叠到单元层：区分单元内冲突（直接不可行）与单元间冲突（着色约束）。
    intra: list[dict] = []
    inter: dict[tuple[int, int], list[dict]] = {}
    for (i, j), sites in primer_edges.items():
        id_a, id_b = ids[i], ids[j]
        ui, uj = unit_idx_of[id_a], unit_idx_of[id_b]
        if ui == uj:
            intra.append(
                {
                    "primer_a": id_a,
                    "primer_b": id_b,
                    "unit_a": units[ui]["unit_id"],
                    "unit_b": units[uj]["unit_id"],
                    "site_count": len(sites),
                    "worst_risk": sites[0]["risk"],
                    "sites": sites,
                }
            )
        else:
            key = (min(ui, uj), max(ui, uj))
            inter.setdefault(key, []).extend(sites)

    intra.sort(key=lambda e: (-_RISK_ORDER[e["worst_risk"]], e["primer_a"], e["primer_b"]))

    adj: dict[int, set[int]] = {i: set() for i in range(len(units))}
    for ui, uj in inter:
        adj[ui].add(uj)
        adj[uj].add(ui)

    # 4) 最少池数：先求最大团作下界，再从下界起逐个尝试精确 k 着色。
    feasible = True
    infeasible_reason: str | None = None
    coloring: list[int] | None = None
    if intra:
        feasible = False
        infeasible_reason = _intra_unit_reason(intra)
    else:
        clique = _max_clique(adj)
        for k in range(max(1, len(clique)), max_pools + 1):
            coloring = _exact_k_coloring(adj, k)
            if coloring is not None:
                break
        else:
            coloring = None
        if coloring is None:
            feasible = False
            if len(clique) > max_pools:
                members = [pid for i in clique for pid in units[i]["members"]]
                infeasible_reason = (
                    f"冲突图至少需要 {len(clique)} 个池：引物 "
                    f"{'、'.join(members)} 两两之间都存在 high/critical 异源二聚体，"
                    f"互为冲突的 {len(clique)} 个单元必须各占一池，超过池数上限 "
                    f"{max_pools}。"
                )
            else:
                infeasible_reason = (
                    f"在池数上限 {max_pools} 内无解：穷尽搜索表明该冲突图至少需要 "
                    f"{max_pools + 1} 个池，请放宽必须同池配对、提高池数上限或调整引物。"
                )

    pools: list[dict] = []
    if feasible and coloring is not None:
        grouped: dict[int, list[int]] = {}
        for ui, color in enumerate(coloring):
            grouped.setdefault(color, []).append(ui)
        for color in sorted(grouped):
            members = [
                pid
                for ui in grouped[color]
                for pid in units[ui]["members"]
            ]
            pools.append(
                {
                    "pool": color + 1,
                    "unit_ids": [units[ui]["unit_id"] for ui in grouped[color]],
                    "members": members,
                }
            )

    def pool_of(ui: int) -> int | None:
        return coloring[ui] + 1 if coloring is not None else None

    conflicts: list[dict] = []
    for edge in intra:
        conflicts.append(
            {
                **edge,
                "same_unit": True,
                "pool_a": None,
                "pool_b": None,
            }
        )
    for (ui, uj), sites in inter.items():
        sites = _sorted_sites(sites)
        conflicts.append(
            {
                "primer_a": ids_in_unit(units, ui),
                "primer_b": ids_in_unit(units, uj),
                "unit_a": units[ui]["unit_id"],
                "unit_b": units[uj]["unit_id"],
                "site_count": len(sites),
                "worst_risk": sites[0]["risk"],
                "sites": sites,
                "same_unit": False,
                "pool_a": pool_of(ui),
                "pool_b": pool_of(uj),
            }
        )
    conflicts.sort(
        key=lambda e: (
            0 if e["same_unit"] else 1,
            -_RISK_ORDER[e["worst_risk"]],
            e["unit_a"],
            e["unit_b"],
        )
    )

    return {
        "feasible": feasible,
        "params": {
            "salt_conc": salt_conc,
            "temperature": temperature,
            "min_bp": min_bp,
        },
        "primer_ids": ids,
        "max_pools": max_pools,
        "units": units,
        "pool_count": len(pools) if feasible else None,
        "pools": pools,
        "conflicts": conflicts,
        "infeasible_reason": infeasible_reason,
    }


def ids_in_unit(units: list[dict], idx: int) -> str:
    """冲突边展示用：单元只有 1 条引物时返回该 id，否则返回 "A/B" 简述。"""

    members = units[idx]["members"]
    return members[0] if len(members) == 1 else "/".join(members)
