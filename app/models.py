"""请求/响应数据模型（Pydantic v2）。"""

from __future__ import annotations

import math
import re
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator

from .screen import (
    MULTIPLEX_MAX_PAIRS,
    MULTIPLEX_MAX_POOLS,
    MULTIPLEX_MAX_PRIMERS,
    MULTIPLEX_MIN_POOLS,
    MULTIPLEX_MIN_PRIMERS,
    VALID_BASES,
    clean_sequence,
)

ID_PATTERN = re.compile(r"^[A-Za-z0-9_.\-]{1,64}$")

# 敏感性分析：盐浓度 / Tm 阈值列表各自最多 8 项，组合不超过 64 格。
MAX_SENSITIVITY_LIST = 8
MAX_SENSITIVITY_CELLS = 64


class PrimerIn(BaseModel):
    id: str = Field(..., description="引物名称，1-64 个字母/数字/._- 字符")
    sequence: str = Field(..., description="5'->3' 引物序列，仅允许 A/C/G/T")

    @field_validator("id")
    @classmethod
    def _check_id(cls, v: str) -> str:
        v = v.strip()
        if not ID_PATTERN.match(v):
            raise ValueError("引物 id 只能包含字母、数字、点、下划线、连字符（1-64 字符）")
        return v

    @field_validator("sequence")
    @classmethod
    def _check_sequence(cls, v: str) -> str:
        seq = clean_sequence(v)
        if not 4 <= len(seq) <= 120:
            raise ValueError("引物长度需在 4-120 nt 之间")
        bad = set(seq) - VALID_BASES
        if bad:
            raise ValueError(f"引物序列含非法碱基：{''.join(sorted(bad))}（仅允许 A/C/G/T）")
        return seq


class ScreenRequest(BaseModel):
    primers: list[PrimerIn] = Field(..., min_length=1, max_length=200)
    salt_conc: float = Field(0.05, ge=0.001, le=1.5, description="一价盐浓度（M Na+ 当量）")
    temperature: float = Field(25.0, ge=0.0, le=95.0, description="退火温度 / Tm 筛选阈值（°C）")
    min_bp: int = Field(4, ge=3, le=30, description="最小连续互补 bp 数（粗筛门槛）")
    name: str | None = Field(None, max_length=120, description="本次筛选任务的备注名")


class JobSummary(BaseModel):
    job_id: str
    name: str | None
    created_at: str
    salt_conc: float
    temperature: float
    min_bp: int
    primer_count: int
    events_total: int
    high_risk_hairpins: int
    high_risk_dimers: int
    discard: list[str]


class ScreenResponse(BaseModel):
    job_id: str
    result: dict[str, Any]


class ErrorResponse(BaseModel):
    detail: str


def _validate_param_list(values: list[float], *, lo: float, hi: float, label: str) -> list[float]:
    """校验敏感性分析的参数列表：长度、有限数值、取值范围。"""

    if not values:
        raise ValueError(f"{label}列表不能为空，至少提供 1 项，最多 {MAX_SENSITIVITY_LIST} 项")
    if len(values) > MAX_SENSITIVITY_LIST:
        raise ValueError(f"{label}列表最多 {MAX_SENSITIVITY_LIST} 项，当前为 {len(values)} 项")
    cleaned: list[float] = []
    for v in values:
        try:
            f = float(v)
        except (TypeError, ValueError):
            raise ValueError(f"{label}必须为数值，收到 {v!r}")
        if not math.isfinite(f):
            raise ValueError(f"{label}必须为有限数值，收到 {v!r}")
        if not lo <= f <= hi:
            raise ValueError(f"{label}需在 {lo}-{hi} 范围内，收到 {f}")
        cleaned.append(f)
    if len(set(cleaned)) != len(cleaned):
        raise ValueError(f"{label}列表中存在重复值，请去重后再提交")
    return cleaned


class SensitivityRequest(BaseModel):
    """条件敏感性分析请求：对历史任务复用同一批引物逐格重算。"""

    salt_concs: list[float] = Field(
        ...,
        description="盐浓度列表（M Na+ 当量），1-8 项，每项 0.001-1.5，不可重复",
    )
    temperatures: list[float] = Field(
        ...,
        description="Tm 阈值列表（°C），1-8 项，每项 0-95，不可重复",
    )

    @field_validator("salt_concs")
    @classmethod
    def _check_salts(cls, v: list[float]) -> list[float]:
        return _validate_param_list(v, lo=0.001, hi=1.5, label="盐浓度")

    @field_validator("temperatures")
    @classmethod
    def _check_temps(cls, v: list[float]) -> list[float]:
        return _validate_param_list(v, lo=0.0, hi=95.0, label="Tm 阈值")

    def validate_grid(self) -> None:
        """组合格数上限校验（放在路由中显式调用，以返回 400）。"""

        cells = len(self.salt_concs) * len(self.temperatures)
        if cells > MAX_SENSITIVITY_CELLS:
            raise ValueError(
                f"盐浓度 × Tm 阈值组合为 {cells} 格，超过上限 {MAX_SENSITIVITY_CELLS} 格"
                f"（两个列表各不超过 {MAX_SENSITIVITY_LIST} 项）"
            )


class MultiplexRequest(BaseModel):
    """多重 PCR 分池请求：对历史任务中选中的引物按任务条件全量重算后分池。"""

    primer_ids: list[str] = Field(
        ...,
        min_length=MULTIPLEX_MIN_PRIMERS,
        max_length=MULTIPLEX_MAX_PRIMERS,
        description=f"选中的引物 id，{MULTIPLEX_MIN_PRIMERS}-{MULTIPLEX_MAX_PRIMERS} 条，需属于该任务",
    )
    force_pairs: list[list[str]] = Field(
        default_factory=list,
        max_length=MULTIPLEX_MAX_PAIRS,
        description=(
            "必须同池的引物对（每对两个 id），各对互不重叠、每条引物至多属于一对，"
            f"最多 {MULTIPLEX_MAX_PAIRS} 对"
        ),
    )
    max_pools: int = Field(
        MULTIPLEX_MIN_POOLS,
        ge=MULTIPLEX_MIN_POOLS,
        le=MULTIPLEX_MAX_POOLS,
        description=f"池数上限（{MULTIPLEX_MIN_POOLS}-{MULTIPLEX_MAX_POOLS}），后端搜索不超过该上限的最少池数",
    )

    @field_validator("primer_ids")
    @classmethod
    def _check_primer_ids(cls, v: list[str]) -> list[str]:
        cleaned = []
        for pid in v:
            if not isinstance(pid, str):
                raise ValueError("引物 id 必须为字符串")
            pid = pid.strip()
            if not ID_PATTERN.match(pid):
                raise ValueError(
                    f"引物 id {pid!r} 不合法（仅允许字母、数字、点、下划线、连字符，1-64 字符）"
                )
            cleaned.append(pid)
        if len(set(cleaned)) != len(cleaned):
            raise ValueError("选中引物存在重复，请去重后再提交")
        return cleaned

    @field_validator("force_pairs")
    @classmethod
    def _check_force_pairs(cls, v: list[list[str]]) -> list[list[str]]:
        cleaned: list[list[str]] = []
        for pair in v:
            if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                raise ValueError("每个必须同池引物对需为恰好两个引物 id 的数组，如 [\"F1\", \"R1\"]")
            ids: list[str] = []
            for pid in pair:
                if not isinstance(pid, str):
                    raise ValueError("引物 id 必须为字符串")
                pid = pid.strip()
                if not ID_PATTERN.match(pid):
                    raise ValueError(
                        f"引物 id {pid!r} 不合法（仅允许字母、数字、点、下划线、连字符，1-64 字符）"
                    )
                ids.append(pid)
            if ids[0] == ids[1]:
                raise ValueError(f"必须同池引物对不能由同一条引物组成：{ids[0]}")
            cleaned.append(ids)
        if len(cleaned) > MULTIPLEX_MAX_PAIRS:
            raise ValueError(f"必须同池引物对最多 {MULTIPLEX_MAX_PAIRS} 对")
        return cleaned
