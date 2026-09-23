"""FastAPI 主应用：PCR 引物互补性筛选台。"""

from __future__ import annotations

from contextlib import asynccontextmanager
from dataclasses import asdict
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from . import db
from .models import PoolingRequest, ScreenRequest, SensitivityRequest
from .screen import multiplex_pooling, screen_primers, sensitivity_scan

BASE_DIR = Path(__file__).resolve().parent
STATIC_DIR = BASE_DIR / "static"


@asynccontextmanager
async def lifespan(_app: FastAPI):
    db.init_db()
    yield


app = FastAPI(
    title="PCR 引物互补性筛选台",
    description=(
        "粘贴引物 JSON，设置盐浓度与退火温度，检测发卡结构与引物二聚体，"
        "按简化 Tm 阈值返回可定位的高风险配对区间与淘汰建议。"
    ),
    version="1.0.0",
    lifespan=lifespan,
)


# ---------------------------------------------------------------------------
# 示例数据
# ---------------------------------------------------------------------------

EXAMPLE_PRIMERS = [
    # 3' 端发卡：14 bp 茎 + TTT 环，3' 臂直接贴在 3' 端
    {
        "id": "Hairpin_F1",
        "sequence": "GCATCGGATCGAATTTTATTCGATCCGATGC",
    },
    # 两条共享 16 nt 3' 尾（含回文 BamHI 位点 GGATCC）的引物：
    # 会形成 3' 锚定异源二聚体，且各自存在自二聚体
    {
        "id": "Primer_F2",
        "sequence": "ATGCGTACGTTAGCGGATCCGGATCCGGAT",
    },
    {
        "id": "Primer_R2",
        "sequence": "TACGATCGATTACGGGATCCGGATCCGGAT",
    },
    # 一条常规、较干净的引物
    {
        "id": "Normal_F3",
        "sequence": "GCTACAGTACGAGCTTGATC",
    },
]

EXAMPLE_PAYLOAD = {
    "name": "示例任务：1 条发卡引物 + 1 对 3' 互补引物 + 1 条正常引物",
    "salt_conc": 0.05,
    "temperature": 25.0,
    "min_bp": 4,
    "primers": EXAMPLE_PRIMERS,
}


@app.get("/api/example")
def get_example() -> dict:
    """返回可直接粘贴/一键载入的示例请求数据。"""

    return EXAMPLE_PAYLOAD


# ---------------------------------------------------------------------------
# 筛选接口
# ---------------------------------------------------------------------------


@app.post("/api/screen")
def post_screen(payload: ScreenRequest) -> dict:
    ids = [p.id for p in payload.primers]
    if len(set(ids)) != len(ids):
        raise HTTPException(status_code=422, detail="引物 id 存在重复，请保证每条引物 id 唯一")

    primers = [(p.id, p.sequence) for p in payload.primers]
    result = screen_primers(
        primers,
        salt_conc=payload.salt_conc,
        temperature=payload.temperature,
        min_bp=payload.min_bp,
    )
    result_dict = asdict(result)

    request_record = payload.model_dump()
    job_id = db.save_job(payload.name, request_record, result_dict)
    return {"job_id": job_id, "result": result_dict}


@app.get("/api/jobs")
def get_jobs() -> dict:
    return {"jobs": db.list_jobs()}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str) -> dict:
    job = db.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="任务不存在或已被删除")
    return job


@app.delete("/api/jobs/{job_id}")
def remove_job(job_id: str) -> dict:
    if not db.delete_job(job_id):
        raise HTTPException(status_code=404, detail="任务不存在或已被删除")
    return {"deleted": job_id}


# ---------------------------------------------------------------------------
# 条件敏感性分析（复用历史任务的同一批引物，盐浓度 × Tm 阈值逐格重算）
# ---------------------------------------------------------------------------


@app.post("/api/jobs/{job_id}/sensitivity")
def post_sensitivity(job_id: str, payload: SensitivityRequest) -> dict:
    job = db.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="任务不存在或已被删除，无法执行敏感性分析")

    # 理论上 Pydantic 的列表长度限制已保证 ≤ 64 格，这里保留显式 400 防线。
    try:
        payload.validate_grid()
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    req = job["request"]
    primers = [(p["id"], p["sequence"]) for p in req["primers"]]
    analysis = sensitivity_scan(
        primers,
        salt_concs=payload.salt_concs,
        temperatures=payload.temperatures,
        min_bp=int(req.get("min_bp", 4)),
        baseline_salt=float(req["salt_conc"]),
        baseline_temp=float(req["temperature"]),
    )
    return {
        "job_id": job_id,
        "job_name": job["name"],
        "sensitivity": analysis,
    }


# ---------------------------------------------------------------------------
# 多重 PCR 分池（复用历史任务的引物与保存的盐浓度 / 温度 / min_bp，全量重算
# 所选引物两两异源二聚体；保留原有筛选与敏感性分析接口）
# ---------------------------------------------------------------------------


@app.post("/api/jobs/{job_id}/pooling")
def post_pooling(job_id: str, payload: PoolingRequest) -> dict:
    job = db.get_job(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="任务不存在或已被删除，无法执行多重 PCR 分池")

    req = job["request"]
    primer_map = {p["id"]: p["sequence"] for p in req["primers"]}

    missing = [pid for pid in payload.primer_ids if pid not in primer_map]
    if missing:
        raise HTTPException(
            status_code=400,
            detail="以下引物 id 不属于该任务：" + "、".join(missing),
        )

    forced_pairs = [tuple(pair) for pair in payload.forced_pairs]
    pair_ids = {pid for pair in forced_pairs for pid in pair}
    not_selected = sorted(pair_ids - set(payload.primer_ids))
    if not_selected:
        raise HTTPException(
            status_code=400,
            detail="必须同池引物对中的以下引物未在所选引物中：" + "、".join(not_selected),
        )

    # 严格按用户选择顺序传入（保持池内成员展示顺序与选择顺序一致）。
    selected = [(pid, primer_map[pid]) for pid in payload.primer_ids]
    analysis = multiplex_pooling(
        selected,
        forced_pairs=forced_pairs,
        max_pools=payload.max_pools,
        salt_conc=float(req["salt_conc"]),
        temperature=float(req["temperature"]),
        min_bp=int(req.get("min_bp", 4)),
    )
    return {
        "job_id": job_id,
        "job_name": job["name"],
        "pooling": analysis,
    }


@app.get("/healthz")
def healthz() -> dict:
    return {"status": "ok"}


# ---------------------------------------------------------------------------
# 静态前端（原生 HTML/CSS/JS，无构建步骤）
# ---------------------------------------------------------------------------


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
