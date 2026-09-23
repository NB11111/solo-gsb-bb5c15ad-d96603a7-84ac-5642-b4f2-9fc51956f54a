"""SQLite 持久化层：保存筛选任务历史。"""

from __future__ import annotations

import json
import os
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path

DB_PATH = os.environ.get("PCR_DB_PATH", str(Path(__file__).resolve().parent.parent / "data" / "jobs.db"))


def _connect() -> sqlite3.Connection:
    Path(DB_PATH).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db() -> None:
    with _connect() as conn:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS jobs (
                job_id        TEXT PRIMARY KEY,
                name          TEXT,
                created_at    TEXT NOT NULL,
                request_json  TEXT NOT NULL,
                result_json   TEXT NOT NULL
            )
            """
        )
        conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_created ON jobs(created_at DESC)")


def _new_id() -> str:
    return uuid.uuid4().hex[:12]


def save_job(name: str | None, request: dict, result: dict) -> str:
    job_id = _new_id()
    created_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with _connect() as conn:
        conn.execute(
            "INSERT INTO jobs (job_id, name, created_at, request_json, result_json) "
            "VALUES (?, ?, ?, ?, ?)",
            (job_id, name, created_at, json.dumps(request, ensure_ascii=False), json.dumps(result, ensure_ascii=False)),
        )
    return job_id


def list_jobs(limit: int = 50) -> list[dict]:
    with _connect() as conn:
        rows = conn.execute(
            "SELECT job_id, name, created_at, result_json FROM jobs "
            "ORDER BY created_at DESC LIMIT ?",
            (limit,),
        ).fetchall()
    jobs = []
    for row in rows:
        result = json.loads(row["result_json"])
        jobs.append(
            {
                "job_id": row["job_id"],
                "name": row["name"],
                "created_at": row["created_at"],
                "salt_conc": result["salt_conc"],
                "temperature": result["temperature"],
                "min_bp": result["min_bp"],
                "primer_count": result["primer_count"],
                "events_total": result["events_total"],
                "high_risk_hairpins": result["summary"]["high_risk_hairpins"],
                "high_risk_dimers": result["summary"]["high_risk_dimers"],
                "discard": result["discard"],
            }
        )
    return jobs


def get_job(job_id: str) -> dict | None:
    with _connect() as conn:
        row = conn.execute(
            "SELECT job_id, name, created_at, request_json, result_json FROM jobs WHERE job_id = ?",
            (job_id,),
        ).fetchone()
    if row is None:
        return None
    return {
        "job_id": row["job_id"],
        "name": row["name"],
        "created_at": row["created_at"],
        "request": json.loads(row["request_json"]),
        "result": json.loads(row["result_json"]),
    }


def delete_job(job_id: str) -> bool:
    with _connect() as conn:
        cur = conn.execute("DELETE FROM jobs WHERE job_id = ?", (job_id,))
        return cur.rowcount > 0
