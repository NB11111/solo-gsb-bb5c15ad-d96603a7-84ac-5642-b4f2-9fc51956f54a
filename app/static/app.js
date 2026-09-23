/* PCR 引物互补性筛选台 —— 原生 JS，无依赖、无构建。 */
"use strict";

const $ = (sel) => document.querySelector(sel);

const els = {
  jobName: $("#jobName"),
  salt: $("#salt"),
  temp: $("#temp"),
  minbp: $("#minbp"),
  json: $("#jsonInput"),
  msg: $("#inputMsg"),
  run: $("#btnRun"),
  example: $("#btnLoadExample"),
  format: $("#btnFormat"),
  exportBtn: $("#btnExport"),
  refreshJobs: $("#btnRefreshJobs"),
  empty: $("#emptyResult"),
  body: $("#resultBody"),
  overview: $("#overview"),
  discardBox: $("#discardBox"),
  hairpins: $("#hairpinList"),
  dimers: $("#dimerList"),
  jobList: $("#jobList"),
  sensBox: $("#sensitivityBox"),
  sensSalts: $("#sensSalts"),
  sensTemps: $("#sensTemps"),
  sensRun: $("#btnSensRun"),
  sensMsg: $("#sensMsg"),
  sensResult: $("#sensResult"),
  sensExport: $("#btnSensExport"),
  poolBox: $("#poolingBox"),
  poolPrimerList: $("#poolPrimerList"),
  poolPairs: $("#poolPairs"),
  poolMax: $("#poolMax"),
  poolRun: $("#btnPoolRun"),
  poolMsg: $("#poolMsg"),
  poolResult: $("#poolResult"),
  poolExport: $("#btnPoolExport"),
};

let lastResponse = null; // { job_id, result }
let lastSensitivity = null; // { job_id, sensitivity }
let lastPooling = null; // { job_id, pooling }
let poolPrimers = []; // 当前任务引物 [{id, sequence, risk?}]
let poolChecked = new Set(); // 勾选的引物 id

/* ---------------- 工具 ---------------- */

function setMsg(text, kind = "") {
  els.msg.textContent = text;
  els.msg.className = "msg" + (kind ? " " + kind : "");
}

const RISK_CN = { critical: "极高", high: "高", medium: "中", low: "低" };

function badge(text, cls) {
  return `<span class="badge ${cls}">${text}</span>`;
}

function esc(s) {
  return String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

/* 接受数组形式或扁平 id->seq 形式的 JSON。 */
function normalizePrimers(parsed) {
  if (Array.isArray(parsed)) {
    return parsed.map((p) => {
      if (!p || typeof p.id !== "string" || typeof p.sequence !== "string") {
        throw new Error("数组中每项需为 {\"id\": ..., \"sequence\": ...}");
      }
      return { id: p.id, sequence: p.sequence };
    });
  }
  if (parsed && typeof parsed === "object") {
    return Object.entries(parsed).map(([id, sequence]) => {
      if (typeof sequence !== "string") throw new Error(`引物 ${id} 的序列必须是字符串`);
      return { id, sequence };
    });
  }
  throw new Error("JSON 需为引物数组或 id->sequence 对象");
}

/* ---------------- 渲染 ---------------- */

function alignmentBlock(al) {
  return `<div class="align">` +
    `<span class="dir">5'→3'</span> ${esc(al.top_5to3)}\n` +
    `             ${esc(al.pairs)}\n` +
    `<span class="dir">3'←5'</span> ${esc(al.bottom_3to5)}` +
    `</div>`;
}

function renderHairpin(e) {
  const edgeTag =
    e.at_edge === "3-prime" ? badge("3′ 端锚定", "tag")
    : e.at_edge === "5-prime" ? badge("5′ 端", "tag")
    : badge("链内部", "tag");
  return `<div class="event ${e.risk}">
    <div class="event-title">
      ${badge(RISK_CN[e.risk] + "风险", e.risk)}
      <span class="who">${esc(e.primer_id)}</span>
      ${edgeTag}
    </div>
    <div class="event-meta">
      <span>5′ 臂 <b>${e.arm_5prime[0]}–${e.arm_5prime[1]}</b></span>
      <span>3′ 臂 <b>${e.arm_3prime[0]}–${e.arm_3prime[1]}</b></span>
      <span>环长 <b>${e.loop_length}</b> nt</span>
      <span>茎区 <b>${e.bp}</b> bp</span>
      <span>GC <b>${e.gc_pairs}</b></span>
      <span>Tm <b>${e.tm}°C</b> <span style="opacity:.7">(无盐 ${e.tm_base}°C)</span></span>
    </div>
    ${alignmentBlock(e.alignment)}
    <div class="advice">💡 ${esc(e.recommendation)}</div>
  </div>`;
}

function renderDimer(e) {
  const aTag = e.a_is_three_prime ? badge("A 侧 3′", "tag") : "";
  const bTag = e.b_is_three_prime ? badge("B 侧 3′", "tag") : "";
  const selfTag = e.same_primer ? badge("自二聚体", "tag") : badge("异源二聚体", "tag");
  const anchorTag = e.anchored_three_prime ? badge("3′ 锚定 · 可延伸", "critical") : "";
  return `<div class="event ${e.risk}">
    <div class="event-title">
      ${badge(RISK_CN[e.risk] + "风险", e.risk)}
      <span class="who">${esc(e.primer_a)}</span>
      <span style="color:var(--muted)">×</span>
      <span class="who">${esc(e.primer_b)}</span>
      ${selfTag}${anchorTag}${aTag}${bTag}
    </div>
    <div class="event-meta">
      <span>A 区间 <b>${e.interval_a[0]}–${e.interval_a[1]}</b></span>
      <span>B 区间 <b>${e.interval_b[0]}–${e.interval_b[1]}</b></span>
      <span>配对 <b>${e.bp}</b> bp</span>
      <span>GC <b>${e.gc_pairs}</b></span>
      <span>Tm <b>${e.tm}°C</b> <span style="opacity:.7">(无盐 ${e.tm_base}°C)</span></span>
    </div>
    ${alignmentBlock(e.alignment)}
    <div class="advice">💡 ${esc(e.recommendation)}</div>
  </div>`;
}

function renderResult(data) {
  lastResponse = data;
  const r = data.result;
  els.empty.hidden = true;
  els.body.hidden = false;
  els.exportBtn.disabled = false;

  // 条件敏感性分析需要可引用的历史任务；纯导出态（无 job_id）时隐藏。
  resetSensitivity();
  resetPooling();
  renderPoolPrimerList();

  const cards = [
    [r.primer_count, "引物数量"],
    [r.summary.hairpin_events, "发卡位点"],
    [r.summary.dimer_events, "二聚体位点"],
    [r.summary.high_risk_hairpins, "高风险发卡"],
    [r.summary.high_risk_dimers, "高风险二聚体"],
    [r.discard.length, "建议淘汰"],
  ];
  els.overview.innerHTML = cards
    .map(([num, lbl]) => `<div class="card"><div class="num">${num}</div><div class="lbl">${lbl}</div></div>`)
    .join("");

  if (r.discard.length) {
    const perPrimer = r.summary.primers.filter((p) => p.discard);
    els.discardBox.className = "discard-box";
    els.discardBox.innerHTML = `
      <h4>🚫 建议淘汰的引物（${r.discard.length}）</h4>
      <div class="ids">${r.discard.map(esc).join("，")}</div>
      <ul>${perPrimer
        .map((p) => `<li><b>${esc(p.primer_id)}</b>：${p.reasons.map(esc).join("；")}</li>`)
        .join("")}</ul>`;
  } else {
    els.discardBox.className = "discard-box safe";
    els.discardBox.innerHTML = `<h4>✅ 未检出高风险组合</h4>
      <div>所有引物在当前盐浓度（${r.salt_conc} M）与 Tm 阈值（${r.temperature}°C）下均可保留。</div>`;
  }

  els.hairpins.innerHTML = r.hairpins.length
    ? r.hairpins.map(renderHairpin).join("")
    : `<div class="empty" style="padding:10px 0">未检出超过阈值的发卡结构。</div>`;
  els.dimers.innerHTML = r.dimers.length
    ? r.dimers.map(renderDimer).join("")
    : `<div class="empty" style="padding:10px 0">未检出超过阈值的引物二聚体。</div>`;

  if (r.capped && (r.capped.hairpins || r.capped.dimers)) {
    els.dimers.insertAdjacentHTML(
      "beforebegin",
      `<p class="hint">注：另有 ${r.capped.hairpins} 个发卡 / ${r.capped.dimers} 个二聚体事件因数量上限被截断。</p>`
    );
  }
}

/* ---------------- 动作 ---------------- */

async function runScreen() {
  let parsed;
  try {
    parsed = JSON.parse(els.json.value);
  } catch (e) {
    setMsg("JSON 解析失败：" + e.message, "err");
    return;
  }
  let primers;
  try {
    primers = normalizePrimers(parsed);
  } catch (e) {
    setMsg(e.message, "err");
    return;
  }
  const payload = {
    name: els.jobName.value.trim() || null,
    salt_conc: parseFloat(els.salt.value),
    temperature: parseFloat(els.temp.value),
    min_bp: parseInt(els.minbp.value, 10),
    primers,
  };
  els.run.disabled = true;
  setMsg("筛选中…");
  try {
    const resp = await fetch("/api/screen", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const body = await resp.json();
    if (!resp.ok) {
      setMsg("请求失败（" + resp.status + "）：" + (body.detail || ""), "err");
      return;
    }
    setMsg(`完成 ✓ 任务 ID ${body.job_id}（已存入 SQLite）`, "ok");
    renderResult(body);
    loadJobs();
  } catch (e) {
    setMsg("网络错误：" + e.message, "err");
  } finally {
    els.run.disabled = false;
  }
}

function exportJson() {
  if (!lastResponse) return;
  const blob = new Blob([JSON.stringify(lastResponse, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `primer-screen-${lastResponse.job_id}.json`;
  a.click();
  URL.revokeObjectURL(a.href);
}

async function loadExample() {
  const ex = await (await fetch("/api/example")).json();
  els.jobName.value = ex.name || "";
  els.salt.value = ex.salt_conc;
  els.temp.value = ex.temperature;
  els.minbp.value = ex.min_bp;
  els.json.value = JSON.stringify(ex.primers, null, 2);
  setMsg("已载入示例数据", "ok");
}

function formatJson() {
  try {
    els.json.value = JSON.stringify(JSON.parse(els.json.value), null, 2);
    setMsg("JSON 已格式化", "ok");
  } catch (e) {
    setMsg("JSON 解析失败：" + e.message, "err");
  }
}

/* ---------------- 历史任务 ---------------- */

async function loadJobs() {
  const body = await (await fetch("/api/jobs")).json();
  if (!body.jobs.length) {
    els.jobList.innerHTML = `<div class="empty" style="padding:10px 0">暂无历史任务。</div>`;
    return;
  }
  els.jobList.innerHTML = body.jobs
    .map((j) => {
      const risk = (j.high_risk_hairpins + j.high_risk_dimers > 0)
        ? `<span class="badge high">高风险 ${j.high_risk_hairpins + j.high_risk_dimers}</span>`
        : `<span class="badge low">通过</span>`;
      return `<div class="job-row">
        <span class="jid" data-id="${j.job_id}" title="点击查看完整结果">${j.job_id}</span>
        ${risk}
        <span class="jname">${esc(j.name || "（未命名任务）")}</span>
        <span class="jmeta">
          <span>${j.created_at.replace("T", " ").replace("+00:00", " UTC")}</span>
          <span>${j.primer_count} 条引物</span>
          <span>[Na⁺] ${j.salt_conc} M</span>
          <span>Tm 阈值 ${j.temperature}°C</span>
          <span>事件 ${j.events_total}</span>
          <span>淘汰 ${j.discard.length}</span>
        </span>
        <button class="btn tiny danger" data-del="${j.job_id}" type="button">删除</button>
      </div>`;
    })
    .join("");

  els.jobList.querySelectorAll(".jid").forEach((node) =>
    node.addEventListener("click", () => openJob(node.dataset.id)));
  els.jobList.querySelectorAll("[data-del]").forEach((node) =>
    node.addEventListener("click", () => removeJob(node.dataset.del)));
}

async function openJob(jobId) {
  const job = await (await fetch(`/api/jobs/${jobId}`)).json();
  const req = job.request;
  els.jobName.value = job.name || "";
  els.salt.value = req.salt_conc;
  els.temp.value = req.temperature;
  els.minbp.value = req.min_bp;
  els.json.value = JSON.stringify(req.primers, null, 2);
  renderResult({ job_id: job.job_id, result: job.result });
  setMsg(`已载入历史任务 ${job.job_id}`, "ok");
  window.scrollTo({ top: 0, behavior: "smooth" });
}

async function removeJob(jobId) {
  if (!confirm(`确定删除任务 ${jobId}？`)) return;
  const resp = await fetch(`/api/jobs/${jobId}`, { method: "DELETE" });
  if (resp.ok) {
    setMsg(`已删除 ${jobId}`, "ok");
    loadJobs();
  }
}

/* ---------------- 条件敏感性分析 ---------------- */

const SENS_COLORS = {
  critical: "255,92,122",
  high: "255,159,67",
  medium: "255,209,102",
  low: "77,213,153",
};

function resetSensitivity() {
  els.sensBox.hidden = !lastResponse || !lastResponse.job_id;
  lastSensitivity = null;
  els.sensResult.hidden = true;
  els.sensResult.innerHTML = "";
  els.sensMsg.textContent = "";
  els.sensMsg.className = "msg";
  els.sensExport.disabled = true;
}

function setSensMsg(text, kind = "") {
  els.sensMsg.textContent = text;
  els.sensMsg.className = "msg" + (kind ? " " + kind : "");
}

/* 解析逗号/空白分隔的数值列表。 */
function parseNumList(raw) {
  const parts = String(raw).split(/[\s,，;；]+/).filter((x) => x !== "");
  const out = [];
  for (const p of parts) {
    const v = Number(p);
    if (!Number.isFinite(v)) throw new Error(`无法解析数值：${p}`);
    out.push(v);
  }
  return out;
}

function fmtNum(v) {
  return Math.abs(v - Math.round(v)) < 1e-9 ? String(Math.round(v)) : String(v);
}

/* 根据格内风险结构给出背景色：优先 critical，其次 high，再次 medium。 */
function cellStyle(cell) {
  if (cell.critical_total > 0) {
    const a = Math.min(0.95, 0.22 + 0.12 * Math.min(cell.critical_total, 6));
    return `background:rgba(${SENS_COLORS.critical},${a.toFixed(2)});color:#fff`;
  }
  if (cell.high_risk_total > 0) {
    const a = Math.min(0.8, 0.16 + 0.1 * Math.min(cell.high_risk_total, 6));
    return `background:rgba(${SENS_COLORS.high},${a.toFixed(2)});color:#fff`;
  }
  if (cell.events_reported > 0) {
    const a = Math.min(0.55, 0.12 + 0.06 * Math.min(cell.events_reported, 6));
    return `background:rgba(${SENS_COLORS.medium},${a.toFixed(2)});color:#1a1a1a`;
  }
  return "background:rgba(77,213,153,.12);color:var(--muted)";
}

function renderSensitivity(data) {
  lastSensitivity = data;
  const s = data.sensitivity;
  els.sensResult.hidden = false;
  els.sensExport.disabled = false;

  const salts = s.salt_concs;
  const temps = s.temperatures;
  const cellAt = (si, ti) => s.cells[si * temps.length + ti];
  const fc = s.first_critical;

  const head = `<tr><th class="corner">盐浓度＼Tm 阈值</th>${
    temps.map((t) => `<th>${fmtNum(t)}°C</th>`).join("")
  }</tr>`;

  const rows = salts.map((salt, si) => {
    const tds = temps.map((_, ti) => {
      const c = cellAt(si, ti);
      const isFirst = fc && fc.salt_index === si && fc.temp_index === ti;
      const isBase = c.is_baseline;
      return `<td class="sens-cell${isFirst ? " first-crit" : ""}${isBase ? " baseline" : ""}"
        style="${cellStyle(c)}" data-si="${si}" data-ti="${ti}" tabindex="0"
        title="点击查看该格详情">
        <span class="cell-num">${c.high_risk_total}</span>
        ${c.critical_total > 0 ? `<span class="cell-crit">⛔${c.critical_total}</span>` : ""}
        <span class="cell-discard">🚫${c.discard_count}</span>
        ${isFirst ? '<span class="first-star" title="首次达到 critical 的条件">★</span>' : ""}
      </td>`;
    }).join("");
    return `<tr><th class="rowhead">${fmtNum(salt)} M</th>${tds}</tr>`;
  }).join("");

  const firstText = fc
    ? `⭐ <b>首次达到 critical</b>：盐浓度 <b>${fmtNum(fc.salt_conc)} M</b>、Tm 阈值 <b>${fmtNum(fc.temperature)}°C</b>（critical ${fc.critical_total} / 高风险 ${fc.high_risk_total}）。`
      + `<br><span class="hint" style="margin:0">判定顺序：从最低盐浓度行起（行内自最高 Tm 阈值向左），第一个出现 critical 事件的格子。</span>`
    : `<span class="safe-text">✅ 所有条件组合下均未出现 critical 事件。</span>`;

  els.sensResult.innerHTML = `
    <div class="sens-legend">
      <span><i class="sw" style="background:rgb(${SENS_COLORS.critical})"></i>critical</span>
      <span><i class="sw" style="background:rgb(${SENS_COLORS.high})"></i>high</span>
      <span><i class="sw" style="background:rgb(${SENS_COLORS.medium})"></i>medium</span>
      <span><i class="sw sw-pass"></i>无中高风险</span>
      <span><b>★</b> 首次 critical</span>
      <span class="base-swatch">加粗边框 = 任务原始条件（baseline）</span>
    </div>
    <p class="sens-first">${firstText}</p>
    <p class="hint">每格格式：<b>高风险数</b>（high+critical）／<b>⛔critical 数</b>／<b>🚫淘汰引物数</b>。点击格子查看淘汰清单与相对 baseline 的风险变化。</p>
    <div class="heat-scroll">
      <table class="heatmap"><thead>${head}</thead><tbody>${rows}</tbody></table>
    </div>
    <div id="sensCellDetail" class="sens-detail"><div class="empty" style="padding:14px 0">点击任意格子查看详情。</div></div>
    <p class="hint">对齐事件明细共 ${s.events_total} 条${s.events_capped ? `（JSON 中仅保留风险最高的 ${s.events_shown} 条；聚合计数不受影响）` : "（全部包含在下载 JSON 中）"}，每条附完整风险矩阵（行=盐浓度、列=Tm 阈值）与各自的首次 critical 坐标。</p>
  `;

  els.sensResult.querySelectorAll(".sens-cell").forEach((node) => {
    node.addEventListener("click", () =>
      renderSensCellDetail(cellAt(Number(node.dataset.si), Number(node.dataset.ti))));
  });
}

function changeListHtml(title, arr, count, capped, colorClass) {
  if (count === 0) return "";
  const items = arr.length
    ? arr.map((k) => `<li class="mono">${esc(k)}</li>`).join("")
    : "";
  return `<div class="change-col ${colorClass}">
    <h5>${title}（${count}${capped ? "，仅列前 " + arr.length : ""}）</h5>
    ${items ? `<ul>${items}</ul>` : `<div class="hint" style="margin:0">事件键过多，已省略列表，计数以括号为准。</div>`}
  </div>`;
}

function renderSensCellDetail(c) {
  const box = $("#sensCellDetail");
  const ch = c.changes_vs_baseline;
  box.innerHTML = `
    <h4>${fmtNum(c.salt_conc)} M × Tm 阈值 ${fmtNum(c.temperature)}°C
      ${c.is_baseline ? '<span class="badge tag">baseline</span>' : ""}</h4>
    <div class="detail-grid">
      <span>高风险（high+critical）<b>${c.high_risk_total}</b></span>
      <span>critical <b>${c.critical_total}</b>（发卡 ${c.critical_hairpins} / 二聚体 ${c.critical_dimers}）</span>
      <span>high <b>${c.high_risk_hairpins + c.high_risk_dimers}</b>（发卡 ${c.high_risk_hairpins} / 二聚体 ${c.high_risk_dimers}）</span>
      <span>medium <b>${c.medium_hairpins + c.medium_dimers}</b>（发卡 ${c.medium_hairpins} / 二聚体 ${c.medium_dimers}）</span>
      <span>进入报告事件 <b>${c.events_reported}</b></span>
      <span>淘汰引物 <b>${c.discard_count}</b></span>
    </div>
    <div class="discard-line">🚫 ${c.discard.length ? c.discard.map(esc).join("，") : "（无）"}</div>
    ${c.is_baseline ? `<p class="hint" style="margin:8px 0 0">该格即任务原始条件，风险变化为空。</p>` : `
      <div class="change-vs-base">
        <h5>相对 baseline（${fmtNum(lastSensitivity.sensitivity.baseline.salt_conc)} M / ${fmtNum(lastSensitivity.sensitivity.baseline.temperature)}°C）的风险变化</h5>
        <div class="change-cols">
          ${changeListHtml("新进入高风险的事件", ch.new_high, ch.new_high_count, ch.new_high_capped, "up")}
          ${changeListHtml("跌出高风险的事件", ch.resolved, ch.resolved_count, ch.resolved_capped, "down")}
        </div>
        <div class="discard-diff">
          <span class="up">新增淘汰：${ch.discard_added.length ? ch.discard_added.map(esc).join("，") : "无"}</span>
          <span class="down">解除淘汰：${ch.discard_removed.length ? ch.discard_removed.map(esc).join("，") : "无"}</span>
        </div>
      </div>`}
  `;
}

async function runSensitivity() {
  if (!lastResponse || !lastResponse.job_id) {
    setSensMsg("请先执行或载入一个已保存的历史任务。", "err");
    return;
  }
  let salts, temps;
  try {
    salts = parseNumList(els.sensSalts.value);
    temps = parseNumList(els.sensTemps.value);
  } catch (e) {
    setSensMsg(e.message, "err");
    return;
  }
  if (!salts.length || !temps.length) {
    setSensMsg("盐浓度列表与 Tm 阈值列表均不能为空（各不超过 8 项）。", "err");
    return;
  }
  if (salts.length * temps.length > 64) {
    setSensMsg(`组合为 ${salts.length * temps.length} 格，超过上限 64 格（各列表 ≤ 8 项）。`, "err");
    return;
  }
  els.sensRun.disabled = true;
  setSensMsg("逐格重算中…");
  try {
    const resp = await fetch(`/api/jobs/${lastResponse.job_id}/sensitivity`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ salt_concs: salts, temperatures: temps }),
    });
    const body = await resp.json();
    if (!resp.ok) {
      setSensMsg("请求失败（" + resp.status + "）：" + (body.detail || ""), "err");
      return;
    }
    setSensMsg("分析完成 ✓ 共 " + body.sensitivity.cell_count + " 格", "ok");
    renderSensitivity(body);
  } catch (e) {
    setSensMsg("网络错误：" + e.message, "err");
  } finally {
    els.sensRun.disabled = false;
  }
}

function exportSensitivityJson() {
  if (!lastSensitivity) return;
  const blob = new Blob([JSON.stringify(lastSensitivity, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `primer-sensitivity-${lastSensitivity.job_id}.json`;
  a.click();
  URL.revokeObjectURL(a.href);
}

/* ---------------- 多重 PCR 分池 ---------------- */

function resetPooling() {
  els.poolBox.hidden = !lastResponse || !lastResponse.job_id;
  lastPooling = null;
  poolPrimers = [];
  poolChecked = new Set();
  els.poolPrimerList.innerHTML = "";
  els.poolPairs.value = "";
  els.poolMax.value = "4";
  els.poolResult.hidden = true;
  els.poolResult.innerHTML = "";
  setPoolMsg("", "");
  els.poolExport.disabled = true;
}

function setPoolMsg(text, kind = "") {
  els.poolMsg.textContent = text;
  els.poolMsg.className = "msg" + (kind ? " " + kind : "");
}

/* 从当前表单 JSON 解析 id->sequence；解析失败时退回结果页 discard 信息。 */
function currentPrimerMap() {
  try {
    const parsed = JSON.parse(els.json.value);
    return new Map(normalizePrimers(parsed).map((p) => [p.id, p.sequence]));
  } catch (e) {
    return new Map();
  }
}

function renderPoolPrimerList() {
  const seqMap = currentPrimerMap();
  const worst = new Map();
  const r = lastResponse && lastResponse.result;
  if (r) {
    for (const p of r.summary.primers) {
      if (p.worst_risk !== "low") worst.set(p.primer_id, p.worst_risk);
    }
  }
  poolPrimers = Array.from(seqMap.keys())
    .filter((id) => r && r.summary.primers.some((p) => p.primer_id === id))
    .map((id) => {
      const p = r.summary.primers.find((x) => x.primer_id === id);
      return { id, sequence: seqMap.get(id) || "", risk: (p && p.worst_risk) || "low" };
    });
  if (!poolPrimers.length) {
    els.poolPrimerList.innerHTML = `<div class="empty" style="padding:8px 0">无法读取该任务的引物列表，请重新载入任务。</div>`;
    return;
  }
  els.poolPrimerList.innerHTML = poolPrimers.map((p) => {
    const checked = poolChecked.has(p.id) ? "checked" : "";
    const riskTag = p.risk !== "low" ? badge(RISK_CN[p.risk] + "风险", p.risk) : "";
    return `<label class="pool-primer" title="${esc(p.sequence || "序列不可用，请重新载入任务")}">
      <input type="checkbox" data-pid="${esc(p.id)}" ${checked} />
      <span class="mono">${esc(p.id)}</span>${riskTag}
    </label>`;
  }).join("");
  els.poolPrimerList.querySelectorAll("input[type=checkbox]").forEach((cb) => {
    cb.addEventListener("change", () => {
      const pid = cb.dataset.pid;
      if (cb.checked) {
        poolChecked.add(pid);
        return;
      }
      poolChecked.delete(pid);
      // 去掉包含该引物的必须同池对（逐行解析，用户正在输入的非法行保留）。
      removePairLineUsing(pid);
    });
  });
}

/* 取消勾选引物时，直接从文本框中删掉引用它的行，保留其余（含非法）行。 */
function removePairLineUsing(pid) {
  const kept = els.poolPairs.value
    .split("\n")
    .filter((line) => line.trim())
    .filter((line) => {
      try {
        const ids = parsePoolPairs(line)[0];
        return !ids.includes(pid);
      } catch (e) {
        return true;
      }
    });
  els.poolPairs.value = kept.join("\n");
}

/* 解析 textarea：每行一对，分隔符为逗号/空白（兼容 JSON 风格 ["A","B"] 等写法）。 */
function parsePoolPairs(raw) {
  const pairs = [];
  const lines = String(raw).split(/\n+/).map((l) => l.trim()).filter(Boolean);
  for (const line of lines) {
    const ids = line.replace(/[\[\]"'`]/g, "").split(/[\s,，;；/]+/).filter(Boolean);
    if (ids.length !== 2) {
      throw new Error(`无法解析必须同池对（每行需恰好 2 个引物 id）：${line}`);
    }
    pairs.push(ids);
  }
  return pairs;
}

function poolSiteRow(site) {
  const anc = site.anchored_three_prime ? badge("3′ 锚定 · 可延伸", "critical") : "";
  return `<div class="pool-site">
    <div class="event-title">
      ${badge(RISK_CN[site.risk] + "风险", site.risk)}
      ${anc}
      <span class="who mono">${esc(site.primer_a)}</span>
      <span style="color:var(--muted)">×</span>
      <span class="who mono">${esc(site.primer_b)}</span>
    </div>
    <div class="event-meta">
      <span>${esc(site.primer_a)} 区间 <b>${site.interval_a[0]}–${site.interval_a[1]}</b></span>
      <span>${esc(site.primer_b)} 区间 <b>${site.interval_b[0]}–${site.interval_b[1]}</b></span>
      <span>配对 <b>${site.bp}</b> bp</span>
      <span>GC <b>${site.gc_pairs}</b></span>
      <span>Tm <b>${site.tm}°C</b> <span style="opacity:.7">(无盐 ${site.tm_base}°C)</span></span>
    </div>
    ${alignmentBlock(site.alignment)}
  </div>`;
}

function renderPooling(data) {
  lastPooling = data;
  const p = data.pooling;
  els.poolResult.hidden = false;
  els.poolExport.disabled = false;

  if (!p.feasible) {
    els.poolResult.innerHTML = `
      <div class="pool-infeasible">
        <h4>⛔ 在上限 ${p.max_pools} 池内无法完成分池</h4>
        <div>${esc(p.infeasible_reason || "存在无法消解的冲突。")}</div>
        <p class="hint" style="margin:8px 0 0">不会输出无效方案。可调整必须同池配对、提高池数上限，或回到筛选结果更换引物后重算。</p>
      </div>
      ${renderPoolConflicts(p.conflicts, false)}`;
    return;
  }

  const unitNote = p.units.some((u) => u.forced)
    ? `<p class="hint" style="margin:6px 0 0">标签 U1、U2… 为按必须同池对折叠后的单元；池内按勾选顺序列出引物。</p>`
    : "";
  const pools = p.pools.map((pool) => `
    <div class="pool-card">
      <h4>池 ${pool.pool} <span class="hint" style="margin:0">（单元 ${pool.unit_ids.map(esc).join("、")}，${pool.members.length} 条）</span></h4>
      <ul class="pool-members">${pool.members.map((m) => `<li class="mono">${esc(m)}</li>`).join("")}</ul>
    </div>`).join("");

  els.poolResult.innerHTML = `
    <div class="pool-ok">
      <h4>✅ 最少需要 ${p.pool_count} 个池（上限 ${p.max_pools}）</h4>
      <div class="pool-cards">${pools}</div>
      ${unitNote}
    </div>
    ${renderPoolConflicts(p.conflicts, true)}`;
}

function renderPoolConflicts(conflicts, feasible) {
  if (!conflicts.length) {
    return `<p class="hint" style="margin-top:10px">所选引物两两之间均无 high/critical 异源二聚体位点，全部引物可同池。</p>`;
  }
  const head = feasible
    ? `<h4>⚔️ 不可同池冲突位点（${conflicts.length} 对单元，已被分配到不同池）</h4>`
    : `<h4>⚔️ 冲突位点明细（${conflicts.length} 对单元）</h4>`;
  const items = conflicts.map((c) => {
    const tag = c.same_unit
      ? badge(`同单元 ${esc(c.unit_a)} · 强制同池冲突`, "critical")
      : badge(`${esc(c.unit_a)} → 池 ${c.pool_a} · ${esc(c.unit_b)} → 池 ${c.pool_b}`, "tag");
    return `<details class="pool-conflict ${c.same_unit ? "critical" : c.worst_risk}"${c.same_unit ? " open" : ""}>
      <summary>
        ${badge(RISK_CN[c.worst_risk] + "风险", c.worst_risk)}
        <span class="mono">${esc(c.primer_a)}</span>
        <span style="color:var(--muted)">×</span>
        <span class="mono">${esc(c.primer_b)}</span>
        ${tag}
        <span class="hint" style="margin:0">${c.site_count} 个位点</span>
      </summary>
      <div class="pool-site-list">${c.sites.map(poolSiteRow).join("")}</div>
    </details>`;
  }).join("");
  return `<div class="pool-conflicts">${head}<p class="hint">冲突边仅统计 high/critical 异源二聚体（发卡与自二聚体不影响分池）；位点均按任务保存条件全量重扫得到。</p>${items}</div>`;
}

async function runPooling() {
  if (!lastResponse || !lastResponse.job_id) {
    setPoolMsg("请先执行或载入一个已保存的历史任务。", "err");
    return;
  }
  const ids = Array.from(poolChecked);
  if (ids.length < 2 || ids.length > 12) {
    setPoolMsg(`请选择 2–12 条引物（当前选了 ${ids.length} 条）。`, "err");
    return;
  }
  let forcedPairs;
  try {
    forcedPairs = parsePoolPairs(els.poolPairs.value);
  } catch (e) {
    setPoolMsg(e.message, "err");
    return;
  }
  const idSet = new Set(ids);
  const badRef = forcedPairs.flat().find((x) => !idSet.has(x));
  if (badRef) {
    setPoolMsg(`必须同池对中的引物 ${badRef} 不在已勾选的引物中。`, "err");
    return;
  }
  const maxPools = parseInt(els.poolMax.value, 10);
  if (!(maxPools >= 2 && maxPools <= 6)) {
    setPoolMsg("池数上限需为 2–6 的整数。", "err");
    return;
  }
  // 勾选顺序 = 任务引物顺序，保证池内成员展示与选择顺序一致。
  const order = new Map(poolPrimers.map((p, i) => [p.id, i]));
  const orderedIds = ids.slice().sort((a, b) => order.get(a) - order.get(b));
  els.poolRun.disabled = true;
  setPoolMsg("全量重扫异源二聚体并搜索最少池数…");
  try {
    const resp = await fetch(`/api/jobs/${lastResponse.job_id}/pooling`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ primer_ids: orderedIds, forced_pairs: forcedPairs, max_pools: maxPools }),
    });
    const body = await resp.json();
    if (!resp.ok) {
      setPoolMsg("请求失败（" + resp.status + "）：" + (body.detail || ""), "err");
      return;
    }
    const p = body.pooling;
    setPoolMsg(
      p.feasible
        ? `分池完成 ✓ 最少 ${p.pool_count} 池（上限 ${p.max_pools}），冲突边 ${p.conflicts.length} 对。`
        : `无法在 ${p.max_pools} 池内分池：已给出原因，未输出无效方案。`,
      p.feasible ? "ok" : "err"
    );
    renderPooling(body);
  } catch (e) {
    setPoolMsg("网络错误：" + e.message, "err");
  } finally {
    els.poolRun.disabled = false;
  }
}

function exportPoolingJson() {
  if (!lastPooling) return;
  const blob = new Blob([JSON.stringify(lastPooling, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `primer-pooling-${lastPooling.job_id}.json`;
  a.click();
  URL.revokeObjectURL(a.href);
}

/* ---------------- 绑定 ---------------- */

els.run.addEventListener("click", runScreen);
els.exportBtn.addEventListener("click", exportJson);
els.example.addEventListener("click", loadExample);
els.format.addEventListener("click", formatJson);
els.refreshJobs.addEventListener("click", loadJobs);
els.sensRun.addEventListener("click", runSensitivity);
els.sensExport.addEventListener("click", exportSensitivityJson);
els.poolRun.addEventListener("click", runPooling);
els.poolExport.addEventListener("click", exportPoolingJson);

loadJobs();
