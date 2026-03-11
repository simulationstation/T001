from __future__ import annotations

import html
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
from astropy.time import Time
from PIL import Image

from .utils import ensure_dir, env_like, git_like, json_list, utc_stamp, write_json, write_progress, write_text


CLASS_ORDER = [
    "STRONG_EXPORT_FAILURE_LIKE",
    "INTERMEDIATE_BRANCH_LIKE",
    "DUST_SURVIVOR_LIKE",
    "SYSTEMATIC_LIKE",
]

CLASS_DESCRIPTIONS = {
    "STRONG_EXPORT_FAILURE_LIKE": "Best current read: a branch-selected export-failure / underluminous-collapse lead that survives the strict detector controls.",
    "INTERMEDIATE_BRANCH_LIKE": "Best current read: a partial-suppression or fallback-heavy lead that looks more collapse-like than ordinary variability, but is not clean enough to call strong.",
    "DUST_SURVIVOR_LIKE": "Best current read: a real fade competitor whose behavior is better explained by obscuration or survival than by true export failure.",
    "SYSTEMATIC_LIKE": "Best current read: subtraction, registration, or coverage systematics dominate the signal.",
}

CLASS_BADGES = {
    "STRONG_EXPORT_FAILURE_LIKE": "Strong export-failure-like",
    "INTERMEDIATE_BRANCH_LIKE": "Intermediate branch-like",
    "DUST_SURVIVOR_LIKE": "Dust-survivor-like",
    "SYSTEMATIC_LIKE": "Systematic-like",
}

CLASS_COLORS = {
    "STRONG_EXPORT_FAILURE_LIKE": "#d96f32",
    "INTERMEDIATE_BRANCH_LIKE": "#c79c2e",
    "DUST_SURVIVOR_LIKE": "#5f7ea3",
    "SYSTEMATIC_LIKE": "#7c536f",
}


def _as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, tuple):
        return [str(item) for item in value]
    if hasattr(value, "tolist"):
        maybe = value.tolist()
        if isinstance(maybe, list):
            return [str(item) for item in maybe]
    text = str(value).strip()
    if not text:
        return []
    try:
        payload = json.loads(text)
    except Exception:
        return [text]
    if isinstance(payload, list):
        return [str(item) for item in payload]
    return [str(payload)]


def _fmt_float(value: Any, digits: int = 3) -> str:
    try:
        num = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(num):
        return "n/a"
    return f"{num:.{digits}f}"


def _fmt_days(value: Any) -> str:
    try:
        num = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(num):
        return "n/a"
    return f"{num:.1f} d"


def _fmt_pct(value: Any) -> str:
    try:
        num = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(num):
        return "n/a"
    return f"{100.0 * num:.1f}%"


def _fmt_date_from_mjd(value: Any) -> str:
    try:
        num = float(value)
    except Exception:
        return "n/a"
    if not math.isfinite(num):
        return "n/a"
    try:
        return Time(num, format="mjd").utc.iso.split()[0]
    except Exception:
        return "n/a"


def _representative_detection(frame: pd.DataFrame) -> pd.Series:
    if frame.empty:
        raise ValueError("cluster has no detections")
    rank = frame.copy()
    status_order = {"PASS": 2, "REVIEW": 1, "BENCHMARK_SIGNAL": 1}
    rank["status_order"] = rank["status"].astype(str).map(status_order).fillna(0)
    rank = rank.sort_values(
        [
            "status_order",
            "diff_sigma",
            "pre_snr",
            "branch_rank_score",
            "registration_residual_px",
        ],
        ascending=[False, False, False, False, True],
    ).reset_index(drop=True)
    return rank.iloc[0]


def _filter_summary(frame: pd.DataFrame) -> str:
    filters: list[str] = []
    for col in ("pre_filter", "post_filter"):
        if col not in frame.columns:
            continue
        for value in frame[col].dropna().astype(str):
            if value and value not in filters:
                filters.append(value)
    return ", ".join(filters) if filters else "n/a"


def _instrument_summary(frame: pd.DataFrame) -> str:
    instruments: list[str] = []
    for col in ("pre_instrument", "post_instrument"):
        if col not in frame.columns:
            continue
        for value in frame[col].dropna().astype(str):
            if value and value not in instruments:
                instruments.append(value)
    return ", ".join(instruments) if instruments else "n/a"


def _observation_window(frame: pd.DataFrame) -> str:
    pre_min = pd.to_numeric(frame["pre_mjd"], errors="coerce").min()
    post_max = pd.to_numeric(frame["post_mjd"], errors="coerce").max()
    return f"{_fmt_date_from_mjd(pre_min)} to {_fmt_date_from_mjd(post_max)}"


def _build_gif_and_panel(
    *,
    cutout_path: Path,
    gif_path: Path,
    panel_path: Path,
    diff_path: Path,
) -> None:
    image = Image.open(cutout_path).convert("RGB")
    width, height = image.size
    x1 = round(width / 3.0)
    x2 = round(2.0 * width / 3.0)

    panel_preview = image.copy()
    panel_preview.thumbnail((1080, 1080), Image.Resampling.LANCZOS)
    panel_preview.save(panel_path, optimize=True)

    # Crop inside the plotted region so the blink does not jump because of panel labels.
    margin_x = max(12, x1 // 20)
    margin_top = max(24, height // 12)
    margin_bottom = max(22, height // 10)
    pre = image.crop((0 + margin_x, margin_top, x1 - margin_x, height - margin_bottom))
    post = image.crop((x1 + margin_x, margin_top, x2 - margin_x, height - margin_bottom))
    diff = image.crop((x2 + margin_x, margin_top, width - margin_x, height - margin_bottom))

    target_height = 340
    target_width = max(200, round(pre.size[0] * (target_height / max(pre.size[1], 1))))
    pre_small = pre.resize((target_width, target_height), Image.Resampling.LANCZOS)
    post_small = post.resize((target_width, target_height), Image.Resampling.LANCZOS)
    diff_small = diff.resize((target_width, target_height), Image.Resampling.LANCZOS)

    diff_small.save(diff_path, optimize=True)
    pre_small.save(
        gif_path,
        save_all=True,
        append_images=[post_small, pre_small, post_small],
        duration=[1800, 1800, 1800, 1800],
        loop=0,
        optimize=True,
        disposal=2,
    )


def _candidate_interpretation(cluster: pd.Series) -> str:
    branch_class = str(cluster["branch_class"])
    base = CLASS_DESCRIPTIONS.get(branch_class, "No interpretation available.")
    support = int(cluster.get("n_support_pairs", 0) or 0)
    late_return = float(cluster.get("forced_late_return_fraction", 0.0) or 0.0)
    dust = float(cluster.get("dust_survivor_score", 0.0) or 0.0)
    systematic = float(cluster.get("systematic_risk", 0.0) or 0.0)
    if branch_class == "STRONG_EXPORT_FAILURE_LIKE" and support >= 4 and late_return <= 0.1:
        return base + " Multi-pair same-location support is already present, and late-time return is weak."
    if branch_class == "INTERMEDIATE_BRANCH_LIKE" and dust < 0.25 and systematic < 0.2:
        return base + " The signal is coherent enough to treat as a real weak-branch lead rather than just an interesting fade."
    if branch_class == "DUST_SURVIVOR_LIKE":
        return base + " The fade is still interesting, but the current scoring says collapse-like export failure is not the clean leading explanation."
    return base


def _top_math_rows(cluster: pd.Series) -> list[tuple[str, str]]:
    return [
        ("Branch rank", f"{int(cluster['branch_rank'])} / {int(cluster['n_clusters_total'])}"),
        ("Rank score", _fmt_float(cluster["branch_rank_score"])),
        ("Export-failure score", _fmt_float(cluster["export_failure_score"])),
        ("Intermediate-branch score", _fmt_float(cluster["partial_branch_score"])),
        ("Dust-survivor score", _fmt_float(cluster["dust_survivor_score"])),
        ("Systematic risk", _fmt_float(cluster["systematic_risk"])),
        ("Detector confidence", _fmt_float(cluster["detector_confidence"])),
        ("Registration confidence", _fmt_float(cluster["registration_confidence"])),
        ("Subtraction cleanliness", _fmt_float(cluster["subtraction_cleanliness"])),
        ("Coverage completeness", _fmt_float(cluster["coverage_completeness"])),
        ("Forced-photometry coherence", _fmt_float(cluster["forced_photometry_coherence"])),
        ("Fade persistence", _fmt_float(cluster["fade_persistence"])),
        ("Rebrightening penalty", _fmt_float(cluster["rebrightening_penalty"])),
        ("Artifact penalty", _fmt_float(cluster["artifact_penalty"])),
        ("Forced ratio median R", _fmt_float(cluster["forced_ratio_median"])),
        ("Forced ratio IQR", _fmt_float(cluster["forced_ratio_iqr"])),
        ("Max diff significance", _fmt_float(cluster["forced_diff_sigma_max"])),
        ("Forced baseline span", _fmt_days(cluster["forced_baseline_span_days"])),
        ("Support pairs", str(int(cluster["n_support_pairs"]))),
        ("Support filters", str(int(cluster["n_support_filters"]))),
        ("Members in cluster", str(int(cluster["n_members"]))),
        ("PASS members", str(int(cluster["n_pass_members"]))),
        ("Forced measured pairs", f"{int(cluster['forced_n_pairs_measured'])} / {int(cluster['forced_n_pairs_considered'])}"),
        ("Partial fraction", _fmt_pct(cluster["forced_partial_fraction"])),
        ("Severe fraction", _fmt_pct(cluster["forced_severe_fraction"])),
        ("Return fraction", _fmt_pct(cluster["forced_return_fraction"])),
        ("Late-fade fraction", _fmt_pct(cluster["forced_late_fade_fraction"])),
        ("Late-return fraction", _fmt_pct(cluster["forced_late_return_fraction"])),
        ("Median registration residual", _fmt_float(cluster["forced_registration_residual_median"]) + " px"),
        ("Centroid scatter", _fmt_float(cluster["member_centroid_scatter_arcsec"]) + " arcsec"),
        ("Pair stable fraction", _fmt_pct(cluster["forced_pair_stable_fraction_median"])),
    ]


def _metrics_table(rows: list[tuple[str, str]]) -> str:
    items = []
    for label, value in rows:
        items.append(
            f"<tr><th>{html.escape(label)}</th><td>{html.escape(value)}</td></tr>"
        )
    return "\n".join(items)


def _card_html(candidate: dict[str, Any]) -> str:
    metrics = _metrics_table(candidate["math_rows"])
    reasons = "".join(f"<li><code>{html.escape(reason)}</code></li>" for reason in candidate["reason_codes"])
    warnings = "".join(f"<li><code>{html.escape(flag)}</code></li>" for flag in candidate["warning_flags"])
    warnings_block = warnings or "<li>None</li>"
    return f"""
<article class="candidate-card" data-class="{html.escape(candidate['branch_class'])}">
  <div class="card-top">
    <div class="card-heading">
      <div class="badge" style="--badge-color:{html.escape(candidate['badge_color'])}">{html.escape(candidate['badge'])}</div>
      <h3 id="{html.escape(candidate['cluster_id'])}">{html.escape(candidate['cluster_id'])}</h3>
      <p class="subhead">{html.escape(candidate['galaxy_name'])} | rank {candidate['branch_rank']} | confidence {html.escape(candidate['branch_confidence'])}</p>
      <p class="interpretation">{html.escape(candidate['interpretation'])}</p>
    </div>
    <div class="hero">
      <img loading="lazy" src="{html.escape(candidate['gif_rel'])}" alt="{html.escape(candidate['cluster_id'])} registered blink">
    </div>
  </div>
  <div class="panel-row">
    <img loading="lazy" class="triptych" src="{html.escape(candidate['panel_rel'])}" alt="{html.escape(candidate['cluster_id'])} pre/post/diff triptych">
    <img loading="lazy" class="diff-panel" src="{html.escape(candidate['diff_rel'])}" alt="{html.escape(candidate['cluster_id'])} difference panel">
  </div>
  <div class="mini-grid">
    <div>
      <h4>Observation context</h4>
      <ul>
        <li>Representative pair: <code>{html.escape(candidate['representative_pair'])}</code></li>
        <li>Filters: <code>{html.escape(candidate['filters'])}</code></li>
        <li>Instruments: <code>{html.escape(candidate['instruments'])}</code></li>
        <li>Date window: <code>{html.escape(candidate['date_window'])}</code></li>
        <li>Sky position: <code>RA {html.escape(candidate['ra_deg'])}</code>, <code>Dec {html.escape(candidate['dec_deg'])}</code></li>
      </ul>
    </div>
    <div>
      <h4>What this most likely is</h4>
      <p>{html.escape(candidate['interpretation'])}</p>
      <p class="small">This page is ranking branch-aware export-failure coherence under strict detector controls. It is not a claim of confirmed failed-supernova status.</p>
    </div>
  </div>
  <div class="table-wrap">
    <table>
      <tbody>
        {metrics}
      </tbody>
    </table>
  </div>
  <div class="details-grid">
    <div>
      <h4>Reason codes</h4>
      <ul>{reasons}</ul>
    </div>
    <div>
      <h4>Warning flags</h4>
      <ul>{warnings_block}</ul>
    </div>
  </div>
</article>
"""


def _style_text() -> str:
    return """
:root {
  --bg: #f3efe6;
  --paper: #fffdf8;
  --ink: #1d1b18;
  --muted: #6d6258;
  --line: #d8cfc0;
  --accent: #8c4b2f;
  --accent-soft: #c96e3a;
  --card-shadow: 0 18px 45px rgba(76, 52, 33, 0.12);
  --mono: "IBM Plex Mono", "SFMono-Regular", ui-monospace, monospace;
  --sans: "Source Serif 4", Georgia, serif;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: var(--sans);
  color: var(--ink);
  background:
    radial-gradient(circle at top left, rgba(210, 150, 94, 0.18), transparent 32%),
    linear-gradient(180deg, #f6f1e8 0%, #efe7d8 100%);
}
a { color: inherit; }
.wrap {
  max-width: 1420px;
  margin: 0 auto;
  padding: 40px 24px 96px;
}
.hero-block {
  border: 1px solid var(--line);
  border-radius: 24px;
  background: rgba(255,255,255,0.72);
  padding: 28px;
  box-shadow: var(--card-shadow);
  margin-bottom: 28px;
}
.eyebrow {
  font-family: var(--mono);
  color: var(--accent);
  font-size: 13px;
  letter-spacing: 0.08em;
  text-transform: uppercase;
}
h1 {
  font-size: clamp(2.2rem, 4vw, 4rem);
  line-height: 1;
  margin: 10px 0 14px;
}
.lede {
  max-width: 900px;
  font-size: 1.08rem;
  line-height: 1.55;
}
.top-grid,
.math-grid,
.section-grid,
.details-grid,
.mini-grid {
  display: grid;
  gap: 18px;
}
.top-grid { grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); margin-top: 22px; }
.math-grid { grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); margin-top: 24px; }
.section-grid { grid-template-columns: 1fr; gap: 20px; }
.mini-grid { grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); margin: 18px 0; }
.details-grid { grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); }
.stat, .math-box, .candidate-card {
  border: 1px solid var(--line);
  border-radius: 22px;
  background: var(--paper);
  box-shadow: var(--card-shadow);
}
.stat, .math-box { padding: 20px; }
.stat strong {
  display: block;
  font-size: 2rem;
  margin-bottom: 6px;
}
.math-box h3, .candidate-card h3, .section h2, .mini-grid h4, .details-grid h4 {
  margin: 0 0 10px;
}
.math-box p, .section p, .candidate-card p, .candidate-card li { line-height: 1.45; }
.section {
  margin-top: 38px;
}
.section h2 {
  font-size: clamp(1.6rem, 2vw, 2.4rem);
  margin-bottom: 6px;
}
.candidate-card {
  padding: 22px;
}
.card-top {
  display: grid;
  gap: 18px;
  grid-template-columns: minmax(260px, 1.4fr) minmax(220px, 360px);
  align-items: start;
}
.card-heading .subhead {
  margin: 6px 0 10px;
  color: var(--muted);
  font-family: var(--mono);
  font-size: 0.92rem;
}
.interpretation {
  font-size: 1rem;
  margin: 0;
}
.badge {
  display: inline-block;
  padding: 6px 10px;
  border-radius: 999px;
  background: color-mix(in srgb, var(--badge-color) 18%, white);
  border: 1px solid color-mix(in srgb, var(--badge-color) 55%, white);
  color: color-mix(in srgb, var(--badge-color) 75%, black);
  font-family: var(--mono);
  font-size: 0.8rem;
  margin-bottom: 12px;
}
.hero img {
  width: 100%;
  border-radius: 16px;
  border: 1px solid var(--line);
  background: #111;
}
.panel-row {
  display: grid;
  gap: 14px;
  grid-template-columns: minmax(320px, 1fr) minmax(180px, 260px);
  margin: 18px 0;
}
.triptych, .diff-panel {
  width: 100%;
  border-radius: 16px;
  border: 1px solid var(--line);
  background: #111;
}
.table-wrap {
  overflow-x: auto;
}
table {
  width: 100%;
  border-collapse: collapse;
}
th, td {
  border-top: 1px solid var(--line);
  padding: 9px 10px;
  text-align: left;
  vertical-align: top;
}
th {
  width: 36%;
  font-weight: 600;
}
code {
  font-family: var(--mono);
  font-size: 0.9em;
}
.small {
  color: var(--muted);
  font-size: 0.95rem;
}
.section-note {
  color: var(--muted);
  margin-bottom: 18px;
}
.toc {
  display: flex;
  gap: 10px;
  flex-wrap: wrap;
  margin-top: 18px;
}
.toc a {
  text-decoration: none;
  padding: 8px 12px;
  border-radius: 999px;
  border: 1px solid var(--line);
  background: rgba(255,255,255,0.72);
  font-family: var(--mono);
  font-size: 0.88rem;
}
footer {
  margin-top: 42px;
  color: var(--muted);
  font-size: 0.98rem;
}
@media (max-width: 960px) {
  .card-top,
  .panel-row {
    grid-template-columns: 1fr;
  }
}
"""


def _build_html(*, candidates: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    cards_by_class: dict[str, list[str]] = {key: [] for key in CLASS_ORDER}
    for candidate in candidates:
        cards_by_class.setdefault(candidate["branch_class"], []).append(_card_html(candidate))
    toc = "".join(
        f'<a href="#{branch_class.lower()}">{html.escape(CLASS_BADGES.get(branch_class, branch_class))} ({summary["class_counts"].get(branch_class, 0)})</a>'
        for branch_class in CLASS_ORDER
        if summary["class_counts"].get(branch_class, 0) > 0
    )
    sections = []
    for branch_class in CLASS_ORDER:
        count = summary["class_counts"].get(branch_class, 0)
        if count <= 0:
            continue
        sections.append(
            f"""
<section class="section" id="{branch_class.lower()}">
  <h2>{html.escape(CLASS_BADGES.get(branch_class, branch_class))}</h2>
  <p class="section-note">{html.escape(CLASS_DESCRIPTIONS.get(branch_class, ""))} Current count: {count}.</p>
  <div class="section-grid">
    {''.join(cards_by_class.get(branch_class, []))}
  </div>
</section>
"""
        )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Uber Nova Search</title>
  <style>{_style_text()}</style>
</head>
<body>
  <div class="wrap">
    <section class="hero-block">
      <div class="eyebrow">SUPERNOVA / hidden operator page</div>
      <h1>Uber Nova Search</h1>
      <p class="lede">
        This page is the strict branch-aware candidate board generated from the empirically registered difference-image rerun and the
        branch-aware export-failure scorer. It is intentionally not linked from the site homepage. The practical question is not
        simply whether a source faded, but whether the residual stays coherent enough, long enough, and cleanly enough to look like
        a branch-selected export-failure collapse rather than dust, ordinary variability, or subtraction damage.
      </p>
      <div class="top-grid">
        <div class="stat"><strong>{summary['n_clusters']}</strong> scored clusters</div>
        <div class="stat"><strong>{summary['n_candidates']}</strong> displayed candidates</div>
        <div class="stat"><strong>{summary['benchmark_recovery_fraction']:.2f}</strong> blind benchmark recovery fraction</div>
        <div class="stat"><strong>{summary['benchmark_median_sep_arcsec']:.2f}\"</strong> benchmark median localization error</div>
      </div>
      <div class="math-grid">
        <div class="math-box">
          <h3>What Ryan likely cares about here</h3>
          <p><code>R_med = median(F_post / F_pre)</code> from same-location forced photometry across scanned pairs, not just a single catalog hit.</p>
          <p><code>S_diff,max = max(F_diff / sigma_diff)</code> over the support set, with empirical star-based registration already baked into the strict rerun.</p>
          <p><code>epsilon_reg</code> is the median registration residual in pixels. High values hurt rank directly.</p>
        </div>
        <div class="math-box">
          <h3>Branch-aware logic</h3>
          <p>The ranking favors persistent, same-location suppression with low rebrightening, low artifact risk, high forced-photometry coherence, and clean provenance.</p>
          <p>It penalizes dipole-like behavior, centroid migration, weak coverage, chip-gap / edge ambiguity, and dust-survivor-like patterns.</p>
        </div>
        <div class="math-box">
          <h3>Interpretive posture</h3>
          <p>These are not publication claims. They are strict detector survivors ranked by export-failure coherence under the current assumption-conditioned framework.</p>
          <p>The standing guardrail is still the blind known-supernova benchmark, which passed at <code>{summary['benchmark_recovery_fraction']:.2f}</code>.</p>
        </div>
      </div>
      <div class="toc">{toc}</div>
    </section>
    {''.join(sections)}
    <footer>
      Built {html.escape(summary['created_utc'])} from <code>{html.escape(summary['source_branch_dir'])}</code>.
      Page source is tracked in the SUPERNOVA repo and deployed separately under <code>/uber_nova_search</code>.
    </footer>
  </div>
</body>
</html>
"""


def run_uber_nova_site(
    *,
    root_dir: Path,
    branch_output_dir: Path,
    output_dir: Path | None = None,
    include_classes: tuple[str, ...] = tuple(CLASS_ORDER),
) -> dict[str, Path]:
    branch_output_dir = branch_output_dir.resolve()
    out_dir = (output_dir or (root_dir / "site" / "uber_nova_search")).resolve()
    ensure_dir(out_dir)
    ensure_dir(out_dir / "assets" / "gifs")
    ensure_dir(out_dir / "assets" / "panels")
    ensure_dir(out_dir / "assets" / "data")

    write_json(out_dir / "config.json", {"branch_output_dir": str(branch_output_dir), "include_classes": list(include_classes)})
    write_json(out_dir / "env.json", env_like())
    write_json(out_dir / "git_like.json", git_like(root_dir))
    write_progress(out_dir / "progress.json", "build", "running", {"output_dir": str(out_dir)})

    cluster_df = pd.read_parquet(branch_output_dir / "branch_cluster_scores.parquet")
    detection_df = pd.read_parquet(branch_output_dir / "branch_detection_scores.parquet")
    branch_summary = json.loads((branch_output_dir / "branch_summary.json").read_text(encoding="utf-8"))

    cluster_df = cluster_df[cluster_df["branch_class"].astype(str).isin(include_classes)].copy()
    cluster_df = cluster_df.sort_values(["branch_rank", "branch_rank_score"], ascending=[True, False]).reset_index(drop=True)
    cluster_df["n_clusters_total"] = int(len(cluster_df))

    candidates: list[dict[str, Any]] = []
    for _, cluster in cluster_df.iterrows():
        cluster_id = str(cluster["cluster_id"])
        members = detection_df[detection_df["cluster_id"].astype(str).eq(cluster_id)].copy()
        if members.empty:
            continue
        rep = _representative_detection(members)
        cutout_path = Path(str(rep["cutout_path"]))
        if not cutout_path.exists():
            continue
        gif_rel = Path("assets") / "gifs" / f"{cluster_id}_blink.gif"
        panel_rel = Path("assets") / "panels" / f"{cluster_id}_panel.png"
        diff_rel = Path("assets") / "panels" / f"{cluster_id}_diff.png"
        _build_gif_and_panel(
            cutout_path=cutout_path,
            gif_path=out_dir / gif_rel,
            panel_path=out_dir / panel_rel,
            diff_path=out_dir / diff_rel,
        )

        reason_codes = _as_list(cluster.get("reason_codes_json"))
        warning_codes = sorted({flag for flags in members["warning_flags_json"] for flag in _as_list(flags)})
        candidate = {
            "cluster_id": cluster_id,
            "branch_class": str(cluster["branch_class"]),
            "badge": CLASS_BADGES.get(str(cluster["branch_class"]), str(cluster["branch_class"])),
            "badge_color": CLASS_COLORS.get(str(cluster["branch_class"]), "#8c4b2f"),
            "branch_confidence": str(cluster["branch_confidence"]),
            "branch_rank": int(cluster["branch_rank"]),
            "branch_rank_score": float(cluster["branch_rank_score"]),
            "galaxy_name": str(cluster["galaxy_name"]),
            "ra_deg": _fmt_float(cluster["ra_deg"], 5),
            "dec_deg": _fmt_float(cluster["dec_deg"], 5),
            "representative_pair": str(rep["pair_id"]),
            "filters": _filter_summary(members),
            "instruments": _instrument_summary(members),
            "date_window": _observation_window(members),
            "gif_rel": str(gif_rel),
            "panel_rel": str(panel_rel),
            "diff_rel": str(diff_rel),
            "reason_codes": reason_codes,
            "warning_flags": warning_codes,
            "interpretation": _candidate_interpretation(cluster),
        }
        cluster_with_total = cluster.copy()
        cluster_with_total["n_clusters_total"] = len(cluster_df)
        candidate["math_rows"] = _top_math_rows(cluster_with_total)
        candidates.append(candidate)

    class_counts = {label: int((cluster_df["branch_class"] == label).sum()) for label in include_classes}
    summary = {
        "created_utc": utc_stamp(),
        "source_branch_dir": str(branch_output_dir),
        "n_clusters": int(len(cluster_df)),
        "n_candidates": int(len(candidates)),
        "class_counts": class_counts,
        "benchmark_recovery_fraction": float(branch_summary.get("benchmark_recovery_fraction", 0.0) or 0.0),
        "benchmark_median_sep_arcsec": float(branch_summary.get("benchmark_median_recovered_sep_arcsec", float("nan")) or float("nan")),
    }
    payload = {"summary": summary, "candidates": candidates}
    write_json(out_dir / "assets" / "data" / "candidates.json", payload)
    write_json(out_dir / "summary.json", summary)
    write_text(out_dir / "index.html", _build_html(candidates=candidates, summary=summary))
    write_progress(out_dir / "progress.json", "build", "completed", {"n_candidates": len(candidates), "output_dir": str(out_dir)})
    return {
        "index_html": out_dir / "index.html",
        "summary": out_dir / "summary.json",
        "candidate_data": out_dir / "assets" / "data" / "candidates.json",
    }
