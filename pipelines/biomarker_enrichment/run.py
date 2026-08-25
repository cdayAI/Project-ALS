#!/usr/bin/env python3
"""
pipelines/biomarker_enrichment/run.py

H-008 (NfL-trajectory ML trial enrichment) -- PRO-ACT-style SIMULATION.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALL DATA GENERATED HERE IS SIMULATED. NO REAL PATIENT RECORDS ARE READ. !!
!! Parameters are calibrated to published summaries; see config.yaml PMIDs.!!
!! Real-data validation is PENDING: data/proact/APPLICATION_DRAFT.md       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

What it does
------------
1. Simulates a screening pool of ALS subjects with latent ALSFRS-R decline
   rates (log-normal, PRO-ACT-calibrated), baseline ALSFRS-R, a noisy 0-3 mo
   clinical slope, and simulated plasma NfL (baseline + early slope).
2. Compares five enrollment arms:
     none            - no enrichment
     noise           - rank on pure noise        (negative control)
     clinical_static - ridge ML on observed early ALSFRS-R slope (+ baseline)
     nfl_ml          - ridge ML on NfL trajectory + clinical features
     oracle          - ranks on TRUE progression rate (positive control)
   Each arm keeps the top enrich_fraction of the pool by predicted
   membership in the 'normal progressor' band.
3. Bootstrap power curves over a per-arm n grid for each treatment-effect
   scenario; interpolates n needed for 80% power; reports % sample-size
   reduction vs the unenriched arm and the screening burden.

Run:    python3 run.py                (from this directory)
Exit 0 = positive/negative controls passed. Exit 1 = controls failed
(rule 5 of AGENTS.md: no trusted output without passing positive control).
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

HERE = Path(__file__).resolve().parent


# --------------------------------------------------------------------------
# 1. Cohort simulation
# --------------------------------------------------------------------------
def simulate_cohort(cfg, rng):
    """Subject-level simulated screening pool. SIMULATED DATA."""
    p = cfg["progression"]
    n = cfg["n_pool"]

    lam = rng.lognormal(p["log_mean"], p["log_sd"], size=n)
    lam = np.clip(lam, p["min"], p["max"])

    b_cfg = cfg["baseline_alsfrsr"]
    r = b_cfg["corr_with_log_lambda"]
    z_lam = (np.log(lam) - np.log(lam).mean()) / np.log(lam).std()
    base = b_cfg["mean"] + b_cfg["sd"] * (
        r * z_lam + np.sqrt(1 - r**2) * rng.normal(size=n)
    )
    base = np.clip(base, 18.0, 46.0)

    # Observed early clinical slope from two noisy visits (month 0 vs month 3)
    meas_sd = cfg["measurement"]["alsfrsr_visit_noise_sd"]
    w = cfg["measurement"]["early_window_months"]
    v0 = base + rng.normal(0, meas_sd, n)
    v1 = base - w * lam + rng.normal(0, meas_sd, n)
    clin_slope = -(v0 - v1) / w  # positive = declining, noisy proxy of lambda

    # Simulated NfL trajectory
    nf = cfg["nfl"]
    log_nfl0 = nf["baseline_log_mean"] + nf["beta_baseline"] * np.log(lam) \
        + rng.normal(0, nf["eps_baseline_sd"], n)
    log_nfl_slope = nf["slope_gamma"] * np.log(lam) \
        + rng.normal(0, nf["eps_slope_sd"], n)

    df = pd.DataFrame({
        "lambda_true_pts_per_month": lam,
        "baseline_alsfrsr": base,
        "obs_clinical_slope_0_3mo": clin_slope,
        "log_nfl_baseline": log_nfl0,
        "log_nfl_early_slope": log_nfl_slope,
    })

    # Sanity gates on calibration (fail loudly rather than drift silently)
    r_nfl = stats.pearsonr(df.log_nfl_baseline, np.log(df.lambda_true_pts_per_month))[0]
    assert 0.55 <= r_nfl <= 0.75, f"NfL-baseline correlation off target: {r_nfl:.2f}"
    med = df.lambda_true_pts_per_month.median()
    assert 0.6 <= med <= 1.1, f"median progression off target: {med:.2f}"
    return df, {"corr_log_nfl_baseline_vs_log_lambda": round(float(r_nfl), 3),
                "median_progression_pts_per_month": round(float(med), 3)}


# --------------------------------------------------------------------------
# 2. Enrichment arms
# --------------------------------------------------------------------------
def ridge_predict(X_train, y_train, X_apply, lmb=1.0):
    """Closed-form ridge on standardized features; adds an unpenalized
    intercept column internally. numpy-only."""
    mu, sd = X_train.mean(0), X_train.std(0)
    sd[sd == 0] = 1.0
    Ztr = (X_train - mu) / sd
    Zap = (X_apply - mu) / sd
    Xtr = np.hstack([np.ones((len(Ztr), 1)), Ztr])
    Xap = np.hstack([np.ones((len(Zap), 1)), Zap])
    A = Xtr.T @ Xtr + lmb * np.eye(Xtr.shape[1])
    A[0, 0] -= lmb  # do not penalize the intercept
    coef = np.linalg.solve(A, Xtr.T @ y_train)
    return Xap @ coef


def band_score(pred_lambda, band):
    """Higher score = more confidently inside the target progressor band."""
    lo, hi = band
    center, half = (lo + hi) / 2, (hi - lo) / 2
    return -np.abs(pred_lambda - center) / half


def build_arm_scores(df, cfg):
    """Return dict arm -> selection score array (higher = selected first)."""
    rng = np.random.default_rng(cfg["seed"] + 777)
    band = cfg["trial"]["target_band"]
    q = cfg["trial"]["enrich_fraction"]
    y = df.lambda_true_pts_per_month.values

    scores = {}

    # none: everyone equally eligible -> selection is random subset
    scores["none"] = rng.uniform(size=len(df))

    # noise negative control: ranking carries zero information
    scores["noise"] = rng.normal(size=len(df))

    # oracle positive control: perfect knowledge of lambda
    scores["oracle"] = band_score(y, band)

    # learned arms: fit predictors of lambda on a training split, apply to all
    idx = rng.permutation(len(df))
    tr, ap = idx[: len(df) // 2], idx[len(df) // 2:]
    feats = {
        "clinical_static": ["obs_clinical_slope_0_3mo", "baseline_alsfrsr"],
        "nfl_ml": ["log_nfl_baseline", "log_nfl_early_slope",
                   "obs_clinical_slope_0_3mo", "baseline_alsfrsr"],
    }
    for arm, cols in feats.items():
        Xtr = df.loc[tr, cols].values
        pred = np.empty(len(df))
        pred[tr] = ridge_predict(Xtr, y[tr], df.loc[tr, cols].values)
        pred[ap] = ridge_predict(Xtr, y[tr], df.loc[ap, cols].values)
        # nested-model diagnostic (H-008 falsification criterion 3):
        # does early NfL slope add beyond baseline NfL alone?
        if arm == "nfl_ml":
            cols_base = ["log_nfl_baseline", "obs_clinical_slope_0_3mo", "baseline_alsfrsr"]
            rss1 = _rss(df, tr, cols_base, y)
            rss2 = _rss(df, tr, cols, y)
            n_tr = len(tr)
            f = ((rss1 - rss2) / 1) / (rss2 / (n_tr - len(cols) - 1))
            p_lrt = float(stats.f.sf(f, 1, n_tr - len(cols) - 1))
            scores.setdefault("_diagnostics", {})
            scores["_diagnostics"]["nfl_slope_beyond_baseline_F_pvalue"] = p_lrt
        scores[arm] = band_score(pred, band)

    diag = scores.pop("_diagnostics", {})

    # apply fraction-q selection mask per arm
    masks, kept = {}, {}
    for arm, s in scores.items():
        order = np.argsort(-s)
        m = np.zeros(len(df), dtype=bool)
        m[order[: int(q * len(df))]] = True
        masks[arm] = m
        kept[arm] = int(m.sum())
    return masks, kept, diag


def _rss(df, tr_idx, cols, y):
    X = np.hstack([np.ones((len(tr_idx), 1)), df.iloc[tr_idx][cols].values])
    beta, *_ = np.linalg.lstsq(X, y[tr_idx], rcond=None)
    resid = y[tr_idx] - X @ beta
    return float(resid @ resid)


# --------------------------------------------------------------------------
# 3. Trial simulation + bootstrap power curves
# --------------------------------------------------------------------------
def simulate_trials(outcomes_ctrl, outcomes_trt, alpha):
    """Vectorized Welch t-tests across bootstrap replicates."""
    diff = outcomes_trt.mean(axis=1) - outcomes_ctrl.mean(axis=1)
    v1 = outcomes_ctrl.var(axis=1, ddof=1)
    v2 = outcomes_trt.var(axis=1, ddof=1)
    n1, n2 = outcomes_ctrl.shape[1], outcomes_trt.shape[1]
    se = np.sqrt(v1 / n1 + v2 / n2)
    t = diff / se
    df_w = (v1 / n1 + v2 / n2) ** 2 / (
        (v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1)
    )
    p = 2 * stats.t.sf(np.abs(t), df_w)
    return (p < alpha).mean()


def power_and_n80(pool_deltas, mask, cfg, effect_pts_month, rng):
    """Bootstrap replicates at each grid n; interpolate n at target power."""
    tcfg = cfg["trial"]
    dur = tcfg["duration_months"]
    sel = pool_deltas[mask]           # control-arm deltas of eligible subjects
    # Treatment shrinks each subject's delta magnitude proportionally so that
    # the mean absolute benefit equals effect_pts_per_month * duration:
    mean_mag = -sel.mean() / dur      # mean |decline| pts/month of eligibles
    frac = max(0.05, min(0.9, effect_pts_month / mean_mag))
    sel_trt = sel * (1.0 - frac)
    noise_sd = tcfg["outcome_noise_sd_per_delta"]
    reps = tcfg["bootstrap_reps_per_n"]
    rows = []
    for n in tcfg["n_grid"]:
        idx_c = rng.integers(0, len(sel), size=(reps, n))
        idx_t = rng.integers(0, len(sel_trt), size=(reps, n))
        oc = sel[idx_c] + rng.normal(0, noise_sd, (reps, n))
        ot = sel_trt[idx_t] + rng.normal(0, noise_sd, (reps, n))
        power = simulate_trials(oc, ot, tcfg["alpha"])
        rows.append({"effect_pts_per_month": effect_pts_month,
                     "n_per_arm": n, "power": round(float(power), 4)})
    pc = pd.DataFrame(rows)
    g = pc.n_per_arm.values
    pw = pc.power.values
    if pw.max() < tcfg["target_power"]:
        n80 = np.nan
    elif pw.min() >= tcfg["target_power"]:
        n80 = float(g.min())
    else:
        j = np.argmax(pw >= tcfg["target_power"])
        x0, x1, y0, y1 = g[j - 1], g[j], pw[j - 1], pw[j]
        n80 = float(x0 + (tcfg["target_power"] - y0) * (x1 - x0) / (y1 - y0))
    return pc, n80


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    cfg = yaml.safe_load((HERE / "config.yaml").read_text())
    outdir = HERE / cfg["outputs_dir"]
    outdir.mkdir(exist_ok=True)

    print("=" * 74)
    print("H-008 SIMULATION FRAMEWORK -- 100% SIMULATED DATA (no PRO-ACT records)")
    print("=" * 74)

    rng = np.random.default_rng(cfg["seed"])
    df, calib = simulate_cohort(cfg, rng)
    print(f"[simulated pool] n={len(df)}  calibration={calib}")

    masks, kept, diag = build_arm_scores(df, cfg)
    print(f"[eligible per arm] {kept}")
    if diag:
        print(f"[diagnostics] nested-model F-test p={diag['nfl_slope_beyond_baseline_F_pvalue']:.3g}"
              "  (criterion 3: early NfL slope must add beyond baseline)")

    dur = cfg["trial"]["duration_months"]
    deltas = -df.lambda_true_pts_per_month.values * dur  # untreated 12-mo delta

    results, curves = [], []
    for eff in cfg["trial"]["treatment_effect_pts_per_month"]:
        n80s = {}
        for arm in cfg["arms"]:
            pc, n80 = power_and_n80(deltas, masks[arm], cfg, eff,
                                    np.random.default_rng(cfg["seed"] + int(round(eff * 1000))))
            pc["arm"], pc["effect_pts_per_month"] = arm, eff
            curves.append(pc)
            n80s[arm] = n80
        ref = n80s["none"]
        for arm, n80 in n80s.items():
            red = 1 - n80 / ref if np.isfinite(n80) else np.nan
            screen_ratio = len(df) / kept[arm] if kept[arm] else np.nan
            results.append({"effect_pts_per_month": eff, "arm": arm,
                            "n_per_arm_for_80pct_power": round(n80, 1) if np.isfinite(n80) else None,
                            "sample_size_reduction_vs_none": round(red, 3) if np.isfinite(red) else None,
                            "patients_screened_per_enrolled": round(screen_ratio, 2)})
        print(f"\n[effect {eff} pts/month] n@80% power per arm:")
        for r_ in results[-len(n80s):]:
            print(f"   {r_['arm']:<16} n={r_['n_per_arm_for_80pct_power']}  "
                  f"reduction={r_['sample_size_reduction_vs_none']}")

    res = pd.DataFrame(results)
    cur = pd.concat(curves, ignore_index=True)

    # ---- Rule-5 controls --------------------------------------------------
    failures = []
    for eff in cfg["trial"]["treatment_effect_pts_per_month"]:
        sub = res[res.effect_pts_per_month == eff].set_index("arm")
        orc = sub.loc["oracle", "sample_size_reduction_vs_none"] or 0
        noi = sub.loc["noise", "sample_size_reduction_vs_none"] or 0
        if orc < 0.25:
            failures.append(f"POSITIVE CONTROL FAILED (eff={eff}): oracle reduction {orc:.0%} < 25%")
        if noi > 0.05 or noi < -0.10:
            failures.append(f"NEGATIVE CONTROL FAILED (eff={eff}): noise reduction {noi:.0%}")
    if failures:
        for f_ in failures:
            print("!!", f_)
        print("Per AGENTS.md rule 5: outputs are NOT trustworthy. Fix before use.")
        sys.exit(1)

    header = (
        "# H-008 NfL-enrichment simulation summary\n\n"
        "> **DATA STATUS: 100% SIMULATED.** No PRO-ACT or other real patient\n"
        "> records were used anywhere in this run. Every number below is a\n"
        "> property of the generative model in `config.yaml`, calibrated to\n"
        "> published summaries (PMIDs 25362243, 29598923, 35585374, 30014505,\n"
        "> 31432691, 34690913, 38674431, 31280619). Real-data validation is\n"
        "> **PENDING**: see `data/proact/APPLICATION_DRAFT.md`.\n\n"
        "## Effect-size assumption\n"
        "Treatment slows ALSFRS-R decline by an absolute "
        f"{cfg['trial']['treatment_effect_pts_per_month']} pts/month scenarios;\n"
        "the headline claim uses 0.30 pts/month (H-008 realistic range 0.3-0.5).\n\n"
        f"- Pool calibration: {calib}\n"
        f"- Eligible subjects per arm (fraction {cfg['trial']['enrich_fraction']}): {kept}\n"
        f"- Nested-model test, NfL slope beyond baseline: p="
        f"{diag.get('nfl_slope_beyond_baseline_F_pvalue', float('nan')):.3g}\n\n"
        "| effect (pts/mo) | arm | n/arm @80% power | reduction vs none | screened/enrolled |\n"
        "|---|---|---|---|---|\n"
    )
    table = "".join(
        f"| {r.effect_pts_per_month} | {r.arm} | {r.n_per_arm_for_80pct_power} "
        f"| {r.sample_size_reduction_vs_none} | {r.patients_screened_per_enrolled} |\n"
        for r in res.itertuples()
    )
    pending = (
        "\n## Pending real data\n"
        "- Re-run with PRO-ACT placebo-arm trajectories once access is granted\n"
        "  (application draft: `data/proact/APPLICATION_DRAFT.md`).\n"
        "- Replace simulated NfL blocks with observed PRO-ACT/NfL values where\n"
        "  available, else keep progression-rate surrogates calibrated to the\n"
        "  NfL literature (per H-008 falsification criterion 1).\n"
        "- ENCALS comparator parameters requested in APPLICATION_DRAFT.md.\n"
    )
    (outdir / "simulation_summary.md").write_text(header + table + pending)
    res.to_csv(outdir / "sample_sizes.csv", index=False)
    cur.to_csv(outdir / "power_curves.csv", index=False)
    (outdir / "run_meta.json").write_text(json.dumps(
        {"seed": cfg["seed"], "data_status": "SIMULATED_ONLY",
         "pool_calibration": calib, "nested_model_p": diag}, indent=2))
    print(f"\n[ok] wrote {outdir}/simulation_summary.md, sample_sizes.csv, "
          f"power_curves.csv, run_meta.json")
    print("[status] end-to-end on SIMULATED data; real-data validation pending")


if __name__ == "__main__":
    main()
