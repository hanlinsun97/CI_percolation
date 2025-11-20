from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLUMN_MAP: Sequence[Tuple[int, str]] = [
    (1, "run"),
    (2, "p_plus"),
    (3, "p_plus_jump"),
    (4, "P_p_plus"),
    (5, "s_p_plus"),
    (6, "nc_p_plus"),
    (7, "nl_p_plus"),
    (8, "q_plus"),
    (9, "P_q_plus"),
    (10, "s_q_plus"),
    (11, "nc_q_plus"),
    (12, "nl_q_plus"),
    (13, "q_zero"),
    (14, "P_q_zero"),
    (15, "s_q_zero"),
    (16, "nc_q_zero"),
    (17, "nl_q_zero"),
    (18, "q_minus"),
    (19, "P_q_minus"),
    (20, "s_q_minus"),
    (21, "nc_q_minus"),
    (22, "nl_q_minus"),
    (23, "p2"),
    (24, "p1"),
]

COLUMN_LOOKUP: Dict[int, str] = {idx: name for idx, name in COLUMN_MAP}

N_PATTERN = re.compile(r"_N(\d+)")
L_PATTERN = re.compile(r"_L(\d+)")
MG_PATTERN = re.compile(r"_Mg(\d+)")


def read_ebe_file(path: Path) -> Tuple[pd.DataFrame, Dict[str, str]]:
    columns = [name for _, name in COLUMN_MAP]
    rows: List[List[float]] = []
    meta: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped.startswith("# N=") or stripped.startswith("# M_graphs="):
                    parts = stripped[1:].split()
                    for part in parts:
                        if "=" in part:
                            key, value = part.split("=", 1)
                            meta[key.strip()] = value.strip()
                elif stripped.startswith("# STATISTICS"):
                    break
                continue
            values = [float(x) for x in stripped.split()]
            if len(values) == len(columns) - 1:
                # Backward compatibility: files without p_plus_jump.
                # Insert placeholder (NaN) right after p_plus.
                values = values[:2] + [float("nan")] + values[2:]
            if len(values) != len(columns):
                raise ValueError(f"Unexpected column count in {path.name}")
            rows.append(values)
    df = pd.DataFrame(rows, columns=columns)
    return df, meta


def welford_std_err(series: pd.Series) -> Tuple[float, float, float, float]:
    values = series.to_numpy(dtype=float)
    values = values[~np.isnan(values)]
    n = len(values)
    mean = float(values.mean())
    std = float(values.std(ddof=1)) if n > 1 else 0.0
    err_mean = std / math.sqrt(n) if n > 1 else 0.0
    if n > 3 and std > 0:
        kurt_excess = float(pd.Series(values).kurtosis())
        err_std = std * math.sqrt(max(kurt_excess + 2.0, 0.0) / (4.0 * n))
    else:
        err_std = 0.0
    return mean, err_mean, std, err_std


def compute_column_stats(df: pd.DataFrame) -> pd.DataFrame:
    records = []
    for idx, name in COLUMN_MAP:
        mean, err_mean, std, err_std = welford_std_err(df[name])
        records.append(
            {
                "row": idx,
                "label": name,
                "mean": mean,
                "err_mean": err_mean,
                "std": std,
                "err_std": err_std,
            }
        )
    return pd.DataFrame(records)


def write_stats_file(
    path: Path, stats: pd.DataFrame, input_file: Path, runs: int
) -> None:
    header = [
        "# Statistical Analysis Results",
        f"# Input file: {input_file.name}",
        f"# Number of runs: {runs}",
        f"# Number of columns: {len(COLUMN_MAP)}",
        "#",
        "# Standard error of std dev uses exact formula: SE(σ) = σ * sqrt[(κ + 2) / (4M)]",
        "#",
        "# Column format:",
        "# (1) Column_ID (2) Mean (3) Error_Mean (4) Std_Dev (5) Error_Std_Dev",
        "#",
    ]
    body = [
        f"{int(row):3d} {mean:15.8f} {err_m:15.8f} {std:15.8f} {err_s:15.8f}"
        for row, mean, err_m, std, err_s in stats[
            ["row", "mean", "err_mean", "std", "err_std"]
        ].to_numpy()
    ]
    footer = [
        "#",
        "# Summary:",
        f"# Successfully processed {runs} runs with {len(COLUMN_MAP)} columns each",
    ]
    content = "\n".join(header + body + footer) + "\n"
    path.write_text(content)


def process_all_ebe(
    data_dir: Path,
    output_dir: Optional[Path],
    pattern: str = "EBE_*.dat",
    write_files: bool = True,
    ci_radius_filter: Optional[int] = None,
) -> pd.DataFrame:
    paths = sorted(data_dir.glob(pattern))
    records = []
    columns = ["file", "runs", "N", "ci_radius", "M_graphs", "stats"]
    for path in paths:
        df, meta = read_ebe_file(path)
        stats = compute_column_stats(df)
        if write_files:
            if output_dir is None:
                raise ValueError("output_dir must be provided when write_files=True")
            stats_path = output_dir / f"stats_{path.name}"
            write_stats_file(stats_path, stats, path, len(df))
        ci_radius = extract_ci_radius(path.name)
        m_graphs = extract_m_graphs(path.name)
        if m_graphs is None:
            mg_meta = meta.get("M_graphs")
            m_graphs = int(mg_meta) if mg_meta else None
        if ci_radius_filter is not None and ci_radius != ci_radius_filter:
            continue
        N_value = int(meta.get("N", len(df)))
        records.append(
            {
                "file": path.name,
                "runs": len(df),
                "N": N_value,
                "ci_radius": ci_radius,
                "M_graphs": m_graphs,
                "stats": stats,
            }
        )
    if not records:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame.from_records(records, columns=columns)


def read_stats_file(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 5:
                continue
            row_id = int(parts[0])
            rows.append(
                {
                    "row": row_id,
                    "label": COLUMN_LOOKUP.get(row_id, f"col_{row_id}"),
                    "mean": float(parts[1]),
                    "err_mean": float(parts[2]),
                    "std": float(parts[3]),
                    "err_std": float(parts[4]),
                }
            )
    return pd.DataFrame(rows)


def extract_system_size(name: str) -> int:
    match = N_PATTERN.search(name)
    if not match:
        raise ValueError(f"Cannot find N in {name}")
    return int(match.group(1))


def extract_ci_radius(name: str) -> Optional[int]:
    match = L_PATTERN.search(name)
    return int(match.group(1)) if match else None


def extract_m_graphs(name: str) -> Optional[int]:
    match = MG_PATTERN.search(name)
    return int(match.group(1)) if match else None


def extract_fss(
    stats_paths: Iterable[Path], row: int, value_kind: str = "mean"
) -> pd.DataFrame:
    assert value_kind in {"mean", "std"}
    records = []
    for path in sorted(stats_paths):
        stats = read_stats_file(path)
        row_data = stats.loc[stats["row"] == row]
        if row_data.empty:
            continue
        N = extract_system_size(path.name)
        if value_kind == "mean":
            value = float(row_data["mean"])
            error = float(row_data["err_mean"])
        else:
            value = float(row_data["std"])
            error = float(row_data["err_std"])
        records.append(
            {"N": N, "value": value, "error": error, "source": path.name}
        )
    return pd.DataFrame(records).sort_values("N").reset_index(drop=True)


def extract_p1_minus_p2(
    stats_paths: Iterable[Path], value_kind: str = "mean"
) -> pd.DataFrame:
    df_p1 = extract_fss(stats_paths, row=24, value_kind=value_kind)
    df_p2 = extract_fss(stats_paths, row=23, value_kind=value_kind)
    merged = pd.merge(df_p1, df_p2, on="N", suffixes=("_p1", "_p2"))
    value = merged["value_p1"] - merged["value_p2"]
    error = np.sqrt(merged["error_p1"] ** 2 + merged["error_p2"] ** 2)
    result = merged[["N"]].copy()
    result["value"] = value
    result["error"] = error
    return result


def weighted_log_fit(fss: pd.DataFrame) -> Dict[str, float]:
    values = fss["value"].to_numpy(dtype=float)
    errors = fss["error"].to_numpy(dtype=float)
    nonzero = values[np.nonzero(values)]
    base_sign = float(np.sign(nonzero[0])) if nonzero.size else 1.0
    used_abs = False
    if np.any(values <= 0):
        unique_signs = np.unique(np.sign(nonzero)) if nonzero.size else np.array([1.0])
        if len(unique_signs) > 1:
            raise ValueError("FSS data crosses zero; cannot perform log-log fit.")
        used_abs = True
        working_values = np.abs(values)
    else:
        working_values = values
    if np.any(working_values <= 0):
        raise ValueError("FSS data contains zeros; cannot take logarithm.")
    logN = np.log(fss["N"].to_numpy(dtype=float))
    logY = np.log(working_values)
    sigma = errors / working_values
    sigma[~np.isfinite(sigma) | (sigma <= 0)] = 1e-9
    w = 1.0 / (sigma**2)
    sw = w.sum()
    sx = (w * logN).sum()
    sy = (w * logY).sum()
    sxx = (w * logN**2).sum()
    sxy = (w * logN * logY).sum()
    delta = sw * sxx - sx * sx
    if abs(delta) < 1e-12:
        raise ValueError("Singular data matrix; check FSS input.")
    a = (sxx * sy - sx * sxy) / delta
    b = (sw * sxy - sx * sy) / delta
    siga = math.sqrt(sxx / delta)
    sigb = math.sqrt(sw / delta)
    A_mag = math.exp(a)
    A_err = A_mag * siga
    alpha = b
    alpha_err = sigb
    residuals = (logY - a - b * logN) / sigma
    chi2 = float(np.sum(residuals**2))
    dof = len(fss) - 2
    A_signed = (-1.0 if used_abs and base_sign < 0 else 1.0) * A_mag
    return {
        "A": A_signed,
        "A_err": A_err,
        "alpha": alpha,
        "alpha_err": alpha_err,
        "chi2_red": chi2 / dof if dof > 0 else np.nan,
        "fit_used_abs": used_abs,
        "fit_value_sign": base_sign if used_abs else 1.0,
    }


def _windowed_weighted_slope(logX: np.ndarray, logY: np.ndarray, sigma: np.ndarray) -> Tuple[float, float]:
    sigma = sigma.copy()
    sigma[~np.isfinite(sigma) | (sigma <= 0)] = 1e-9
    w = 1.0 / (sigma**2)
    sw = w.sum()
    sx = (w * logX).sum()
    sy = (w * logY).sum()
    sxx = (w * logX**2).sum()
    sxy = (w * logX * logY).sum()
    delta = sw * sxx - sx * sx
    if abs(delta) < 1e-16:
        return math.nan, math.nan
    slope = (sw * sxy - sx * sy) / delta
    slope_err = math.sqrt(sw / delta)
    return slope, slope_err


def compute_effective_exponents(fss: pd.DataFrame, window: int = 2) -> pd.DataFrame:
    """Sliding-window effective exponents via weighted log-log fits."""
    if window < 2:
        raise ValueError("window must be at least 2")
    if len(fss) < window:
        return pd.DataFrame(columns=["N_left", "N_right", "N_geom", "exponent", "exponent_err"])
    values = fss["value"].to_numpy(dtype=float)
    nonzero = values[np.nonzero(values)]
    if np.any(values <= 0):
        unique_signs = np.unique(np.sign(nonzero)) if nonzero.size else np.array([1.0])
        if len(unique_signs) > 1:
            raise ValueError(
                "Effective exponent series crosses zero; cannot log-transform."
            )
        working_values = np.abs(values)
    else:
        working_values = values
    if np.any(working_values <= 0):
        raise ValueError("Effective exponent series contains zeros; cannot log-transform.")
    logN = np.log(fss["N"].to_numpy(dtype=float))
    logY = np.log(working_values)
    sigma = (fss["error"] / working_values).to_numpy(dtype=float)

    rows = []
    for start in range(0, len(fss) - window + 1):
        end = start + window
        local_logN = logN[start:end]
        local_logY = logY[start:end]
        local_sigma = sigma[start:end]
        slope, slope_err = _windowed_weighted_slope(local_logN, local_logY, local_sigma)
        if not np.isfinite(slope):
            continue
        rows.append(
            {
                "N_left": float(fss["N"].iloc[start]),
                "N_right": float(fss["N"].iloc[end - 1]),
                "N_geom": float(np.exp(local_logN.mean())),
                "exponent": slope,
                "exponent_err": slope_err,
            }
        )
    return pd.DataFrame(rows)


def plot_fss_with_fit(fss: pd.DataFrame, fit: Dict[str, float], title: str = "") -> None:
    values = fss["value"].to_numpy(dtype=float)
    errors = fss["error"].to_numpy(dtype=float)
    using_abs = bool(np.any(values <= 0))
    display_values = np.abs(values) if using_abs else values
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        fss["N"],
        display_values,
        yerr=errors,
        fmt="o",
        label="data" if not using_abs else "data (|value|)",
    )
    xs = np.linspace(fss["N"].min(), fss["N"].max(), 200)
    amp = abs(fit["A"]) if using_abs else fit["A"]
    ys = amp * xs ** fit["alpha"]
    ax.plot(xs, ys, label=f"fit: α={fit['alpha']:.4f}±{fit['alpha_err']:.4f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("N")
    ylabel = "Observable" if not using_abs else "|Observable|"
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title + (" (|values| shown)" if using_abs else ""))
    elif using_abs:
        ax.set_title("|values| shown")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.3)
    plt.show()


def plot_effective_exponent(eff: pd.DataFrame, title: str = "") -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.errorbar(eff["N_geom"], eff["exponent"], yerr=eff["exponent_err"], fmt="o-")
    ax.set_xscale("log")
    ax.set_xlabel("N (geometric mean)")
    ax.set_ylabel("Effective exponent")
    if title:
        ax.set_title(title)
    ax.grid(True, which="both", ls="--", alpha=0.3)
    plt.show()


def read_nr_file(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 12:
                continue
            nodes_remaining = int(parts[0])
            frac_removed = float(parts[1])
            gc_mean = float(parts[2])
            rows.append(
                {
                    "nodes_remaining": nodes_remaining,
                    "p": frac_removed,
                    "gc": gc_mean,
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df["P_inf"] = df["gc"]
    return df


def derive_output_prefix(prefix: str) -> str:
    if "_EBE_" in prefix:
        return prefix.split("_EBE_")[0]
    if prefix.endswith("_"):
        return prefix[:-1]
    return prefix


def plot_giant_component_vs_p(
    data_dir: Path, file_prefix: str, ci_radius_filter: Optional[int] = None
) -> None:
    output_prefix = derive_output_prefix(file_prefix)
    nr_pattern = f"{output_prefix}_NR_ER_*.dat"
    nr_files = sorted(data_dir.glob(nr_pattern))
    if not nr_files:
        print("No NR files found matching pattern", nr_pattern)
        return

    def key_fn(p: Path) -> Tuple[int, int]:
        size = extract_system_size(p.name)
        radius = extract_ci_radius(p.name)
        if ci_radius_filter is None:
            return (size, radius if radius is not None else -1)
        if radius != ci_radius_filter:
            return (-1, -1)
        return (size, radius if radius is not None else -1)

    candidate = max(nr_files, key=key_fn)
    radius = extract_ci_radius(candidate.name)
    if ci_radius_filter is not None and radius != ci_radius_filter:
        print(f"No NR files found with L={ci_radius_filter}. Showing largest N overall.")
    df = read_nr_file(candidate)
    if df.empty:
        print("Selected file has no data:", candidate.name)
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(df["p"], df["P_inf"], marker="o")
    ax.set_xlabel("Fraction of removed nodes p")
    ax.set_ylabel("P^∞ (giant component fraction)")
    ax.set_title(f"Largest N={extract_system_size(candidate.name)} from {candidate.name}")
    ax.grid(True, ls="--", alpha=0.3)
    plt.show()


fit_power_law = weighted_log_fit
effective_exponents = compute_effective_exponents
