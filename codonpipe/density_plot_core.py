"""codonpipe.density_plot_core

Core plotting routines for cluster overlays in 2D embedding space.
"""

# codonpipe/density_plot_core.py
import os
import math
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

from .colormaps import get_cmap_any
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from scipy.ndimage import gaussian_filter


_LT_INLINE = re.compile(r"\[locus_tag=([^\]]+)\]")


def _canon_tag(s):
    if s is None:
        return ""
    t = str(s).strip()
    if not t or t.lower() == "nan":
        return ""
    m = _LT_INLINE.search(t)
    if m:
        return m.group(1).strip()
    if t.startswith(">"):
        t = t[1:].strip()
    if " " in t:
        t = t.split()[0].strip()
    return t


def _pick_excel():
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.update()
        path = filedialog.askopenfilename(
            title="Select the Excel workbook",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        root.destroy()
        return path
    except Exception:
        return input("Excel file path: ").strip('"').strip()


def _auto_axis_limits(x, y, pad_frac=0.05):
    xmin, xmax = np.nanmin(x), np.nanmax(x)
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    xpad = (xmax - xmin) * pad_frac if xmax > xmin else 1.0
    ypad = (ymax - ymin) * pad_frac if ymax > ymin else 1.0
    return (xmin - xpad, xmax + xpad), (ymin - ypad, ymax + ypad)


def _nice_step(range_val, target_intervals=5):
    if range_val <= 0:
        return 1.0
    raw = range_val / max(target_intervals, 1)
    power = 10 ** math.floor(math.log10(raw))
    norm = raw / power
    if norm < 1.5:
        nice = 1.0
    elif norm < 3.0:
        nice = 2.0
    elif norm < 7.0:
        nice = 5.0
    else:
        nice = 10.0
    return nice * power


def _draw_grid(ax, xlim, ylim, major_step=None, major_alpha=0.55, major_lw=1.0, major_color="black"):
    xmin, xmax = xlim
    ymin, ymax = ylim
    x_range = xmax - xmin
    y_range = ymax - ymin

    if major_step is None or major_step <= 0:
        x_step = _nice_step(x_range, target_intervals=5)
        y_step = _nice_step(y_range, target_intervals=5)
    else:
        x_step = major_step
        y_step = major_step

    def _ticks(vmin, vmax, step):
        if step <= 0:
            return []
        start = math.ceil(vmin / step) * step
        end = math.floor(vmax / step) * step
        if end < start:
            return []
        n = int(round((end - start) / step)) + 1
        return [start + i * step for i in range(n)]

    for x in _ticks(xmin, xmax, x_step):
        ax.axvline(x, ls=(0, (4, 4)), color=major_color, lw=major_lw, alpha=major_alpha, zorder=5)
    for y in _ticks(ymin, ymax, y_step):
        ax.axhline(y, ls=(0, (4, 4)), color=major_color, lw=major_lw, alpha=major_alpha, zorder=5)


def _wrap_cluster_title(name, n_value, enabled=True, width=35, max_name_lines=2):
    base = str(name)
    if not enabled:
        return f"{base}\n n = {n_value}"

    wrapper = textwrap.TextWrapper(width=width, break_long_words=False, break_on_hyphens=False)
    lines = wrapper.wrap(base)

    if max_name_lines is not None and max_name_lines > 0 and len(lines) > max_name_lines:
        head = lines[:max_name_lines - 1]
        tail = " ".join(lines[max_name_lines - 1:])
        lines = head + [tail]

    wrapped = "\n".join(lines)  # avoid backslash inside f-string expression
    return f"{wrapped}\n n = {n_value}"


def _load_custom_colormaps(xlsx_path):
    maps = {}
    if not xlsx_path or (not os.path.isfile(xlsx_path)):
        if xlsx_path:
            print(f"[WARN] Custom cmap Excel not found: {xlsx_path}")
        return maps
    try:
        xls = pd.ExcelFile(xlsx_path)
        for sheet in xls.sheet_names:
            df = xls.parse(sheet, header=None)
            if df.shape[1] < 4:
                continue

            # Robust: ignore non-numeric rows (headers like 'R','G','B')
            rgb_raw = df.iloc[:, 1:4].apply(pd.to_numeric, errors="coerce")
            rgb = rgb_raw.values
            rgb = rgb[~np.isnan(rgb).any(axis=1)]
            if rgb.shape[0] < 2:
                continue

            if np.nanmax(rgb) > 1.0:
                rgb = np.clip(rgb / 255.0, 0.0, 1.0)
            else:
                rgb = np.clip(rgb, 0.0, 1.0)

            maps[sheet] = LinearSegmentedColormap.from_list(sheet, [tuple(c) for c in rgb], N=256)

        if maps:
            print(f"[INFO] Loaded custom colormaps: {', '.join(maps.keys())}")
    except Exception as e:
        print(f"[WARN] Failed to parse custom colormaps: {e}")
    return maps


def _choose_cmap(name, custom_maps):
    # Supports built-in 'parula' / 'parula_r' via codonpipe.colormaps
    return get_cmap_any(name, custom_maps=custom_maps, fallback="plasma")


def _compute_density_colors(Y, nbins=150, sigma=4.0, use_log=True, min_rel=0.0, max_rel=1.0, cmap=None):
    if cmap is None:
        cmap = plt.get_cmap("plasma")
    if Y is None or Y.shape[0] == 0:
        return None

    x = Y[:, 0]
    y = Y[:, 1]
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    x_pad = 0.05 * x_range if x_range > 0 else 1.0
    y_pad = 0.05 * y_range if y_range > 0 else 1.0
    x_min = x.min() - x_pad
    x_max = x.max() + x_pad
    y_min = y.min() - y_pad
    y_max = y.max() + y_pad

    xedges = np.linspace(x_min, x_max, nbins + 1)
    yedges = np.linspace(y_min, y_max, nbins + 1)
    H, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])

    Hs = gaussian_filter(H, sigma=sigma, mode="nearest")
    if use_log:
        Hs = np.log1p(Hs)

    if not np.isfinite(Hs).any() or Hs.max() <= 0:
        return cmap(np.zeros_like(x))

    ix = np.searchsorted(xedges, x, side="right") - 1
    iy = np.searchsorted(yedges, y, side="right") - 1
    ix = np.clip(ix, 0, nbins - 1)
    iy = np.clip(iy, 0, nbins - 1)
    dens = Hs[ix, iy]

    dmax = dens.max()
    min_rel = max(0.0, min(1.0, float(min_rel)))
    max_rel = max(min_rel, min(1.0, float(max_rel)))

    vmin = min_rel * dmax
    vmax = max(vmin + 1e-12, max_rel * dmax)

    dens_clipped = np.clip(dens, vmin, vmax)
    norm = (dens_clipped - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)
    return cmap(norm)


# ---------------- Enrichment math (parameterized by cfg) ----------------
def _edges_from_limits(xlim, ylim, nbins: int):
    return np.linspace(xlim[0], xlim[1], nbins + 1), np.linspace(ylim[0], ylim[1], nbins + 1)


def _density_field(points_xy, xedges, yedges, sigma):
    if points_xy is None or points_xy.size == 0:
        return np.zeros((len(xedges) - 1, len(yedges) - 1), dtype=float)
    H, _, _ = np.histogram2d(points_xy[:, 0], points_xy[:, 1], bins=(xedges, yedges))
    return gaussian_filter(H, sigma=sigma, mode="nearest")


def _sample_field(field, xedges, yedges, points_xy):
    if points_xy is None or points_xy.size == 0:
        return np.array([], dtype=float)
    x = points_xy[:, 0]
    y = points_xy[:, 1]
    ix = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, field.shape[0] - 1)
    iy = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, field.shape[1] - 1)
    return field[ix, iy].astype(float)


def _rgb_luma(rgb3):
    r, g, b = float(rgb3[0]), float(rgb3[1]), float(rgb3[2])
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _looks_diverging(cmap):
    """
    Heuristic: diverging maps typically have a bright/neutral middle (near-white)
    and two darker, different-colored ends.
    """
    c0 = np.array(cmap(0.0)[:3], float)
    c1 = np.array(cmap(1.0)[:3], float)
    cm = np.array(cmap(0.5)[:3], float)

    l0, l1, lm = _rgb_luma(c0), _rgb_luma(c1), _rgb_luma(cm)
    ends_diff = float(np.linalg.norm(c0 - c1))

    # middle noticeably brighter than both ends + ends clearly different
    return (lm > l0 + 0.15) and (lm > l1 + 0.15) and (ends_diff > 0.40)


def _unilateral_cmap(base_cmap, side="high", mode="auto"):
    """
    For unilateral enrichment:
      - If cmap is diverging: crop to one half (center->end) and renormalize (uses full 0..1).
      - If cmap is sequential: DO NOT crop (cropping makes the bar look "partial").
    mode:
      - "auto" (default): crop only if it looks diverging
      - "crop": always crop (even sequential; generally NOT recommended)
      - "full": never crop
    side:
      - "high": keep 0.5..1.0
      - "low" : keep 0.0..0.5
    """
    mode = (mode or "auto").lower().strip()
    side = (side or "high").lower().strip()

    if mode == "full":
        return base_cmap

    do_crop = (mode == "crop") or (mode == "auto" and _looks_diverging(base_cmap))
    if not do_crop:
        return base_cmap

    if side == "low":
        xs = np.linspace(0.0, 0.5, 256)
        tag = "lowhalf"
    else:
        xs = np.linspace(0.5, 1.0, 256)
        tag = "highhalf"

    cols = base_cmap(xs)
    return LinearSegmentedColormap.from_list(f"{base_cmap.name}_{tag}", cols, N=256)


def _auto_enrichment_vmax(vals, scale, percentile, style):
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 1.0

    style = (style or "").lower().strip()

    if scale == "log2":
        if style == "diverging":
            vmax = float(np.nanpercentile(np.abs(vals), percentile))
            return max(vmax, 1e-6)
        # unilateral: only care about >0
        vpos = vals[vals > 0]
        if vpos.size == 0:
            return 1e-3
        vmax = float(np.nanpercentile(vpos, percentile))
        return max(vmax, 1e-6)

    # ratio scale: neutral is 1.0
    if style == "diverging":
        v = np.abs(np.log(vals[vals > 0]))  # symmetric in log-space
        if v.size == 0:
            return 1.0 + 1e-6
        vmax_log = float(np.nanpercentile(v, percentile))
        vmax = float(np.exp(vmax_log))
        return max(vmax, 1.0 + 1e-6)

    # unilateral: only values > 1
    vpos = vals[vals > 1.0]
    if vpos.size == 0:
        return 1.0 + 1e-6
    vmax = float(np.nanpercentile(vpos, percentile))
    return max(vmax, 1.0 + 1e-6)


def _enrichment_norm_and_cmap(
    vmax,
    base_cmap,
    scale,
    style,
    symmetric,
    unilateral_side="high",
    unilateral_cmap_mode="auto",
):
    """
    style must be either:
      - "diverging"
      - "unilateral_diverging"  (anything other than "diverging" is treated as unilateral)
    """
    style = (style or "").lower().strip()

    if scale == "log2":
        center = 0.0
        if style == "diverging":
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=center, vmax=vmax)
            return norm, base_cmap, center

        # unilateral: only show positive log2 enrichment
        cmap_u = _unilateral_cmap(base_cmap, side=unilateral_side, mode=unilateral_cmap_mode)
        norm = Normalize(vmin=0.0, vmax=vmax)  # 0..vmax uses full cmap_u range
        return norm, cmap_u, center

    # ratio scale
    center = 1.0
    if style == "diverging":
        # symmetric around 1 in log-space -> vmin = 1/vmax
        vmin = 1.0 / vmax if vmax > 0 else 1.0
        norm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
        return norm, base_cmap, center

    # unilateral: show only enrichment > 1
    cmap_u = _unilateral_cmap(base_cmap, side=unilateral_side, mode=unilateral_cmap_mode)
    norm = Normalize(vmin=1.0, vmax=vmax)  # 1..vmax uses full cmap_u range
    return norm, cmap_u, center


def _compute_enrichment_values(cluster_xy, all_xy, xedges, yedges, sigma, eps, scale, field_all_cache=None):
    if cluster_xy is None or cluster_xy.size == 0:
        return np.array([], dtype=float)

    n_all = max(int(all_xy.shape[0]), 1)
    n_clu = max(int(cluster_xy.shape[0]), 1)

    field_all = field_all_cache if field_all_cache is not None else _density_field(all_xy, xedges, yedges, sigma=sigma)
    field_clu = _density_field(cluster_xy, xedges, yedges, sigma=sigma)

    dens_all = _sample_field(field_all, xedges, yedges, cluster_xy) / n_all
    dens_clu = _sample_field(field_clu, xedges, yedges, cluster_xy) / n_clu

    ratio = (dens_clu + eps) / (dens_all + eps)
    if scale == "log2":
        ratio = np.log2(ratio)
    return ratio


def _map_enrichment_to_colors(vals, norm, cmap, center_value, force_neutral, neutral_rgba):
    """
    IMPORTANT: norm+cmap MUST match the colorbar.
    If force_neutral is True, values <= center_value are shown as neutral.
    """
    vals = np.asarray(vals, float)
    cols = cmap(norm(vals))

    if force_neutral:
        mask = vals <= center_value
        cols[mask] = np.array(neutral_rgba, dtype=float)
    return cols


def _get_cluster_preset(n, presets):
    for preset in presets:
        min_n = preset.get("MIN_N", 0)
        max_n = preset.get("MAX_N", None)
        if n < min_n:
            continue
        if max_n is None or n < max_n:
            return preset
    return presets[-1] if presets else None


# =========================================================
# Entry point called from the short plotting script
# =========================================================
def run_density_dotplot(cfg: dict):
    # Pull settings (same names as the calling script)
    EXCEL_PATH = cfg.get("EXCEL_PATH", "")
    UMAP_SHEET = cfg.get("UMAP_SHEET", "tsne coordinates")
    DATA_SHEET = cfg.get("DATA_SHEET", "Binary")
    USE_FILE_DIALOG = bool(cfg.get("USE_FILE_DIALOG", True))

    INCLUDE_COLUMNS = cfg.get("INCLUDE_COLUMNS", None)
    HIGHLIGHT_COLUMNS = cfg.get("HIGHLIGHT_COLUMNS", None)

    PLOT_MODE = cfg.get("PLOT_MODE", "grid")
    MAX_NROWS = int(cfg.get("MAX_NROWS", 4))
    PANEL_W_IN = float(cfg.get("PANEL_W_IN", 5.0))
    PANEL_H_IN = float(cfg.get("PANEL_H_IN", 5.0))
    X_LABEL = cfg.get("X_LABEL", "t-SNE 1")
    Y_LABEL = cfg.get("Y_LABEL", "t-SNE 2")

    SHOW_BACKGROUND_SMALL = bool(cfg.get("SHOW_BACKGROUND_SMALL", True))
    BACKGROUND_COLOR = cfg.get("BACKGROUND_COLOR", "lightgrey")
    BACKGROUND_ALPHA = float(cfg.get("BACKGROUND_ALPHA", 0.3))
    BACKGROUND_SIZE = float(cfg.get("BACKGROUND_SIZE", 8))
    BACKGROUND_LW = float(cfg.get("BACKGROUND_LW", 0.0))

    HIGHLIGHT_MARKER = cfg.get("HIGHLIGHT_MARKER", "o")
    HIGHLIGHT_EDGE_COLOR = cfg.get("HIGHLIGHT_EDGE_COLOR", "black")
    HIGHLIGHT_EDGE_WIDTH = float(cfg.get("HIGHLIGHT_EDGE_WIDTH", 0.6))
    HIGHLIGHT_ALPHA = float(cfg.get("HIGHLIGHT_ALPHA", 0.95))

    ALL_SCATTER_POINT_SIZE = float(cfg.get("ALL_SCATTER_POINT_SIZE", 3.5))
    ALL_SCATTER_ALPHA = float(cfg.get("ALL_SCATTER_ALPHA", 0.7))
    ALL_SCATTER_EDGE_WIDTH = float(cfg.get("ALL_SCATTER_EDGE_WIDTH", 0.15))
    ALL_SCATTER_DENSITY_NBINS = int(cfg.get("ALL_SCATTER_DENSITY_NBINS", 150))
    ALL_SCATTER_DENSITY_SIGMA = float(cfg.get("ALL_SCATTER_DENSITY_SIGMA", 4.0))
    ALL_SCATTER_DENSITY_USE_LOG = bool(cfg.get("ALL_SCATTER_DENSITY_USE_LOG", True))
    ALL_SCATTER_DENSITY_MIN_REL = float(cfg.get("ALL_SCATTER_DENSITY_MIN_REL", 0.00))
    ALL_SCATTER_DENSITY_MAX_REL = float(cfg.get("ALL_SCATTER_DENSITY_MAX_REL", 1.00))

    CLUSTER_SCATTER_PRESETS = cfg.get("CLUSTER_SCATTER_PRESETS", [])

    DRAW_GRID = bool(cfg.get("DRAW_GRID", True))
    GRID_MAJOR_STEP = cfg.get("GRID_MAJOR_STEP", None)
    GRID_MAJOR_ALPHA = float(cfg.get("GRID_MAJOR_ALPHA", 0.55))
    GRID_MAJOR_LW = float(cfg.get("GRID_MAJOR_LW", 1.0))
    GRID_MAJOR_COLOR = cfg.get("GRID_MAJOR_COLOR", "black")

    CAPTION_SIZE = float(cfg.get("CAPTION_SIZE", 12))
    CAPTION_WEIGHT = cfg.get("CAPTION_WEIGHT", "bold")
    CAPTION_PAD = float(cfg.get("CAPTION_PAD", 8))
    CAPTION_WRAP_ENABLED = bool(cfg.get("CAPTION_WRAP_ENABLED", True))
    CAPTION_WRAP_CHAR_WIDTH = int(cfg.get("CAPTION_WRAP_CHAR_WIDTH", 35))
    CAPTION_MAX_NAME_LINES = cfg.get("CAPTION_MAX_NAME_LINES", 2)

    SUBPLOT_WSPACE = float(cfg.get("SUBPLOT_WSPACE", 0.20))
    SUBPLOT_HSPACE = float(cfg.get("SUBPLOT_HSPACE", 0.30))

    SHOW_COLORBAR = bool(cfg.get("SHOW_COLORBAR", True))
    CUSTOM_CMAPS_XLSX = cfg.get("CUSTOM_CMAPS_XLSX", "")

    DENSITY_CMAP_NAME = cfg.get("DENSITY_CMAP_NAME", "plasma_r")

    COLOR_MODE = (cfg.get("COLOR_MODE", "density") or "density").lower().strip()
    ENRICHMENT_SCALE = (cfg.get("ENRICHMENT_SCALE", "ratio") or "ratio").lower().strip()
    ENRICHMENT_STYLE = (cfg.get("ENRICHMENT_STYLE", "unilateral_diverging") or "unilateral_diverging").lower().strip()
    ENRICHMENT_CMAP_NAME = cfg.get("ENRICHMENT_CMAP_NAME", "RdBu_r")
    ENRICHMENT_VMAX = cfg.get("ENRICHMENT_VMAX", None)
    ENRICHMENT_PERCENTILE = float(cfg.get("ENRICHMENT_PERCENTILE", 99.0))
    ENRICHMENT_SYMMETRIC = bool(cfg.get("ENRICHMENT_SYMMETRIC", True))
    ENRICHMENT_EPS = float(cfg.get("ENRICHMENT_EPS", 1e-12))
    ALL_GENES_PANEL_MODE = (cfg.get("ALL_GENES_PANEL_MODE", "background") or "background").lower().strip()

    # Grey-below-neutral toggle (keep it, but you can set False in the pipeline script)
    ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER = bool(cfg.get("ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER", True))
    ENRICHMENT_NEUTRAL_RGBA = cfg.get("ENRICHMENT_NEUTRAL_RGBA", (0.85, 0.85, 0.85, 1.0))

    # NEW (optional): controls unilateral colormap handling
    ENRICHMENT_UNILATERAL_CMAP_MODE = (cfg.get("ENRICHMENT_UNILATERAL_CMAP_MODE", "auto") or "auto").lower().strip()
    ENRICHMENT_UNILATERAL_SIDE = (cfg.get("ENRICHMENT_UNILATERAL_SIDE", "high") or "high").lower().strip()

    SAVE_FIG = bool(cfg.get("SAVE_FIG", True))
    PNG_DPI = int(cfg.get("PNG_DPI", 220))

    # ---- open workbook
    xlsx = EXCEL_PATH
    if not xlsx or not os.path.isfile(xlsx):
        if USE_FILE_DIALOG:
            xlsx = _pick_excel()
        if not xlsx or not os.path.isfile(xlsx):
            print("No valid Excel selected. Aborting.")
            return

    print(f"[INFO] Workbook: {xlsx}")
    excel = pd.ExcelFile(xlsx)

    # Coordinates
    umap_df = excel.parse(UMAP_SHEET, header=None)
    if umap_df.shape[1] < 3:
        raise ValueError(f'"{UMAP_SHEET}" must have 3 columns: [locus_tag, X, Y]')
    umap_df = umap_df.iloc[:, :3]
    umap_df.columns = ["locus_tag", "X", "Y"]
    umap_df["locus_tag"] = umap_df["locus_tag"].astype(str).str.strip()
    umap_df["X"] = pd.to_numeric(umap_df["X"], errors="coerce")
    umap_df["Y"] = pd.to_numeric(umap_df["Y"], errors="coerce")
    umap_df = umap_df.dropna(subset=["X", "Y"])

    # Binary
    data_df = excel.parse(DATA_SHEET, dtype=object)
    if data_df.shape[1] < 3:
        raise ValueError(f'"{DATA_SHEET}" must have ≥3 columns (locus, gene, data...)"')
    locus_col = data_df.columns[0]
    data_cols = list(data_df.columns[2:])

    if INCLUDE_COLUMNS is not None:
        wanted = {c.lower() for c in INCLUDE_COLUMNS}
        data_cols = [c for c in data_cols if c.lower() in wanted]
    if not data_cols:
        print("[WARN] No data columns selected; nothing to plot.")
        return

    data_df[locus_col] = data_df[locus_col].astype(str).str.strip()

    # canonical merge key
    umap_df["tag_key"] = umap_df["locus_tag"].map(_canon_tag)
    data_df["tag_key"] = data_df[locus_col].map(_canon_tag)

    umap_df = umap_df.drop_duplicates(subset=["tag_key"])
    data_df = data_df.drop_duplicates(subset=["tag_key"])

    merged = pd.merge(
        umap_df[["tag_key", "X", "Y"]],
        data_df[["tag_key"] + data_cols],
        on="tag_key",
        how="inner",
    )

    n_all_coords = umap_df.shape[0]
    n_merged = merged.shape[0]
    print(f"[INFO] Matched {n_merged}/{n_all_coords} coordinate points to data rows ({n_all_coords - n_merged} coords without binary data).")

    # global cloud
    Xg = umap_df["X"].to_numpy(float)
    Yg = umap_df["Y"].to_numpy(float)
    all_xy = np.column_stack([Xg, Yg])
    xlim, ylim = _auto_axis_limits(Xg, Yg)

    # cluster subset
    x_all = merged["X"].to_numpy(float)
    y_all = merged["Y"].to_numpy(float)

    if HIGHLIGHT_COLUMNS is None:
        highlight_cols = data_cols
    else:
        wanted = {c.lower() for c in HIGHLIGHT_COLUMNS}
        highlight_cols = [c for c in data_cols if c.lower() in wanted]
    if not highlight_cols:
        print("[WARN] No highlighted columns found (check HIGHLIGHT_COLUMNS).")
        return

    # cluster sizes/points
    cluster_sizes, cluster_points = {}, {}
    for col in highlight_cols:
        s = pd.to_numeric(merged[col], errors="coerce").fillna(0.0).to_numpy(float)
        mask = s > 0.0
        n = int(mask.sum())
        cluster_sizes[col] = n
        cluster_points[col] = np.column_stack([x_all[mask], y_all[mask]]) if n > 0 else np.zeros((0, 2))
        print(f"[INFO] Cluster '{col}': n = {n}")

    # colormaps
    custom_maps = _load_custom_colormaps(CUSTOM_CMAPS_XLSX)
    cmap_density = _choose_cmap(DENSITY_CMAP_NAME, custom_maps)
    cmap_enrich_base = _choose_cmap(ENRICHMENT_CMAP_NAME, custom_maps)

    # enrichment scaling cache
    enrich_norm = enrich_cmap = enrich_center = sm_for_colorbar = None
    dens_all_cache = {}
    pooled_vals = []

    if COLOR_MODE == "enrichment":
        for col in highlight_cols:
            n = cluster_sizes[col]
            if n <= 0:
                continue
            preset = _get_cluster_preset(n, CLUSTER_SCATTER_PRESETS)
            if preset is None:
                continue
            nbins = int(preset["DENSITY_NBINS"])
            sigma = float(preset["DENSITY_SIGMA"])
            key = (nbins, sigma)
            if key not in dens_all_cache:
                xedges, yedges = _edges_from_limits(xlim, ylim, nbins)
                field_all = _density_field(all_xy, xedges, yedges, sigma=sigma)
                dens_all_cache[key] = (xedges, yedges, field_all)

            xedges, yedges, field_all = dens_all_cache[key]
            vals = _compute_enrichment_values(
                cluster_xy=cluster_points[col],
                all_xy=all_xy,
                xedges=xedges, yedges=yedges,
                sigma=sigma,
                eps=ENRICHMENT_EPS,
                scale=ENRICHMENT_SCALE,
                field_all_cache=field_all,
            )
            if vals.size:
                pooled_vals.append(vals)

        pooled = np.concatenate(pooled_vals) if pooled_vals else np.array([], dtype=float)

        if ENRICHMENT_VMAX is None:
            vmax = _auto_enrichment_vmax(pooled, ENRICHMENT_SCALE, ENRICHMENT_PERCENTILE, ENRICHMENT_STYLE)
        else:
            vmax = float(ENRICHMENT_VMAX)

        enrich_norm, enrich_cmap, enrich_center = _enrichment_norm_and_cmap(
            vmax=vmax,
            base_cmap=cmap_enrich_base,
            scale=ENRICHMENT_SCALE,
            style=ENRICHMENT_STYLE,
            symmetric=ENRICHMENT_SYMMETRIC,
            unilateral_side=ENRICHMENT_UNILATERAL_SIDE,
            unilateral_cmap_mode=ENRICHMENT_UNILATERAL_CMAP_MODE,
        )

        # Colorbar mappable must use the SAME norm+cmap
        sm_for_colorbar = ScalarMappable(norm=enrich_norm, cmap=enrich_cmap)
        # Give it a tiny array so matplotlib fully initializes the gradient
        sm_for_colorbar.set_array(np.array([getattr(enrich_norm, "vmin", 0.0), getattr(enrich_norm, "vmax", 1.0)], dtype=float))

        print(
            f"[INFO] Enrichment scaling: scale={ENRICHMENT_SCALE}, style={ENRICHMENT_STYLE}, "
            f"vmax={vmax:.4g}, unilateral_cmap_mode={ENRICHMENT_UNILATERAL_CMAP_MODE}, unilateral_side={ENRICHMENT_UNILATERAL_SIDE}"
        )

    # layout
    n_panels = (1 + len(highlight_cols)) if PLOT_MODE == "grid" else 1
    if PLOT_MODE == "grid":
        nrows = min(MAX_NROWS, int(np.ceil(np.sqrt(n_panels))))
        nrows = max(1, nrows)
        ncols = int(np.ceil(n_panels / nrows))
    elif PLOT_MODE == "overlay":
        nrows, ncols = 1, 1
    else:
        raise ValueError("PLOT_MODE must be 'grid' or 'overlay'.")

    fig_w = PANEL_W_IN * ncols
    fig_h = PANEL_H_IN * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), squeeze=False)
    plt.subplots_adjust(wspace=SUBPLOT_WSPACE, hspace=SUBPLOT_HSPACE)

    def _panel_axes(i):
        r, c = divmod(i, ncols)
        return axes[r, c]

    # plot: grid mode
    if PLOT_MODE == "grid":
        ax0 = _panel_axes(0)

        # all genes panel
        if all_xy.shape[0] > 0:
            if COLOR_MODE == "enrichment" and ALL_GENES_PANEL_MODE == "background":
                ax0.scatter(all_xy[:, 0], all_xy[:, 1], s=BACKGROUND_SIZE, c=BACKGROUND_COLOR, alpha=0.35, linewidths=0.0)
            else:
                colors_all = _compute_density_colors(
                    all_xy,
                    nbins=ALL_SCATTER_DENSITY_NBINS,
                    sigma=ALL_SCATTER_DENSITY_SIGMA,
                    use_log=ALL_SCATTER_DENSITY_USE_LOG,
                    min_rel=ALL_SCATTER_DENSITY_MIN_REL,
                    max_rel=ALL_SCATTER_DENSITY_MAX_REL,
                    cmap=cmap_density
                )
                ax0.scatter(all_xy[:, 0], all_xy[:, 1], s=ALL_SCATTER_POINT_SIZE, c=colors_all,
                            alpha=ALL_SCATTER_ALPHA, linewidths=ALL_SCATTER_EDGE_WIDTH)

        ax0.set_xlim(xlim); ax0.set_ylim(ylim)
        ax0.set_xlabel(X_LABEL); ax0.set_ylabel(Y_LABEL)
        if DRAW_GRID:
            _draw_grid(ax0, xlim, ylim, major_step=GRID_MAJOR_STEP, major_alpha=GRID_MAJOR_ALPHA,
                       major_lw=GRID_MAJOR_LW, major_color=GRID_MAJOR_COLOR)
        ax0.set_title(
            _wrap_cluster_title("All genes", all_xy.shape[0], enabled=CAPTION_WRAP_ENABLED,
                                width=CAPTION_WRAP_CHAR_WIDTH, max_name_lines=CAPTION_MAX_NAME_LINES),
            fontsize=CAPTION_SIZE, fontweight=CAPTION_WEIGHT, pad=CAPTION_PAD
        )

        # cluster panels
        for i, col in enumerate(highlight_cols, start=1):
            ax = _panel_axes(i)
            n = cluster_sizes[col]
            cluster_xy = cluster_points[col]

            if SHOW_BACKGROUND_SMALL:
                ax.scatter(all_xy[:, 0], all_xy[:, 1], s=BACKGROUND_SIZE, c=BACKGROUND_COLOR,
                           alpha=BACKGROUND_ALPHA, linewidths=BACKGROUND_LW, edgecolors="none", zorder=1)

            if n > 0:
                preset = _get_cluster_preset(n, CLUSTER_SCATTER_PRESETS) or {}
                nbins = int(preset.get("DENSITY_NBINS", ALL_SCATTER_DENSITY_NBINS))
                sigma = float(preset.get("DENSITY_SIGMA", ALL_SCATTER_DENSITY_SIGMA))

                if COLOR_MODE == "density":
                    colors = _compute_density_colors(
                        cluster_xy,
                        nbins=nbins,
                        sigma=sigma,
                        use_log=bool(preset.get("DENSITY_USE_LOG", True)),
                        min_rel=float(preset.get("DENSITY_MIN_REL", 0.0)),
                        max_rel=float(preset.get("DENSITY_MAX_REL", 1.0)),
                        cmap=cmap_density,
                    )
                else:
                    key = (nbins, sigma)
                    if key not in dens_all_cache:
                        xedges, yedges = _edges_from_limits(xlim, ylim, nbins)
                        field_all = _density_field(all_xy, xedges, yedges, sigma=sigma)
                        dens_all_cache[key] = (xedges, yedges, field_all)
                    xedges, yedges, field_all = dens_all_cache[key]

                    vals = _compute_enrichment_values(
                        cluster_xy=cluster_xy,
                        all_xy=all_xy,
                        xedges=xedges, yedges=yedges,
                        sigma=sigma,
                        eps=ENRICHMENT_EPS,
                        scale=ENRICHMENT_SCALE,
                        field_all_cache=field_all
                    )

                    colors = _map_enrichment_to_colors(
                        vals=vals,
                        norm=enrich_norm,
                        cmap=enrich_cmap,
                        center_value=enrich_center,
                        force_neutral=ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER,
                        neutral_rgba=ENRICHMENT_NEUTRAL_RGBA
                    )

                ax.scatter(
                    cluster_xy[:, 0], cluster_xy[:, 1],
                    s=float(preset.get("POINT_SIZE", 10.0)),
                    c=colors,
                    alpha=float(preset.get("ALPHA", HIGHLIGHT_ALPHA)),
                    marker=HIGHLIGHT_MARKER,
                    edgecolors=HIGHLIGHT_EDGE_COLOR,
                    linewidths=float(preset.get("EDGE_WIDTH", HIGHLIGHT_EDGE_WIDTH)),
                    zorder=2,
                )

            ax.set_xlim(xlim); ax.set_ylim(ylim)
            ax.set_xlabel(X_LABEL); ax.set_ylabel(Y_LABEL)
            if DRAW_GRID:
                _draw_grid(ax, xlim, ylim, major_step=GRID_MAJOR_STEP, major_alpha=GRID_MAJOR_ALPHA,
                           major_lw=GRID_MAJOR_LW, major_color=GRID_MAJOR_COLOR)
            ax.set_title(
                _wrap_cluster_title(col, n, enabled=CAPTION_WRAP_ENABLED,
                                    width=CAPTION_WRAP_CHAR_WIDTH, max_name_lines=CAPTION_MAX_NAME_LINES),
                fontsize=CAPTION_SIZE, fontweight=CAPTION_WEIGHT, pad=CAPTION_PAD
            )

        # hide extra axes
        total_slots = nrows * ncols
        for j in range(n_panels, total_slots):
            _panel_axes(j).axis("off")

    else:
        # overlay mode
        ax = axes[0, 0]
        if SHOW_BACKGROUND_SMALL:
            ax.scatter(all_xy[:, 0], all_xy[:, 1], s=BACKGROUND_SIZE, c=BACKGROUND_COLOR,
                       alpha=BACKGROUND_ALPHA, linewidths=BACKGROUND_LW, edgecolors="none", zorder=1)

        for col in highlight_cols:
            n = cluster_sizes[col]
            if n <= 0:
                continue
            cluster_xy = cluster_points[col]
            preset = _get_cluster_preset(n, CLUSTER_SCATTER_PRESETS) or {}
            nbins = int(preset.get("DENSITY_NBINS", ALL_SCATTER_DENSITY_NBINS))
            sigma = float(preset.get("DENSITY_SIGMA", ALL_SCATTER_DENSITY_SIGMA))

            if COLOR_MODE == "density":
                colors = _compute_density_colors(
                    cluster_xy,
                    nbins=nbins,
                    sigma=sigma,
                    use_log=bool(preset.get("DENSITY_USE_LOG", True)),
                    min_rel=float(preset.get("DENSITY_MIN_REL", 0.0)),
                    max_rel=float(preset.get("DENSITY_MAX_REL", 1.0)),
                    cmap=cmap_density,
                )
            else:
                key = (nbins, sigma)
                if key not in dens_all_cache:
                    xedges, yedges = _edges_from_limits(xlim, ylim, nbins)
                    field_all = _density_field(all_xy, xedges, yedges, sigma=sigma)
                    dens_all_cache[key] = (xedges, yedges, field_all)
                xedges, yedges, field_all = dens_all_cache[key]

                vals = _compute_enrichment_values(
                    cluster_xy=cluster_xy,
                    all_xy=all_xy,
                    xedges=xedges, yedges=yedges,
                    sigma=sigma,
                    eps=ENRICHMENT_EPS,
                    scale=ENRICHMENT_SCALE,
                    field_all_cache=field_all
                )

                colors = _map_enrichment_to_colors(
                    vals=vals,
                    norm=enrich_norm,
                    cmap=enrich_cmap,
                    center_value=enrich_center,
                    force_neutral=ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER,
                    neutral_rgba=ENRICHMENT_NEUTRAL_RGBA
                )

            ax.scatter(
                cluster_xy[:, 0], cluster_xy[:, 1],
                s=float(preset.get("POINT_SIZE", 10.0)),
                c=colors,
                alpha=float(preset.get("ALPHA", HIGHLIGHT_ALPHA)),
                marker=HIGHLIGHT_MARKER,
                edgecolors=HIGHLIGHT_EDGE_COLOR,
                linewidths=float(preset.get("EDGE_WIDTH", HIGHLIGHT_EDGE_WIDTH)),
                zorder=2,
            )

        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.set_xlabel(X_LABEL); ax.set_ylabel(Y_LABEL)
        if DRAW_GRID:
            _draw_grid(ax, xlim, ylim, major_step=GRID_MAJOR_STEP, major_alpha=GRID_MAJOR_ALPHA,
                       major_lw=GRID_MAJOR_LW, major_color=GRID_MAJOR_COLOR)
        ax.set_title(
            f"Cluster scatter overlay — mode='overlay'\nCOLOR_MODE={COLOR_MODE} | highlighted: {', '.join(highlight_cols)}",
            fontsize=CAPTION_SIZE, fontweight=CAPTION_WEIGHT, pad=CAPTION_PAD
        )

    # colorbar
    if SHOW_COLORBAR:
        if COLOR_MODE == "enrichment":
            label = "Fold enrichment (cluster vs all)" if ENRICHMENT_SCALE == "ratio" else "log2 enrichment (cluster vs all)"
            cbar = fig.colorbar(sm_for_colorbar, ax=axes, fraction=0.02, pad=0.02)
            cbar.set_label(label, fontsize=CAPTION_SIZE, fontweight=CAPTION_WEIGHT)
        else:
            sm = ScalarMappable(cmap=cmap_density, norm=Normalize(vmin=0.0, vmax=1.0))
            sm.set_array(np.array([0.0, 1.0], dtype=float))
            cbar = fig.colorbar(sm, ax=axes, fraction=0.02, pad=0.02)
            cbar.set_label("relative local density", fontsize=CAPTION_SIZE, fontweight=CAPTION_WEIGHT)

    # export
    if SAVE_FIG:
        base = os.path.splitext(os.path.basename(xlsx))[0]
        mode_tag = f"_ENRICH_{ENRICHMENT_SCALE}_{ENRICHMENT_STYLE}" if COLOR_MODE == "enrichment" else "_DENSITY"
        out = os.path.join(os.path.dirname(xlsx), f"{base}_UMAP_local{mode_tag}.png")
        fig.savefig(out, dpi=PNG_DPI, bbox_inches="tight")
        print(f"[INFO] Saved: {out}")

    plt.show()
