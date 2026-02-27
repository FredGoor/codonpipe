"""Cluster overlay plotting for 2D embeddings (UMAP / t-SNE / PCA).

Inputs:
  - Excel workbook produced by Clustering_Pipeline.py
  - Coordinate sheet (default: "coordinates")
  - Binary membership sheet (default: "Binary")

Outputs:
  - A multi-panel figure showing one panel per cluster,
    with optional density or enrichment coloring.
  - Optionally, a genome-axis gene-cluster localization heatmap.
"""


# ======================== USER SETTINGS ===========================

EXCEL_PATH   = r""
UMAP_SHEET   = "coordinates"
DATA_SHEET   = "Binary"

USE_FILE_DIALOG = True

INCLUDE_COLUMNS  = None
HIGHLIGHT_COLUMNS = None

PLOT_MODE = "grid"

MAX_NROWS    = 4
PANEL_W_IN   = 5.0
PANEL_H_IN   = 5.0

X_LABEL = "t-SNE 1"
Y_LABEL = "t-SNE 2"

SHOW_BACKGROUND_SMALL   = True
BACKGROUND_COLOR        = "grey"
BACKGROUND_ALPHA        = 0.18
BACKGROUND_SIZE         = 8
BACKGROUND_LW           = 0.0

HIGHLIGHT_MARKER         = "o"
HIGHLIGHT_EDGE_COLOR     = "black"
HIGHLIGHT_EDGE_WIDTH     = 0.6
HIGHLIGHT_ALPHA          = 0.95

ALL_SCATTER_POINT_SIZE      = 3.5
ALL_SCATTER_ALPHA           = 0.7
ALL_SCATTER_EDGE_WIDTH      = 0.15
ALL_SCATTER_DENSITY_NBINS   = 150
ALL_SCATTER_DENSITY_SIGMA   = 4.0
ALL_SCATTER_DENSITY_USE_LOG = True
ALL_SCATTER_DENSITY_MIN_REL = 0.00
ALL_SCATTER_DENSITY_MAX_REL = 1.00

CLUSTER_SCATTER_PRESETS = [
    {"MIN_N": 0,    "MAX_N": 100,  "POINT_SIZE": 40.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.6,  "DENSITY_NBINS": 50,  "DENSITY_SIGMA": 2.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 100,  "MAX_N": 250,  "POINT_SIZE": 35.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.6,  "DENSITY_NBINS": 80,  "DENSITY_SIGMA": 2.5, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 250,  "MAX_N": 500,  "POINT_SIZE": 20.0, "ALPHA": 0.95, "EDGE_WIDTH": 0.3,  "DENSITY_NBINS": 110, "DENSITY_SIGMA": 3.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 500,  "MAX_N": 1500, "POINT_SIZE": 14,   "ALPHA": 0.9,  "EDGE_WIDTH": 0.2,  "DENSITY_NBINS": 150, "DENSITY_SIGMA": 4.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 1500, "MAX_N": 2500, "POINT_SIZE": 6.0,  "ALPHA": 0.85, "EDGE_WIDTH": 0.10, "DENSITY_NBINS": 180, "DENSITY_SIGMA": 4.5, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
    {"MIN_N": 2500, "MAX_N": None, "POINT_SIZE": 3.5,  "ALPHA": 0.8,  "EDGE_WIDTH": 0.08, "DENSITY_NBINS": 150, "DENSITY_SIGMA": 4.0, "DENSITY_USE_LOG": True, "DENSITY_MIN_REL": 0.00, "DENSITY_MAX_REL": 1.00},
]

DRAW_GRID        = True
GRID_MAJOR_STEP  = None
GRID_MAJOR_ALPHA = 0.55
GRID_MAJOR_LW    = 1.0
GRID_MAJOR_COLOR = "black"

CAPTION_SIZE         = 12
CAPTION_WEIGHT       = "bold"
CAPTION_PAD          = 8
CAPTION_WRAP_ENABLED     = True
CAPTION_WRAP_CHAR_WIDTH  = 35
CAPTION_MAX_NAME_LINES   = 2

SUBPLOT_WSPACE = 0.20
SUBPLOT_HSPACE = 0.30

SHOW_COLORBAR = True

CUSTOM_CMAPS_XLSX = ""
DENSITY_CMAP_NAME = "plasma_r"

# ===================== Gene-cluster localization heatmap =====================
# Uses the workbook sheet produced by Clustering_Pipeline.py (default: 'Locus Tags').
# The sheet must be:
#   - Column 1: reordered genome locus tags (genome axis)
#   - Other columns: each column lists locus tags belonging to that cluster
RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP = True
LOCUS_TAGS_SHEET = "Locus Tags"

GCHM_COLORMAP = "plasma"          # can be 'plasma', 'viridis', or a custom colormap sheet name
GCHM_CUSTOM_CMAPS_XLSX = CUSTOM_CMAPS_XLSX  # "" disables custom colormaps

GCHM_SIGMA = 10
GCHM_SPREAD_FACTOR = 5
GCHM_HEIGHT_PER_CLUSTER = 0.3
GCHM_LABEL_FONTSIZE = 10
GCHM_DPI = 300

# Colormap range in relative density units (0..1)
GCHM_CMAP_MIN_REL = 0.2
GCHM_CMAP_MAX_REL = 1.0

# Output file (saved next to the workbook)
GCHM_OUTPUT_FILENAME = "gene_cluster_heatmap_KS.png"

# If running from Clustering_Pipeline, you may prefer show=False to avoid popping a second window
GCHM_SHOW_FIG = True

# ---- Density vs enrichment coloring ----
COLOR_MODE = "enrichment"   # "density" or "enrichment"

# ---- Enrichment coloring options (used only when COLOR_MODE = "enrichment") ----
ENRICHMENT_SCALE = "log2"  # "ratio" or "log2"
ENRICHMENT_STYLE = "unilateral_diverging"
ENRICHMENT_CMAP_NAME = "plasma_r"
ENRICHMENT_VMAX = 3
ENRICHMENT_PERCENTILE = 99.0
ENRICHMENT_SYMMETRIC = True
ENRICHMENT_EPS = 1e-12

# CHANGED: force the top-left "All genes" panel to use DENSITY coloring,
# even when COLOR_MODE="enrichment".
ALL_GENES_PANEL_MODE = "density"  # "density" or "background"

ENRICHMENT_FORCE_NEUTRAL_BELOW_CENTER = False
ENRICHMENT_NEUTRAL_RGBA = (0.85, 0.85, 0.85, 1.0)

# ---- Output ----
SAVE_FIG = True
PNG_DPI = 220

# ---- Debugging ----
DEBUG_PRINT_CONFIG = False  # Set False to silence debug prints

# ================================================================

import os, sys

# Optional: allow Clustering_Pipeline.py to enforce consistent colormaps across the whole pipeline
# (it sets CODONPIPE_DENSITY_CMAP / CODONPIPE_ENRICHMENT_CMAP before importing this script).
_env_density = os.environ.get("CODONPIPE_DENSITY_CMAP", "").strip()
_env_enrich  = os.environ.get("CODONPIPE_ENRICHMENT_CMAP", "").strip()

if _env_density:
    DENSITY_CMAP_NAME = _env_density
if _env_enrich:
    ENRICHMENT_CMAP_NAME = _env_enrich
elif _env_density:
    # If only one cmap is specified, use it for both density + enrichment
    ENRICHMENT_CMAP_NAME = _env_density


# Optional: allow Clustering_Pipeline.py to drive gene-cluster heatmap parameters via env vars
# (keeps the bridge minimal; Plotting_Pipeline remains runnable standalone).
def _env_bool(key: str, default: bool) -> bool:
    v = os.environ.get(key, "").strip().lower()
    if v == "":
        return default
    return v in ("1", "true", "yes", "y", "on")

def _env_float(key: str, default: float) -> float:
    v = os.environ.get(key, "").strip()
    if v == "":
        return default
    try:
        return float(v)
    except Exception:
        return default

def _env_int(key: str, default: int) -> int:
    v = os.environ.get(key, "").strip()
    if v == "":
        return default
    try:
        return int(float(v))
    except Exception:
        return default

RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP = _env_bool("CODONPIPE_GCHM_ENABLE", RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP)
LOCUS_TAGS_SHEET = os.environ.get("CODONPIPE_GCHM_SHEET", LOCUS_TAGS_SHEET).strip() or LOCUS_TAGS_SHEET
GCHM_COLORMAP = os.environ.get("CODONPIPE_GCHM_COLORMAP", GCHM_COLORMAP).strip() or GCHM_COLORMAP
GCHM_CUSTOM_CMAPS_XLSX = os.environ.get("CODONPIPE_GCHM_CUSTOM_CMAPS_XLSX", GCHM_CUSTOM_CMAPS_XLSX).strip()
GCHM_SIGMA = _env_float("CODONPIPE_GCHM_SIGMA", GCHM_SIGMA)
GCHM_SPREAD_FACTOR = _env_float("CODONPIPE_GCHM_SPREAD_FACTOR", GCHM_SPREAD_FACTOR)
GCHM_HEIGHT_PER_CLUSTER = _env_float("CODONPIPE_GCHM_HEIGHT_PER_CLUSTER", GCHM_HEIGHT_PER_CLUSTER)
GCHM_LABEL_FONTSIZE = _env_int("CODONPIPE_GCHM_LABEL_FONTSIZE", GCHM_LABEL_FONTSIZE)
GCHM_DPI = _env_int("CODONPIPE_GCHM_DPI", GCHM_DPI)
GCHM_CMAP_MIN_REL = _env_float("CODONPIPE_GCHM_CMAP_MIN_REL", GCHM_CMAP_MIN_REL)
GCHM_CMAP_MAX_REL = _env_float("CODONPIPE_GCHM_CMAP_MAX_REL", GCHM_CMAP_MAX_REL)
GCHM_OUTPUT_FILENAME = os.environ.get("CODONPIPE_GCHM_OUTPUT_FILENAME", GCHM_OUTPUT_FILENAME).strip() or GCHM_OUTPUT_FILENAME
GCHM_SHOW_FIG = _env_bool("CODONPIPE_GCHM_SHOW", GCHM_SHOW_FIG)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from codonpipe.density_plot_core import run_density_dotplot
from codonpipe.legend import append_plotting_sections
from codonpipe.gene_cluster_heatmap import run_gene_cluster_localization_heatmap_from_workbook


def main():
    # ------------------------------------------------------------------
    # Auto-sync axis labels with the dimensionality-reduction method used
    # upstream (Clustering_Pipeline sets CODONPIPE_DIMRED_METHOD).
    # ------------------------------------------------------------------
    global X_LABEL, Y_LABEL

    dimred = os.environ.get("CODONPIPE_DIMRED_METHOD", "").strip().lower()
    x_from_env = os.environ.get("CODONPIPE_X_LABEL", "").strip()
    y_from_env = os.environ.get("CODONPIPE_Y_LABEL", "").strip()

    # Allow explicit axis labels via env vars (highest priority)
    if x_from_env:
        X_LABEL = x_from_env
    if y_from_env:
        Y_LABEL = y_from_env

    # Otherwise, auto-label based on dimred method, but only if the labels
    # are still the default t-SNE ones (so we don't clobber custom labels).
    if dimred and (not x_from_env) and (not y_from_env):
        default_x = str(X_LABEL).strip().lower() in {"t-sne 1", "tsne 1", "tsne1", "tsne_1"}
        default_y = str(Y_LABEL).strip().lower() in {"t-sne 2", "tsne 2", "tsne2", "tsne_2"}
        if default_x and default_y:
            if dimred == "tsne":
                X_LABEL, Y_LABEL = "tSNE1", "tSNE2"
            elif dimred == "umap":
                X_LABEL, Y_LABEL = "UMAP1", "UMAP2"
            elif dimred == "pca":
                X_LABEL, Y_LABEL = "PC1", "PC2"

    # Collect all UPPERCASE globals as config for the core plotter
    cfg = {k: v for k, v in globals().items() if k.isupper()}

    if cfg.get("DEBUG_PRINT_CONFIG", False):
        print("\n========== [PLOT DEBUG] Plotting_Pipeline.py ==========")
        print(f"[PLOT DEBUG] Script path: {os.path.abspath(__file__)}")
        print(f"[PLOT DEBUG] EXCEL_PATH={cfg.get('EXCEL_PATH')}")
        print(f"[PLOT DEBUG] UMAP_SHEET={cfg.get('UMAP_SHEET')} | DATA_SHEET={cfg.get('DATA_SHEET')}")
        print(f"[PLOT DEBUG] COLOR_MODE={cfg.get('COLOR_MODE')}")
        print(f"[PLOT DEBUG] DENSITY_CMAP_NAME={cfg.get('DENSITY_CMAP_NAME')}")
        print(
            "[PLOT DEBUG] ENRICH settings: "
            f"SCALE={cfg.get('ENRICHMENT_SCALE')}, "
            f"STYLE={cfg.get('ENRICHMENT_STYLE')}, "
            f"CMAP={cfg.get('ENRICHMENT_CMAP_NAME')}, "
            f"VMAX={cfg.get('ENRICHMENT_VMAX')}, "
            f"PCTL={cfg.get('ENRICHMENT_PERCENTILE')}, "
            f"SYM={cfg.get('ENRICHMENT_SYMMETRIC')}, "
            f"EPS={cfg.get('ENRICHMENT_EPS')}"
        )
        print(f"[PLOT DEBUG] ALL_GENES_PANEL_MODE={cfg.get('ALL_GENES_PANEL_MODE')}")
        print("=======================================================\n")

    # Run plotting
    run_density_dotplot(cfg)

    # Optional: gene-cluster localization heatmap along genome axis
    if RUN_GENE_CLUSTER_LOCALIZATION_HEATMAP and cfg.get("EXCEL_PATH"):
        try:
            wb = cfg.get("EXCEL_PATH")
            out_png = os.path.join(os.path.dirname(wb), GCHM_OUTPUT_FILENAME)
            print("[INFO] Generating gene-cluster localization heatmap (KS)...")
            run_gene_cluster_localization_heatmap_from_workbook(
                workbook_path=wb,
                sheet_name=LOCUS_TAGS_SHEET,
                output_file=out_png,
                colormap=GCHM_COLORMAP,
                custom_cmaps_xlsx=GCHM_CUSTOM_CMAPS_XLSX,
                sigma=GCHM_SIGMA,
                spread_factor=GCHM_SPREAD_FACTOR,
                height_per_cluster=GCHM_HEIGHT_PER_CLUSTER,
                label_fontsize=GCHM_LABEL_FONTSIZE,
                dpi=GCHM_DPI,
                cmap_min_rel=GCHM_CMAP_MIN_REL,
                cmap_max_rel=GCHM_CMAP_MAX_REL,
                show=GCHM_SHOW_FIG,
            )
            print(f"[INFO] Gene-cluster heatmap saved to:\n  {out_png}")
        except Exception as e:
            print(f"[WARN] Gene-cluster heatmap failed (skipping): {e}")


    # Append to the three text files (if Clustering_Pipeline provided CODONPIPE_OUT_BASE)
    out_base = os.environ.get("CODONPIPE_OUT_BASE", "")
    if out_base and cfg.get("EXCEL_PATH"):
        try:
            append_plotting_sections(out_base=out_base, cfg=cfg, workbook_path=cfg.get("EXCEL_PATH"))
            print(f"[INFO] Appended plotting sections to TXT files with base:\n  {out_base}")
        except Exception as e:
            print(f"[WARN] Could not append plotting sections: {e}")


if __name__ == "__main__":
    main()
