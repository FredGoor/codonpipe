#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File: Clustering_Pipeline.py
"""
Codon-usage clustering pipeline (Excel + FASTA inputs).

This script:
  1) Reads a per-gene codon usage table (Excel).
  2) Computes a 2D embedding (UMAP/t-SNE/PCA) and clusters genes in the embedding.
  3) Produces an output Excel workbook with coordinates, curated cluster membership, and FASTA-derived metrics.
  4) Optionally runs the plotting script (Plotting_Pipeline.py) on the generated workbook.

Repository layout (recommended):
  - Clustering_Pipeline.py
  - Plotting_Pipeline.py
  - codonpipe/ (modules)
"""


import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))  # repository root (this file location)
ROOT = HERE
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from codonpipe.clustering import run_codon_clustering, choose_cluster_file
from codonpipe.fasta_metrics import (
    auto_find_fasta, choose_fasta, build_locus_index,
    canonicalize_id, enforce_unique_after_canon,
    canonicalize_cluster_df, canonicalize_generic_map,
    build_metrics_table, build_summary, build_quantitative_bis
)
from codonpipe.excel_outputs import (
    read_locus_clusters, build_coordinates_df,
    build_locus_tags_sheet, build_binary_sheet,
    build_meta_table, write_per_cluster_workbook,
    write_codon_usage_by_cluster_workbook
)
from codonpipe.ks2d import compute_2d_ks_for_clusters
from codonpipe.density_bridge import run_density_plot_script
from codonpipe.legend import write_all_texts


# =============================================================================
# =========================
# ==== USER SETTINGS ======
# =========================
# =============================================================================

SET = dict(
    # ----------------------------
    # Paths / dataset context
    default_root="",

    # ----------------------------
    # Input mode
    # ----------------------------
    #   'fasta' : user selects a CDS FASTA; codon table is computed on the fly (recommended)
    #   'excel' : user selects a precomputed codon-usage Excel table (+ GeneIDs Excel)
    input_mode='fasta',

    # FASTA-mode options (used only when input_mode='fasta')
    fasta_path='',
    fasta_row_id_mode='primary',         # 'primary' | 'old' | 'locus'
    fasta_trim_to_multiple_of_3=True,
    fasta_include_stops=None,            # None -> infer from codon_set (64_all keeps stops)
    fasta_codon_table_mode='freq',       # 'freq' (default) or 'abs'
    fasta_write_intermediate_excels=False,

    # ----------------------------
    # Codon usage basis + codon set
    # ----------------------------
    usage_basis='acu',  #   Choices : 'RCU', 'ACU', 'AA'. 'RCU' = relative codon usage among synonyms (per gene) - 'ACU' = absolute codon usage / counts-like metrics (per gene) - 'AA'  = amino-acid usage features (per gene)
    codon_set='61_noSTOP',  #   '59_noSTOP_MW', '61_noSTOP', '64_withSTOP', etc.

    # ----------------------------
    # Dimensionality reduction + clustering
    # ----------------------------
    dimred_method='umap',   #   Choices: 'umap', 'tsne', 'pca'
    cluster_method='kmedoids',   #   Choices: 'hierarchical', 'kmedoids', 'spectral', 'dbscan'
    do_2d_ks=False ,   #   If True: compute 2D KS comparisons among clusters in the embedding space

    # ----------------------------
    # UMAP parameters
    # ----------------------------
    umap_neighbors=10,     # 10   - Range 5-100
    umap_min_dist=0.01,    # 0.01 - Range 0-0.5
    umap_metric='cosine',  # cosine - Choices (typical UMAP metrics): 'euclidean', 'cosine', 'manhattan', 'chebyshev'
    umap_components=2,     # usually 2 for plotting; can be >2 if you only use downstream clustering
    umap_randomize=False,  # True/False: if True, re-randomize init/seed behavior inside the implementation
    umap_clip_abs=0.0,     # 0.0 disables; otherwise clips feature magnitudes (implementation-dependent)

    # ----------------------------
    # t-SNE parameters
    # ----------------------------
    tsne_perplexity=10,    # Range 5-75   - opti 10
    tsne_distance='cosine',  #   Choices (typical): 'euclidean', 'cosine'
    tsne_dims=2,
    tsne_exaggeration=10,  #• Range 4-20   - opti 10
    tsne_learnrate=100,     #♦ Range 20-1000    - opti 100

    # ----------------------------
    # Scatter plot parameters
    # ----------------------------
    scatter_fig_size=(4, 4),
    scatter_point_size=2.5,
    scatter_point_alpha=0.7,
    scatter_edge_width=0.1,
    scatter_color_mode='enrichment',  #   Choices (as used in the plotting pipeline): 'enrichment', 'density', 'cluster'

    # ----------------------------
    # Density/enrichment params
    # ----------------------------
    density_nbins=150,
    density_sigma=4.0,
    density_use_log=True,
    density_min_rel=0.00,
    density_max_rel=1.00,
    density_cmap_name='plasma_r',  #   Examples: 'plasma', 'viridis', 'magma', 'inferno', 'cividis', 'turbo', 'gray', ...

    # density_enrichment_use_log:
    density_enrichment_eps=1e-12,
    density_enrichment_use_log=True,

    # ----------------------------
    # Gene-cluster localization heatmap along genome axis - as flattened heatmaps.
    #   This plot is generated by Plotting_Pipeline.py from the workbook sheet defined by PIPELINE['sheet_locus_tags'] (default: 'Locus Tags').
    # ----------------------------
    gchm_enable=True,
    gchm_colormap='plasma',  # any matplotlib colormap name (e.g. 'plasma', 'viridis').
    gchm_custom_cmaps_xlsx="",  # "" disables custom colormaps
    gchm_sigma=10,
    gchm_spread_factor=5,
    gchm_height_per_cluster=0.3,
    gchm_label_fontsize=10,
    gchm_dpi=300,
    gchm_cmap_min_rel=0.2,
    gchm_cmap_max_rel=1.0,
    gchm_output_filename='gene_cluster_heatmap_KS.png',
    gchm_show_fig=True,  # set False to avoid a second window popping up during the full pipeline run

    # ----------------------------
    # PCA parameters
    # ----------------------------
    pca_npcs=3,
    pca_center=True,
    pca_scale=True,

    # ----------------------------
    # Gene ordering / clustering along genome
    # ----------------------------
    gene_dist_metric='euclidean',   #   Choices (typical): 'euclidean', 'cosine', 'correlation', 'spearman'
    gene_linkage='single',          #   Choices: 'single', 'complete', 'average', 'ward'

    # ----------------------------
    # Cluster method parameters
    # ----------------------------
    #   Choices (typical): 'euclidean', 'cosine', 'manhattan'
    kmedoids_k=12,
    kmedoids_dist='euclidean',

    spectral_k=12,
    spectral_dist='euclidean',

    dbscan_eps=1.5,
    dbscan_minpts=10,
    dbscan_dist='euclidean',

    # ----------------------------
    # Feature normalization
    # ----------------------------
    center_features=True,
    scale_features=True,

    # ----------------------------
    # Feature clustering
    # ----------------------------
    feature_dist_metric='spearman',   #   Choices (typical): 'euclidean', 'cosine', 'correlation', 'spearman', 'pearson'
    feature_linkage='single',         #   Choices: 'single', 'complete', 'average', 'ward'

    # ----------------------------
    # Smoothing / binning along genome axis
    # ----------------------------
    apply_smoothing=True,
    smooth_window_genes=6,
    apply_binning=False,
    bin_size_genes=50,

    # ----------------------------
    # Heatmap aesthetics
    # ----------------------------
    heatmap_colormap_name='parula',
    heatmap_caxis_limits=(-0.5, 2.5),  # Tuple(min,max) or None to autoscale
    heatmap_fig_size=(18, 4),
    xtick_every_genes=500,

    # ----------------------------
    # Typography / export
    # ----------------------------
    font_name='Arial',
    font_size_axes=7,
    font_size_xticks=8,
    font_size_yticks=3,
    font_size_titles=10,
    colorbar_title_size=11,
    figure_dpi=300,

    # ----------------------------
    # Titles / labels
    # ----------------------------
    heatmap_title_template='Heatmap — {USAGE} | {CSET} | {DIMRED} | {CLUSTER}',
    scatter_title_template='{DIMRED} scatter — {USAGE} | {CSET}',
    colorbar_title_string=r'$\sigma_{\mathrm{codon}}$',

    # ----------------------------
    # Save formats
    # ----------------------------
    save_png=True,
    save_pdf=False,
    save_jpeg=False,

    # ----------------------------
    # Custom colormaps (optional)
    # ----------------------------
    use_custom_colormaps=False,
    custom_cmap_excel="",
    custom_cmap_sheets={},

    # ----------------------------
    # Export per-gene codon usage tables by cluster (1 sheet per cluster)
    #   Output workbook includes an "all genes" sheet named 'Genome locus tags'
    #   and then one sheet per cluster (cluster file columns).
    #
    #   export_cluster_codon_usage_mode:
    #       'ACU'      : absolute codon usage table (counts-like / original ACU input values)
    #       'RCU'      : relative codon usage among synonyms (NaN when AA is absent in gene)
    #       'RCU_devZ' : Z-scores of per-gene RCU vs genome baseline (baseline = ordered_genes)
    # ----------------------------
    export_cluster_codon_usage_enable=True,
    export_cluster_codon_usage_mode='RCU_devZ',   # 'ACU', 'RCU', 'RCU_devZ'
    export_cluster_codon_usage_codon_set='',      # '' -> reuse SET['codon_set']
    export_cluster_codon_usage_include_genome_sheet=True,
    export_cluster_codon_usage_genome_sheet_name='Genome locus tags',
    export_cluster_codon_usage_round_decimals=6,

    # ----------------------------
    # Text outputs
    # ----------------------------
    #   Master toggle for the 3 methods/legend txt files produced by write_all_texts()
    write_text_outputs=False,
)

PIPELINE = dict(
    # Sheet names inside the main pipeline Excel workbook
    sheet_coordinates="coordinates",
    sheet_locus_tags="Locus Tags",
    sheet_binary="Binary",

    # auto_run_plotting_pipeline:
    auto_run_plotting_pipeline=True,

    # plotting_pipeline_script_path:
    #   Path to Plotting_Pipeline.py (used by the plotting bridge)
    plotting_pipeline_script_path=os.path.join(HERE, "Plotting_Pipeline.py"),

    # per_cluster_suffix:
    #   Suffix for the "per cluster gene lists" workbook
    per_cluster_suffix="_PerClusterGeneLists.xlsx",
)

KS_SETTINGS = dict(
    alpha=0.01,       # alpha: significance threshold used downstream for interpretation/labeling
    method="binned",  # method: 'binned'
    bins=151,         # bins: integer number of bins for the 'binned' method
    n_perm=2000,      # n_perm: number of permutations if using a permutation-based p-value estimation
    random_seed=42,   # random_seed:  integer seed for reproducibility
)



def build_2dks_padj_matrix(ks_df: pd.DataFrame, cluster_order=None) -> pd.DataFrame:
    """Build an n×n matrix of BH-adjusted p-values (p_adj_BH) from the long-form 2D KS table.

    Preserves cluster order if cluster_order is provided (e.g. cluster_df.columns).
    """
    if ks_df is None or ks_df.empty:
        return pd.DataFrame()

    required = {"cluster_A", "cluster_B", "p_adj_BH"}
    missing = required - set(ks_df.columns)
    if missing:
        raise ValueError(f"2D KS raw table missing required columns: {sorted(missing)}")

    present = list(set(ks_df["cluster_A"].astype(str).tolist()) | set(ks_df["cluster_B"].astype(str).tolist()))
    present_set = set(present)

    if cluster_order is not None:
        ordered = [str(c) for c in list(cluster_order) if str(c) in present_set]
        extras = [c for c in present if c not in set(ordered)]
        names = ordered + extras
    else:
        names = present

    if not names:
        return pd.DataFrame()

    mat = pd.DataFrame(np.ones((len(names), len(names))), index=names, columns=names, dtype=float)

    for _, r in ks_df.iterrows():
        a = str(r["cluster_A"])
        b = str(r["cluster_B"])
        try:
            p = float(r["p_adj_BH"])
        except Exception:
            p = np.nan
        if a in mat.index and b in mat.columns:
            mat.loc[a, b] = p
            mat.loc[b, a] = p

    np.fill_diagonal(mat.values, 1.0)
    mat.index.name = "cluster"
    return mat


def build_neglog10_matrix(p_mat: pd.DataFrame, min_p: float = 1e-300) -> pd.DataFrame:
    """Return -log10(p) for a p-value matrix (clips finite values to [min_p, 1] to avoid inf)."""
    if p_mat is None or p_mat.empty:
        return pd.DataFrame()

    arr = p_mat.to_numpy(dtype=float)
    out = np.full_like(arr, np.nan, dtype=float)

    finite = np.isfinite(arr)
    clipped = np.clip(arr[finite], min_p, 1.0)
    out[finite] = -np.log10(clipped)

    neg = pd.DataFrame(out, index=p_mat.index, columns=p_mat.columns)
    neg.index.name = p_mat.index.name
    return neg


def main():
    clustering_results = run_codon_clustering(SET)

    codon_file = clustering_results["codon_file"]
    prefix_hint = clustering_results.get("prefix_hint", "")
    gene_symbol_map = clustering_results.get("gene_symbol_map", {})
    gene_desc_map = clustering_results.get("gene_desc_map", {})

    row_names = clustering_results["RowNames"]
    Y = clustering_results["Y"]
    ordered_genes = clustering_results["ordered_genes"]

    AA_df = clustering_results["AA_df"]
    C_abs_df = clustering_results["C_abs_df"]
    C_rel_df = clustering_results["C_rel_df"]
    features_reorder = clustering_results["features_reorder"]

    output_dir = clustering_results["output_dir"]
    output_subfolder = clustering_results["output_subfolder"]
    strain_prefix = clustering_results["strain_prefix"]

    base_folder = os.path.dirname(codon_file)

    cluster_file_path = choose_cluster_file(base_folder)
    cluster_df = read_locus_clusters(cluster_file_path)

    fasta_path = clustering_results.get("fasta_path")
    if fasta_path is None or (not os.path.isfile(str(fasta_path))):
        fasta_path = auto_find_fasta(base_folder, prefix_hint=prefix_hint)
        if fasta_path is None:
            fasta_path = choose_fasta(base_folder)

    locus_index, alias_map, id_map_df, missing_locus_headers, dup_counts = build_locus_index(fasta_path)

    row_names_can = np.array([canonicalize_id(t, alias_map) for t in row_names], dtype=object)
    ordered_genes_can = np.array([canonicalize_id(t, alias_map) for t in ordered_genes], dtype=object)

    row_names_can = enforce_unique_after_canon(row_names, row_names_can)
    ordered_genes_can = enforce_unique_after_canon(ordered_genes, ordered_genes_can)

    n_changed = int(np.sum([1 for a, b in zip(list(row_names), list(row_names_can)) if str(a) != str(b)]))
    print(f"[INFO] ID canonicalization: changed {n_changed}/{len(row_names_can)} codon-table IDs to canonical IDs.")

    cluster_df = canonicalize_cluster_df(cluster_df, alias_map)
    gene_symbol_map = canonicalize_generic_map(gene_symbol_map, alias_map)
    gene_desc_map = canonicalize_generic_map(gene_desc_map, alias_map)

    row_names = row_names_can
    ordered_genes = ordered_genes_can

    coords_df = build_coordinates_df(row_names, Y)
    locus_tags_df = build_locus_tags_sheet(ordered_genes, cluster_df, genome_colname="Genome locus tags")
    binary_df = build_binary_sheet(ordered_genes, cluster_df, gene_symbol_map, locus_index)

    metrics_df, missing_list = build_metrics_table(ordered_genes, locus_index, gene_symbol_map)
    summary_df = build_summary(metrics_df, total_requested=len(ordered_genes), missing_count=len(missing_list))
    quantitative_bis_df = build_quantitative_bis(metrics_df)

    meta_df = build_meta_table(
        SET, PIPELINE, clustering_results,
        fasta_path=fasta_path,
        locustags_path=cluster_file_path,
        summary_df=summary_df,
        output_subfolder=output_subfolder
    )

    ks_df = None
    if SET.get("do_2d_ks", True):
        ks_df = compute_2d_ks_for_clusters(Y=Y, row_names=row_names, cluster_txt_df=cluster_df, ks_settings=KS_SETTINGS)
    else:
        print("[INFO] 2D KS disabled by USER SETTINGS (SET['do_2d_ks']=False); skipping.")

    # Output base name: keep prefix, add ClusteringAnalysis
    out_base = os.path.join(output_dir, f"{strain_prefix}ClusteringAnalysis")
    out_xlsx = out_base + ".xlsx"

    with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
        AA_df.to_excel(writer, sheet_name="AA usage")
        C_abs_df.to_excel(writer, sheet_name="abs codon usage")
        C_rel_df.to_excel(writer, sheet_name="relative codon usage")

        coords_df.to_excel(writer, sheet_name=PIPELINE["sheet_coordinates"], index=False, header=False)

        feat_colname = "AA" if str(SET["usage_basis"]).upper() == "AA" else "Codon"
        pd.DataFrame({feat_colname: features_reorder}).to_excel(writer, sheet_name="Codons reordered", index=False)

        locus_tags_df.to_excel(writer, sheet_name=PIPELINE["sheet_locus_tags"], index=False)
        binary_df.to_excel(writer, sheet_name=PIPELINE["sheet_binary"], index=False)
        metrics_df.to_excel(writer, sheet_name="Quantitative", index=False)

        # Rename Quantitative_bis -> Quantitative Locus tags
        if quantitative_bis_df is None or quantitative_bis_df.empty:
            pd.DataFrame([{"info": "No binary (0/1) columns detected in 'Quantitative'."}]).to_excel(
                writer, sheet_name="Quantitative Locus tags", index=False
            )
        else:
            quantitative_bis_df.to_excel(writer, sheet_name="Quantitative Locus tags", index=False)

        meta_df.to_excel(writer, sheet_name="Meta", index=False)

        # 2D KS sheets
        ks_raw_sheet = "2D KS - raw data"
        ks_cmp_sheet = "2D KS - comparison"

        if ks_df is not None and (not ks_df.empty):
            ks_df.to_excel(writer, sheet_name=ks_raw_sheet, index=False)

            p_mat = build_2dks_padj_matrix(ks_df, cluster_order=cluster_df.columns)
            nlog_mat = build_neglog10_matrix(p_mat)

            header1 = "Matrix of BH-adjusted p-values (p_adj_BH) for all cluster-pair comparisons (diagonal = 1)."
            pd.DataFrame([[header1]]).to_excel(
                writer, sheet_name=ks_cmp_sheet, index=False, header=False, startrow=0, startcol=0
            )
            start_row_mat1 = 2
            p_mat.to_excel(writer, sheet_name=ks_cmp_sheet, startrow=start_row_mat1, startcol=0, index=True)

            header2 = "-log10(p_adj_BH) matrix (BH-adjusted p-values; diagonal = 0)."
            start_row_header2 = start_row_mat1 + len(p_mat.index) + 3
            pd.DataFrame([[header2]]).to_excel(
                writer, sheet_name=ks_cmp_sheet, index=False, header=False, startrow=start_row_header2, startcol=0
            )
            start_row_mat2 = start_row_header2 + 2
            nlog_mat.to_excel(writer, sheet_name=ks_cmp_sheet, startrow=start_row_mat2, startcol=0, index=True)

        else:
            msg = "2D KS disabled by user settings." if (not SET.get("do_2d_ks", True)) \
                  else "2D KS analysis could not be computed (see console for details)."
            pd.DataFrame([[msg]]).to_excel(writer, sheet_name=ks_raw_sheet, index=False, header=False)
            pd.DataFrame([[msg]]).to_excel(writer, sheet_name=ks_cmp_sheet, index=False, header=False)

        if missing_list:
            pd.DataFrame({"missing_locus_tag": missing_list}).to_excel(writer, sheet_name="Missing", index=False)

    print(f"\n[INFO] Pipeline Excel saved:\n  {out_xlsx}")

    out_cluster_xlsx = out_base + PIPELINE.get("per_cluster_suffix", "_PerClusterGeneLists.xlsx")
    write_per_cluster_workbook(
        out_path=out_cluster_xlsx,
        cluster_df=cluster_df,
        ordered_genes=ordered_genes,
        gene_symbol_map=gene_symbol_map,
        gene_desc_map=gene_desc_map,
        row_names=row_names,
        Y=Y,
        dimred_method=SET.get("dimred_method", "tsne"),
    )
    print(f"[INFO] Per-cluster Excel saved:\n  {out_cluster_xlsx}")

    # ---- NEW: Export per-gene codon usage tables by cluster (ACU / RCU / RCU_devZ) ----
    if bool(SET.get("export_cluster_codon_usage_enable", False)):
        try:
            mode = str(SET.get("export_cluster_codon_usage_mode", "RCU_devZ")).strip()
            mode_u = mode.strip().upper().replace(" ", "")
            # File naming consistent with the exporter convention
            if mode_u in ("RCUDEVZ", "RCU_DEVZ"):
                out_usage = out_base + "__RCUdevZ_by_cluster.xlsx"
            elif mode_u == "RCU":
                out_usage = out_base + "__RCU_by_cluster.xlsx"
            elif mode_u == "ACU":
                out_usage = out_base + "__ACU_by_cluster.xlsx"
            else:
                out_usage = out_base + "__RCU_by_cluster.xlsx"

            codon_set_export = str(SET.get("export_cluster_codon_usage_codon_set", "") or "").strip()
            if not codon_set_export:
                codon_set_export = str(SET.get("codon_set", "59_noSTOP_MW"))

            # Canonicalize index of C_abs_df to match canonical IDs used for clusters/order
            c_abs_idx_can = np.array([canonicalize_id(t, alias_map) for t in C_abs_df.index], dtype=object)
            c_abs_idx_can = enforce_unique_after_canon(C_abs_df.index, c_abs_idx_can)
            C_abs_df_can = C_abs_df.copy()
            C_abs_df_can.index = c_abs_idx_can

            out_usage = write_codon_usage_by_cluster_workbook(
                out_path=out_usage,
                ordered_genes=list(ordered_genes),
                cluster_df=cluster_df,
                c_abs_df=C_abs_df_can,
                mode=mode,
                codon_set=codon_set_export,
                include_genome_sheet=bool(SET.get("export_cluster_codon_usage_include_genome_sheet", True)),
                genome_sheet_name=str(SET.get("export_cluster_codon_usage_genome_sheet_name", "Genome locus tags")),
                round_decimals=SET.get("export_cluster_codon_usage_round_decimals", 6),
            )
            print(f"[INFO] Cluster codon-usage workbook saved:\n  {out_usage}")
        except Exception as e:
            print(f"[WARN] Could not export cluster codon-usage workbook: {e}")

    # ---- Write the 3 TXT files ----
    if bool(SET.get("write_text_outputs", True)):
        try:
            paths = write_all_texts(
                out_base=out_base,
                SET=SET,
                PIPELINE=PIPELINE,
                KS_SETTINGS=KS_SETTINGS,
                clustering_results=clustering_results,
                codon_file=codon_file,
                cluster_file_path=cluster_file_path,
                fasta_path=fasta_path,
                out_xlsx=out_xlsx,
                out_cluster_xlsx=out_cluster_xlsx,
            )
            print("[INFO] Text outputs saved:")
            for k, p in paths.items():
                print(f"  - {k}: {p}")
        except Exception as e:
            print(f"[WARN] Could not write text outputs: {e}")

    if PIPELINE.get("auto_run_plotting_pipeline", False):
        script_path = PIPELINE.get("plotting_pipeline_script_path", "")
        print(f"[INFO] Auto-run plotting enabled. Using script:\n  {script_path}")

        old_ob = os.environ.get("CODONPIPE_OUT_BASE")
        old_dm = os.environ.get("CODONPIPE_DIMRED_METHOD")

        # Gene-cluster heatmap env vars (save old values for restoration)
        _gchm_keys = [
            "CODONPIPE_GCHM_ENABLE",
            "CODONPIPE_GCHM_SHEET",
            "CODONPIPE_GCHM_COLORMAP",
            "CODONPIPE_GCHM_CUSTOM_CMAPS_XLSX",
            "CODONPIPE_GCHM_SIGMA",
            "CODONPIPE_GCHM_SPREAD_FACTOR",
            "CODONPIPE_GCHM_HEIGHT_PER_CLUSTER",
            "CODONPIPE_GCHM_LABEL_FONTSIZE",
            "CODONPIPE_GCHM_DPI",
            "CODONPIPE_GCHM_CMAP_MIN_REL",
            "CODONPIPE_GCHM_CMAP_MAX_REL",
            "CODONPIPE_GCHM_OUTPUT_FILENAME",
            "CODONPIPE_GCHM_SHOW",
        ]
        old_gchm = {k: os.environ.get(k) for k in _gchm_keys}

        os.environ["CODONPIPE_OUT_BASE"] = out_base
        os.environ["CODONPIPE_DIMRED_METHOD"] = str(SET.get("dimred_method", "tsne"))

        # Request the gene-cluster localization heatmap from Plotting_Pipeline.py
        os.environ["CODONPIPE_GCHM_ENABLE"] = "1" if bool(SET.get("gchm_enable", True)) else "0"
        os.environ["CODONPIPE_GCHM_SHEET"] = str(PIPELINE.get("sheet_locus_tags", "Locus Tags"))
        os.environ["CODONPIPE_GCHM_COLORMAP"] = str(SET.get("gchm_colormap", "plasma"))
        os.environ["CODONPIPE_GCHM_CUSTOM_CMAPS_XLSX"] = str(SET.get("gchm_custom_cmaps_xlsx", "") or "")
        os.environ["CODONPIPE_GCHM_SIGMA"] = str(SET.get("gchm_sigma", 10))
        os.environ["CODONPIPE_GCHM_SPREAD_FACTOR"] = str(SET.get("gchm_spread_factor", 5))
        os.environ["CODONPIPE_GCHM_HEIGHT_PER_CLUSTER"] = str(SET.get("gchm_height_per_cluster", 0.3))
        os.environ["CODONPIPE_GCHM_LABEL_FONTSIZE"] = str(SET.get("gchm_label_fontsize", 10))
        os.environ["CODONPIPE_GCHM_DPI"] = str(SET.get("gchm_dpi", 300))
        os.environ["CODONPIPE_GCHM_CMAP_MIN_REL"] = str(SET.get("gchm_cmap_min_rel", 0.2))
        os.environ["CODONPIPE_GCHM_CMAP_MAX_REL"] = str(SET.get("gchm_cmap_max_rel", 1.0))
        os.environ["CODONPIPE_GCHM_OUTPUT_FILENAME"] = str(SET.get("gchm_output_filename", "gene_cluster_heatmap_KS.png"))
        os.environ["CODONPIPE_GCHM_SHOW"] = "1" if bool(SET.get("gchm_show_fig", True)) else "0"

        try:
            run_density_plot_script(out_xlsx, PIPELINE)
        finally:
            if old_ob is None:
                os.environ.pop("CODONPIPE_OUT_BASE", None)
            else:
                os.environ["CODONPIPE_OUT_BASE"] = old_ob

            if old_dm is None:
                os.environ.pop("CODONPIPE_DIMRED_METHOD", None)
            else:
                os.environ["CODONPIPE_DIMRED_METHOD"] = old_dm

            # Restore gene-cluster heatmap env vars
            for k, v in old_gchm.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

    print("\n[DONE] Full codon-usage → clusters → quantitative metrics pipeline finished.")
    plt.show()


if __name__ == "__main__":
    main()
