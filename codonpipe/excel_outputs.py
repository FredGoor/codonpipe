"""codonpipe.excel_outputs

Excel I/O utilities for the codon-usage clustering pipeline.
"""

# codonpipe/excel_outputs.py
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd


def read_locus_clusters(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Cluster file not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".xlsx", ".xls", ".xlsm"):
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, sep=None, engine="python", dtype=str)

    for col in df.columns:
        df[col] = df[col].replace({np.nan: ""}).astype(str).map(lambda x: x.strip().strip('"').strip("'") if x != "nan" else "")
    return df


def clusters_from_df(df_clusters):
    clusters = {}
    for col in df_clusters.columns:
        vals = [v for v in df_clusters[col].tolist() if isinstance(v, str) and v != ""]
        clusters[col] = set(vals)
    return clusters


def build_locus_tags_sheet(ordered_tags, df_clusters, genome_colname="Genome locus tags"):
    ordered_tags = list(ordered_tags)
    n_genome = len(ordered_tags)
    n_rows_clusters = len(df_clusters)
    n_rows = max(n_genome, n_rows_clusters)

    data = {genome_colname: [""] * n_rows}
    for i, tag in enumerate(ordered_tags):
        data[genome_colname][i] = tag

    for col in df_clusters.columns:
        col_vals = [v if isinstance(v, str) else "" for v in df_clusters[col].tolist()]
        if len(col_vals) < n_rows:
            col_vals += [""] * (n_rows - len(col_vals))
        else:
            col_vals = col_vals[:n_rows]
        data[col] = col_vals

    return pd.DataFrame(data)


def build_binary_sheet(ordered_tags, df_clusters, gene_symbol_map, locus_index):
    ordered_tags = list(ordered_tags)
    clusters = clusters_from_df(df_clusters)

    gene_names = []
    missing_gene = 0
    for lt in ordered_tags:
        g = (gene_symbol_map.get(lt, "") or "") if gene_symbol_map else ""
        if not g:
            rec = locus_index.get(lt)
            if rec is not None:
                g = rec.get("gene_name", "") or ""
        if not g:
            missing_gene += 1
        gene_names.append(g)

    data = {"Genome locus tags": ordered_tags, "Gene name": gene_names}
    for cname, members in clusters.items():
        data[cname] = [1 if lt in members else 0 for lt in ordered_tags]

    if missing_gene:
        print(f"[INFO] Gene names missing for {missing_gene} of {len(ordered_tags)} locus tags (left blank in 'Binary').")

    return pd.DataFrame(data)


def build_coordinates_df(row_names, Y):
    if Y is None or Y.shape[0] == 0:
        raise ValueError("No DR coordinates available (Y is None or empty).")
    Y2 = np.column_stack([Y[:, 0], np.zeros_like(Y[:, 0])]) if Y.shape[1] < 2 else Y[:, :2]
    return pd.DataFrame(np.column_stack([np.asarray(row_names, dtype=object), Y2]))


def build_meta_table(SET, PIPELINE, clustering_results,
                     fasta_path, locustags_path,
                     summary_df, output_subfolder):
    # Meta sheet writer (short implementation; behavior unchanged).
    rows = []
    add = rows.append

    add({"Category": "Run", "Key": "Timestamp", "Value": datetime.now().isoformat(timespec="seconds")})
    add({"Category": "Run", "Key": "Output subfolder", "Value": output_subfolder})

    add({"Category": "Files", "Key": "Codon Excel", "Value": clustering_results["codon_file"]})
    add({"Category": "Files", "Key": "GeneIDs Excel", "Value": clustering_results["gene_file"]})
    add({"Category": "Files", "Key": "CDS FASTA", "Value": fasta_path})
    add({"Category": "Files", "Key": "Cluster file", "Value": locustags_path})

    add({"Category": "Dataset", "Key": "Strain prefix", "Value": clustering_results["strain_prefix"]})
    add({"Category": "Dataset", "Key": "N genes (codon table)", "Value": len(clustering_results["RowNames"])})
    add({"Category": "Dataset", "Key": "N genes (ordered)", "Value": len(clustering_results["ordered_genes"])})

    for k in ("usage_basis", "codon_set"):
        add({"Category": "Usage", "Key": k, "Value": SET.get(k)})

    dimred_keys = ["dimred_method", "tsne_perplexity", "tsne_exaggeration", "tsne_learnrate",
                   "umap_neighbors", "umap_min_dist", "umap_metric"]
    for k in dimred_keys:
        add({"Category": "DimRed", "Key": k, "Value": SET.get(k)})

    clu_keys = ["cluster_method", "kmedoids_k", "spectral_k", "dbscan_eps", "dbscan_minpts"]
    for k in clu_keys:
        add({"Category": "Clustering", "Key": k, "Value": SET.get(k)})

    for k in ("apply_smoothing", "smooth_window_genes"):
        add({"Category": "Smoothing", "Key": k, "Value": SET.get(k)})
    for k in ("apply_binning", "bin_size_genes"):
        add({"Category": "Binning", "Key": k, "Value": SET.get(k)})

    if summary_df is not None and not summary_df.empty:
        srow = summary_df.iloc[0]
        for col in summary_df.columns:
            add({"Category": "Quantitative summary", "Key": col, "Value": srow[col]})

    return pd.DataFrame(rows, columns=["Category", "Key", "Value"])


# =========================================================
# BIG FUNCTION MOVED OUT: write_per_cluster_workbook
# =========================================================
_INVALID_SHEET_CHARS = re.compile(r"[:\\/?*\[\]]")


def _safe_sheet_name(name: str, existing: set) -> str:
    base = str(name) if name is not None else "cluster"
    base = base.strip() or "cluster"
    base = _INVALID_SHEET_CHARS.sub("_", base)[:31]

    candidate = base
    if candidate not in existing:
        existing.add(candidate)
        return candidate

    i = 2
    while True:
        suffix = f"_{i}"
        max_len = 31 - len(suffix)
        candidate = (base[:max_len] if max_len > 0 else "") + suffix
        if candidate not in existing:
            existing.add(candidate)
            return candidate
        i += 1


def _coord_colnames(dimred_method: str):
    m = (dimred_method or "").lower().strip()
    if m == "umap":
        return "UMAP-X", "UMAP-Y"
    if m == "tsne":
        return "tSNE-X", "tSNE-Y"
    if m == "pca":
        return "PC1", "PC2"
    return "Dim1", "Dim2"


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def write_per_cluster_workbook(out_path,
                               cluster_df: pd.DataFrame,
                               ordered_genes: np.ndarray,
                               gene_symbol_map: dict,
                               gene_desc_map: dict,
                               row_names: np.ndarray,
                               Y: np.ndarray,
                               dimred_method: str):
    xcol, ycol = _coord_colnames(dimred_method)

    coords = {}
    if Y is not None and len(row_names) == Y.shape[0]:
        for tag, xy in zip(row_names.tolist(), Y[:, :2].tolist()):
            coords[str(tag)] = (float(xy[0]), float(xy[1]))

    order_index = {str(tag): i for i, tag in enumerate(ordered_genes.tolist())}

    existing = set()
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        # All Genes
        all_rows = []
        for lt in row_names.tolist():
            gs = (gene_symbol_map.get(lt, "") or "").strip() if gene_symbol_map else ""
            desc = (gene_desc_map.get(lt, "") or "").strip() if gene_desc_map else ""
            x, y = coords.get(str(lt), (np.nan, np.nan))
            all_rows.append({"Locus tags": lt, "GeneSymbol": gs, xcol: x, ycol: y, "Description": desc})
        pd.DataFrame(all_rows, columns=["Locus tags", "GeneSymbol", xcol, ycol, "Description"]).to_excel(
            writer, sheet_name=_safe_sheet_name("All Genes", existing), index=False
        )

        # Clusters
        for col in cluster_df.columns:
            raw = cluster_df[col].replace({np.nan: ""}).astype(str).tolist()
            members, seen = [], set()
            for v in raw:
                v = _clean_text(v)
                if not v or v.lower() == "nan" or v in seen:
                    continue
                seen.add(v)
                members.append(v)

            members.sort(key=lambda t: order_index.get(str(t), 10**12))
            rows = []
            for lt in members:
                gs = (gene_symbol_map.get(lt, "") or "").strip() if gene_symbol_map else ""
                desc = (gene_desc_map.get(lt, "") or "").strip() if gene_desc_map else ""
                x, y = coords.get(str(lt), (np.nan, np.nan))
                rows.append({"Locus tags": lt, "GeneSymbol": gs, xcol: x, ycol: y, "Description": desc})

            pd.DataFrame(rows, columns=["Locus tags", "GeneSymbol", xcol, ycol, "Description"]).to_excel(
                writer, sheet_name=_safe_sheet_name(col, existing), index=False
            )


# =========================================================
# NEW: Export per-gene codon usage tables by cluster (ACU / RCU / RCU_devZ)
# =========================================================

def _split_label_to_aa_codon(label: str):
    """
    Split a column label like "Ala_GCA" into ("Ala", "GCA").
    Returns ("", "") if parsing fails.
    """
    if label is None:
        return "", ""
    s = str(label).strip()
    if not s or "_" not in s:
        return "", ""
    aa, cod = s.split("_", 1)
    return aa.strip(), cod.strip().upper()


def _dedup_keep_order(values):
    seen = set()
    out = []
    for v in values:
        v = _clean_text(v) if isinstance(v, str) else v
        if not v or str(v).lower() == "nan":
            continue
        if v in seen:
            continue
        seen.add(v)
        out.append(v)
    return out


def _cluster_to_gene_list(cluster_df: pd.DataFrame) -> dict:
    """Return dict: cluster_name -> ordered list of locus tags (deduplicated, blanks removed)."""
    out = {}
    for col in cluster_df.columns:
        raw = cluster_df[col].replace({np.nan: ""}).astype(str).tolist()
        out[str(col)] = _dedup_keep_order(raw)
    return out


def _filter_codon_labels_by_set(labels, codon_set: str):
    """
    Apply codon_set filtering to codon labels.
    Supported conventions:
      - '59_noSTOP_MW' (exclude STOP, Met, Trp)
      - '61_noSTOP'    (exclude STOP)
      - '64_withSTOP'  (include all)
    Unknown codon_set => no filtering (safe default).
    """
    cs = (codon_set or "").lower().strip().replace("with_stop", "withstop")

    aa_list = []
    for lab in labels:
        aa, _ = _split_label_to_aa_codon(lab)
        aa_list.append(aa)

    aa_arr = np.array(aa_list, dtype=object)
    mask = np.ones(len(labels), dtype=bool)

    if cs in ("59_nostop_mw", "59_nostopmw", "59"):
        mask &= (aa_arr != "STOP") & (aa_arr != "Met") & (aa_arr != "Trp")
    elif cs in ("61_nostop", "61"):
        mask &= (aa_arr != "STOP")
    elif cs in ("64_withstop", "64", "all", "64_all"):
        mask &= np.ones(len(labels), dtype=bool)
    else:
        mask &= np.ones(len(labels), dtype=bool)

    return [lab for lab, keep in zip(labels, mask.tolist()) if keep]


def _ordered_codon_labels(labels):
    """
    Stable feature order: alphabetical by AA label, then by codon triplet.
    """
    groups = {}
    for lab in labels:
        aa, cod = _split_label_to_aa_codon(lab)
        if not aa or not cod:
            continue
        groups.setdefault(aa, []).append((cod, lab))

    ordered = []
    for aa in sorted(groups.keys()):
        groups[aa].sort(key=lambda t: t[0])
        ordered.extend([lab for _, lab in groups[aa]])
    return ordered


def compute_rcu_from_abs_df(c_abs_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-gene RCU from an absolute codon-usage table.

    For each AA group, codon values are divided by the AA total for that gene.
    If an AA total is 0, the codons for that AA are set to NaN (matching the exporter behavior).
    """
    df = c_abs_df.copy()
    labels = df.columns.tolist()

    aa_to_cols = {}
    for lab in labels:
        aa, _ = _split_label_to_aa_codon(lab)
        if aa:
            aa_to_cols.setdefault(aa, []).append(lab)

    out = pd.DataFrame(index=df.index, columns=labels, dtype=float)

    for aa, cols in aa_to_cols.items():
        denom = df[cols].sum(axis=1).replace(0, np.nan)
        out[cols] = df[cols].div(denom, axis=0)

    return out


def compute_rcu_devz(rcu_df: pd.DataFrame, baseline_genes: list) -> pd.DataFrame:
    """
    Compute RCU_devZ = (RCU_gene - mean_baseline) / std_baseline, codon-wise.

    - baseline_genes missing from the table are ignored
    - std==0 -> Z=0.0 by convention for genes where the codon value is defined
    - NaNs are ignored in mean/std.
    """
    base = [g for g in baseline_genes if g in rcu_df.index]
    if not base:
        raise ValueError("No baseline genes found in the provided RCU table; cannot compute RCU_devZ.")

    baseline = rcu_df.loc[base]
    mu = baseline.mean(axis=0, skipna=True)
    sd = baseline.std(axis=0, skipna=True, ddof=0)

    sd2 = sd.replace(0.0, np.nan)
    z = (rcu_df - mu) / sd2

    # Where sd was 0 but gene value exists, set to 0.0
    zero_var_cols = sd.index[sd == 0.0].tolist()
    if zero_var_cols:
        for c in zero_var_cols:
            mask = rcu_df[c].notna()
            z.loc[mask, c] = 0.0

    return z


def _build_cluster_codon_sheet_df(genes_in_cluster: list,
                                 values_df: pd.DataFrame,
                                 ordered_labels: list,
                                 round_decimals: int = 6) -> pd.DataFrame:
    """
    Build a sheet DataFrame with 3 header rows:
      Row 1: AA
      Row 2: Codon
      Row 3: Label
      Row 4+: gene rows (first column = locus tag)
    """
    aa_row = ["AA"] + [_split_label_to_aa_codon(lab)[0] for lab in ordered_labels]
    cod_row = ["Codon"] + [_split_label_to_aa_codon(lab)[1] for lab in ordered_labels]
    lab_row = ["Label"] + ordered_labels

    rows = [aa_row, cod_row, lab_row]

    for lt in genes_in_cluster:
        if lt in values_df.index:
            vals = values_df.loc[lt, ordered_labels].to_numpy(dtype=float)
        else:
            vals = np.full(len(ordered_labels), np.nan, dtype=float)

        if round_decimals is not None:
            vals2 = []
            for v in vals.tolist():
                if v is None or (isinstance(v, float) and np.isnan(v)):
                    vals2.append(np.nan)
                else:
                    try:
                        vals2.append(round(float(v), int(round_decimals)))
                    except Exception:
                        vals2.append(np.nan)
            vals = np.array(vals2, dtype=float)

        rows.append([lt] + vals.tolist())

    return pd.DataFrame(rows)


def write_codon_usage_by_cluster_workbook(out_path: str,
                                         ordered_genes: list,
                                         cluster_df: pd.DataFrame,
                                         c_abs_df: pd.DataFrame,
                                         mode: str = "RCU",
                                         codon_set: str = "59_noSTOP_MW",
                                         include_genome_sheet: bool = True,
                                         genome_sheet_name: str = "Genome locus tags",
                                         round_decimals: int = 6) -> str:
    """
    Write a workbook with 1 sheet per cluster, containing per-gene codon usage values.

    Parameters
    ----------
    out_path:
        Output .xlsx path.
    ordered_genes:
        List of "genome" locus tags (defines baseline gene set for RCU_devZ).
    cluster_df:
        DataFrame where each column is a cluster (cells = locus tags; blanks allowed).
    c_abs_df:
        Per-gene absolute codon usage table (rows=genes, columns like 'Ala_GCA').
    mode:
        'ACU' (absolute), 'RCU' (relative among synonyms with NaN when AA absent),
        or 'RCU_devZ' (Z vs genome baseline).
    codon_set:
        '59_noSTOP_MW', '61_noSTOP', '64_withSTOP', ...
    include_genome_sheet:
        If True, adds a first sheet containing ALL genes (ordered_genes).
    genome_sheet_name:
        Sheet name for the "All genes" sheet (default: 'Genome locus tags').
    round_decimals:
        Number of decimals for numeric output (None keeps full precision).

    Returns
    -------
    out_path (str)
    """
    mode_u = (mode or "RCU").strip().upper().replace(" ", "")
    if mode_u not in ("ACU", "RCU", "RCU_DEVZ", "RCUDEVZ"):
        raise ValueError("mode must be one of: 'ACU', 'RCU', 'RCU_devZ'.")

    # Keep only requested codon set and impose stable AA->codon ordering
    all_labels = c_abs_df.columns.tolist()
    labels_f = _filter_codon_labels_by_set(all_labels, codon_set)
    labels_o = _ordered_codon_labels(labels_f)

    # Compute values table
    if mode_u == "ACU":
        values_df = c_abs_df.copy()
    else:
        rcu_df = compute_rcu_from_abs_df(c_abs_df)
        if mode_u in ("RCU_DEVZ", "RCUDEVZ"):
            values_df = compute_rcu_devz(rcu_df, baseline_genes=list(ordered_genes))
        else:
            values_df = rcu_df

    cluster_to_genes = _cluster_to_gene_list(cluster_df)

    existing = set()
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        if include_genome_sheet:
            sheet_name = _safe_sheet_name(genome_sheet_name, existing)
            df_sheet = _build_cluster_codon_sheet_df(
                genes_in_cluster=list(ordered_genes),
                values_df=values_df,
                ordered_labels=labels_o,
                round_decimals=round_decimals
            )
            df_sheet.to_excel(writer, sheet_name=sheet_name, index=False, header=False)

        for cname, genes in cluster_to_genes.items():
            sheet_name = _safe_sheet_name(cname, existing)
            df_sheet = _build_cluster_codon_sheet_df(
                genes_in_cluster=genes,
                values_df=values_df,
                ordered_labels=labels_o,
                round_decimals=round_decimals
            )
            df_sheet.to_excel(writer, sheet_name=sheet_name, index=False, header=False)

    return out_path
