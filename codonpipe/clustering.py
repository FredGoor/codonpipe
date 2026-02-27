"""codonpipe.clustering

Core computations for codon/AA usage embedding and clustering.
"""

# codonpipe/clustering.py
import os
import re
import math
from datetime import datetime
from tkinter import Tk
from tkinter.filedialog import askopenfilename

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list
from scipy.ndimage import gaussian_filter1d, gaussian_filter
from scipy.stats import rankdata

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.cluster import SpectralClustering, DBSCAN, KMeans

from matplotlib.colors import ListedColormap

from .colormaps import get_cmap_any

try:
    from sklearn_extra.cluster import KMedoids
except Exception as e:
    print("[WARN] sklearn-extra KMedoids unavailable; falling back to KMeans:", e)
    KMedoids = None

try:
    import umap
except ImportError:
    umap = None


# -----------------------------
# Small file-selection helpers
# -----------------------------
def choose_excel(initialdir, title):
    """Open a file dialog to select an Excel workbook.

    Falls back to a console prompt if a GUI is not available.
    """
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title=title,
            initialdir=initialdir or None,
            filetypes=[("Excel files", "*.xlsx")]
        )
        root.destroy()
    except Exception:
        path = input(f"{title} (.xlsx): ").strip().strip('"').strip("'")
    if not path:
        raise SystemExit("No file selected.")
    return path


def choose_cluster_file(initialdir):
    """Open a file dialog to select a locus-tag cluster file.

    Falls back to a console prompt if a GUI is not available.
    """
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title="Select locus-tag cluster file (TXT/CSV/TSV/Excel)",
            initialdir=initialdir or None,
            filetypes=[
                ("Cluster files", "*.txt *.csv *.tsv *.xlsx *.xls *.xlsm"),
                ("Text files", "*.txt *.csv *.tsv"),
                ("Excel files", "*.xlsx *.xls *.xlsm"),
                ("All files", "*.*"),
            ]
        )
        root.destroy()
    except Exception:
        path = input("Cluster file path: ").strip().strip('"').strip("'")
    if not path:
        raise SystemExit("No cluster file selected.")
    return path


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def _normalize_colname(c):
    return str(c).strip().lower().replace(" ", "").replace("-", "").replace(".", "").replace("_", "")


def infer_prefix_from_codon_basename(base_name):
    s = str(base_name)
    s = re.sub(r'(?i)codonusage', '', s)
    s = re.sub(r'(?i)codons?', '', s)
    s = re.sub(r'__+', '_', s)
    s = s.strip('_').strip()
    return s


def auto_find_geneids_excel(codon_file):
    folder = os.path.dirname(codon_file)
    base = os.path.splitext(os.path.basename(codon_file))[0]

    for pat in (r'(?i)codons', r'(?i)codon'):
        new_base = re.sub(pat, 'geneIDs', base)
        if new_base != base:
            candidate = os.path.join(folder, new_base + ".xlsx")
            if os.path.exists(candidate):
                return candidate

    for suffix in ("_geneIDs", "_geneids", " geneIDs", " geneids"):
        candidate = os.path.join(folder, base + suffix + ".xlsx")
        if os.path.exists(candidate):
            return candidate

    hint = infer_prefix_from_codon_basename(base).lower().strip('_')
    candidates = []
    for fname in os.listdir(folder):
        if not fname.lower().endswith(".xlsx"):
            continue
        if os.path.abspath(os.path.join(folder, fname)) == os.path.abspath(codon_file):
            continue
        lname = fname.lower()
        if ("geneid" in lname) or ("gene_ids" in lname) or ("geneids" in lname):
            score = 0
            if hint and hint in lname:
                score -= 10
            score -= len(os.path.commonprefix([hint, lname]))
            candidates.append((score, fname))

    if candidates:
        candidates.sort()
        return os.path.join(folder, candidates[0][1])

    return None


def _infer_locus_tag_series(df, cols_l):
    cols = list(df.columns)
    locus_col = None
    for i, c in enumerate(cols_l):
        key = _normalize_colname(c)
        if key in ("locustag", "locustags", "locus"):
            locus_col = cols[i]
            break
    if locus_col is not None:
        return df[locus_col].astype(str).map(_clean_text)

    try:
        idx = pd.Series(df.index.astype(str))
        frac_like = (idx.str.contains(r"[A-Za-z]", regex=True)).mean()
        if frac_like > 0.5 and not idx.str.fullmatch(r"\d+").all():
            return idx.map(_clean_text)
    except Exception:
        pass

    if len(cols) > 0:
        return df[cols[0]].astype(str).map(_clean_text)

    return pd.Series([], dtype=str)


def load_geneids_maps(gene_file):
    if not gene_file or (not os.path.exists(gene_file)):
        print("[WARN] GeneIDs Excel not found; GeneSymbol/Description maps will be empty.")
        return {}, {}

    try:
        df = pd.read_excel(gene_file, dtype=str)
    except Exception:
        df = pd.read_excel(gene_file)

    if df is None or df.empty:
        print("[WARN] GeneIDs Excel is empty; GeneSymbol/Description maps will be empty.")
        return {}, {}

    cols = list(df.columns)
    cols_l = [str(c).strip().lower() for c in cols]
    locus_series = _infer_locus_tag_series(df, cols_l)

    gene_col = None
    for i, c in enumerate(cols_l):
        if _normalize_colname(c) == "genesymbol":
            gene_col = cols[i]
            break

    desc_col = None
    preferred = ["description", "product", "annotation", "function", "geneproduct", "genedescription"]
    for pref in preferred:
        for i, c in enumerate(cols_l):
            key = _normalize_colname(c)
            if key == pref:
                desc_col = cols[i]
                break
        if desc_col is not None:
            break

    if desc_col is None:
        for i, c in enumerate(cols_l):
            key = _normalize_colname(c)
            if ("desc" in key) or ("product" in key) or ("annot" in key) or ("function" in key):
                desc_col = cols[i]
                break

    gene_series = df[gene_col].astype(str).map(_clean_text) if gene_col is not None else pd.Series([""] * len(df), dtype=str)
    desc_series = df[desc_col].astype(str).map(_clean_text) if desc_col is not None else pd.Series([""] * len(df), dtype=str)

    gene_map, desc_map = {}, {}
    for lt, gs, ds in zip(locus_series.tolist(), gene_series.tolist(), desc_series.tolist()):
        if not isinstance(lt, str):
            continue
        lt = lt.strip()
        if lt == "" or lt.lower() == "nan":
            continue
        gs = "" if (not isinstance(gs, str) or gs.lower() == "nan") else gs.strip()
        ds = "" if (not isinstance(ds, str) or ds.lower() == "nan") else ds.strip()

        if lt not in gene_map or ((not gene_map[lt]) and gs):
            gene_map[lt] = gs
        if lt not in desc_map or ((not desc_map[lt]) and ds):
            desc_map[lt] = ds

    return gene_map, desc_map


# -----------------------------
# Codon usage tables
# -----------------------------
def compute_AA_ACU_RCU(codon_df):
    var_names = codon_df.columns.tolist()
    row_names = codon_df.index.to_numpy()
    var_names = [v.replace('END', 'STOP') for v in var_names]
    codon_df.columns = var_names

    aa_names = np.array([v.split('_')[0] for v in var_names], dtype=object)
    unique_all = pd.unique(aa_names)
    unique20 = [a for a in unique_all if a != 'STOP']

    C_abs = codon_df.to_numpy(dtype=float)
    n_genes, n_codons = C_abs.shape

    AA = np.zeros((n_genes, len(unique20)), dtype=float)
    for i, aa in enumerate(unique20):
        cols = (aa_names == aa)
        AA[:, i] = C_abs[:, cols].sum(axis=1)
    AA_df = pd.DataFrame(AA, index=row_names, columns=unique20)

    C_rel = np.zeros_like(C_abs)
    for j in range(n_codons):
        aa_j = aa_names[j]
        if aa_j == 'STOP':
            continue
        idxAA = unique20.index(aa_j)
        denom = AA[:, idxAA]
        denom_safe = denom.copy()
        denom_safe[denom_safe == 0] = np.nan
        C_rel[:, j] = C_abs[:, j] / denom_safe
    C_rel = np.nan_to_num(C_rel)
    C_rel_df = pd.DataFrame(C_rel, index=row_names, columns=var_names)

    return C_abs, C_abs.shape[1], AA_df, C_rel_df, aa_names, unique20, row_names, var_names


def subset_usage(SET, usage_basis, C_abs, C_rel_df, AA_df, aa_names, var_names):
    n_codons_total = C_abs.shape[1]
    mask64 = np.ones(n_codons_total, dtype=bool)
    mask61 = (aa_names != 'STOP')
    mask59 = mask61 & (aa_names != 'Met') & (aa_names != 'Trp')

    codon_set = SET['codon_set'].lower()
    if codon_set == '64_all':
        feat_mask = mask64
    elif codon_set == '61_nostop':
        feat_mask = mask61
    elif codon_set == '59_nostop_mw':
        feat_mask = mask59
    else:
        raise ValueError(f"Unknown codon_set: {SET['codon_set']}")

    usage_basis = usage_basis.upper()
    if usage_basis == 'AA':
        Usage = AA_df.to_numpy(dtype=float)
        feature_labels = AA_df.columns.to_numpy()
        isAA = True
    elif usage_basis == 'ACU':
        Usage = C_abs
        feature_labels = np.array(var_names, dtype=object)
        isAA = False
    elif usage_basis == 'RCU':
        Usage = C_rel_df.to_numpy(dtype=float)
        feature_labels = np.array(var_names, dtype=object)
        isAA = False
    else:
        raise ValueError(f"Unknown usage_basis: {usage_basis}")

    if not isAA:
        Usage = Usage[:, feat_mask]
        feature_labels = feature_labels[feat_mask]

    return Usage, feature_labels, isAA


def normalize_values(SET, Usage):
    values = Usage.copy()
    if SET['center_features']:
        values -= values.mean(axis=0, keepdims=True)
    if SET['scale_features']:
        sd = values.std(axis=0, ddof=0, keepdims=True)
        sd[sd == 0] = 1.0
        values /= sd
    return values


# -----------------------------
# DimRed
# -----------------------------
def run_dimred(SET, values, feature_labels):
    method = SET['dimred_method'].lower()

    if method == 'umap':
        if umap is None:
            raise ImportError("umap-learn not installed but dimred_method='umap'.")
        X_umap = values.copy()
        clip = SET.get('umap_clip_abs', 0.0)
        if clip > 0:
            X_umap = np.clip(X_umap, -clip, clip)
        random_state = None if SET.get('umap_randomize', False) else 0
        reducer = umap.UMAP(
            n_neighbors=SET['umap_neighbors'],
            min_dist=SET['umap_min_dist'],
            n_components=SET['umap_components'],
            metric=SET['umap_metric'],
            random_state=random_state,
        )
        return reducer.fit_transform(X_umap)

    if method == 'tsne':
        tsne = TSNE(
            n_components=SET['tsne_dims'],
            perplexity=SET['tsne_perplexity'],
            metric=SET['tsne_distance'],
            learning_rate=SET['tsne_learnrate'],
            early_exaggeration=SET['tsne_exaggeration'],
            init='random',
            random_state=0,
        )
        return tsne.fit_transform(values)

    if method == 'pca':
        n_features = values.shape[1]
        k = max(2, min(SET['pca_npcs'], n_features))
        Xpca = values.copy()
        if SET['pca_center'] and not SET['center_features']:
            Xpca -= Xpca.mean(axis=0, keepdims=True)
        if SET['pca_scale'] and not SET['scale_features']:
            s = Xpca.std(axis=0, ddof=0, keepdims=True)
            s[s == 0] = 1.0
            Xpca /= s
        pca = PCA(n_components=k, svd_solver='auto', random_state=0)
        return pca.fit_transform(Xpca)

    if method == 'none':
        return None

    raise ValueError(f"Unknown dimred_method: {SET['dimred_method']}")


# -----------------------------
# Clustering
# -----------------------------
def _order_by_cluster_then_hier(Z, labels, dist_metric):
    unique_labels = np.unique(labels)
    k = len(unique_labels)
    centroids = np.zeros((k, Z.shape[1]), dtype=float)
    for i, c in enumerate(unique_labels):
        centroids[i, :] = Z[labels == c, :].mean(axis=0)
    seq = np.argsort(centroids[:, 0])

    order = []
    for i in seq:
        c = unique_labels[i]
        members = np.where(labels == c)[0]
        if len(members) > 2:
            Dloc = pdist(Z[members, :], metric=dist_metric)
            tree = linkage(Dloc, method='single')
            tree_opt = optimal_leaf_ordering(tree, Dloc)
            sub_order = members[leaves_list(tree_opt)]
        else:
            sub_order = members
        order.extend(sub_order.tolist())
    return np.array(order, dtype=int), labels


def cluster_genes(SET, Z):
    n = Z.shape[0]
    method = SET['cluster_method'].lower()

    if method == 'hierarchical':
        Dg = pdist(Z, metric=SET['gene_dist_metric'])
        tree = linkage(Dg, method=SET['gene_linkage'])
        tree_opt = optimal_leaf_ordering(tree, Dg)
        order = leaves_list(tree_opt)
        labels = np.zeros(n, dtype=int)
        return order, labels

    if method == 'kmedoids':
        k = max(2, min(SET['kmedoids_k'], n - 1))
        if KMedoids is not None:
            model = KMedoids(
                n_clusters=k,
                metric=SET['kmedoids_dist'],
                random_state=0,
                init='k-medoids++',
                max_iter=300
            )
            labels = model.fit_predict(Z)
        else:
            print("[WARN] KMedoids unavailable; falling back to KMeans for 'kmedoids'.")
            km = KMeans(n_clusters=k, random_state=0, n_init='auto')
            labels = km.fit_predict(Z)
        return _order_by_cluster_then_hier(Z, labels, SET['gene_dist_metric'])

    if method == 'spectral':
        k = max(2, min(SET['spectral_k'], n - 1))
        spec = SpectralClustering(n_clusters=k, affinity='nearest_neighbors', random_state=0)
        labels = spec.fit_predict(Z)
        return _order_by_cluster_then_hier(Z, labels, SET['gene_dist_metric'])

    if method == 'dbscan':
        db = DBSCAN(eps=SET['dbscan_eps'], min_samples=SET['dbscan_minpts'], metric=SET['dbscan_dist'])
        labels = db.fit_predict(Z)

        unique_labels = sorted(l for l in np.unique(labels) if l != -1)
        order = []
        for c in unique_labels:
            members = np.where(labels == c)[0]
            if len(members) > 2:
                Dloc = pdist(Z[members, :], metric=SET['gene_dist_metric'])
                tree = linkage(Dloc, method='single')
                tree_opt = optimal_leaf_ordering(tree, Dloc)
                sub_order = members[leaves_list(tree_opt)]
            else:
                sub_order = members
            order.extend(sub_order.tolist())

        noise = np.where(labels == -1)[0]
        if len(noise) > 0:
            idx = np.argsort(Z[noise, 0])
            order.extend(noise[idx].tolist())
        return np.array(order, dtype=int), labels

    raise ValueError(f"Unknown cluster_method: {SET['cluster_method']}")


# -----------------------------
# Feature reordering / smoothing
# -----------------------------
def reorder_features(SET, values, feature_labels):
    metric = SET['feature_dist_metric'].lower()
    if metric == 'spearman':
        ranks = np.apply_along_axis(rankdata, 1, values.T)
        Df = pdist(ranks, metric='correlation')
    else:
        Df = pdist(values.T, metric=metric)
    tree = linkage(Df, method=SET['feature_linkage'])
    tree_opt = optimal_leaf_ordering(tree, Df)
    order = leaves_list(tree_opt)
    return feature_labels[order], order


def smooth_and_bin(SET, V):
    if SET['apply_smoothing']:
        window = SET['smooth_window_genes']
        sigma = max(window / 3.0, 1.0)
        V = gaussian_filter1d(V, sigma=sigma, axis=1, mode='nearest')

    bin_size = max(1, int(round(SET['bin_size_genes'])))
    if SET['apply_binning']:
        ncols = V.shape[1]
        nbins = int(np.ceil(ncols / bin_size))
        Vb = np.zeros((V.shape[0], nbins), dtype=float)
        for i in range(nbins):
            L = i * bin_size
            R = min((i + 1) * bin_size, ncols)
            Vb[:, i] = V[:, L:R].mean(axis=1)
        V = Vb

    return V, bin_size


# -----------------------------
# Custom colormaps
# -----------------------------
def load_colormap_sheet(excel_path, sheet_name):
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    if df.empty:
        raise ValueError(f"Sheet '{sheet_name}' in '{excel_path}' is empty.")

    cols_lower = [str(c).strip().lower() for c in df.columns]

    def find_channel(ch):
        for i, name in enumerate(cols_lower):
            if name.startswith(ch):
                return df.columns[i]
        raise ValueError(f"No column starting with '{ch}' in sheet '{sheet_name}'.")

    r_col = find_channel('r')
    g_col = find_channel('g')
    b_col = find_channel('b')

    df = df[[r_col, g_col, b_col]].dropna()
    colors = df.to_numpy(dtype=float)
    if colors.max() > 1.0:
        colors = colors / 255.0

    return ListedColormap(colors, name=str(sheet_name))


def load_custom_colormaps(SET):
    cmap_dict = {}
    if not SET.get('use_custom_colormaps', False):
        return cmap_dict

    excel_path = SET.get('custom_cmap_excel', '')
    if not excel_path:
        print("Custom colormaps enabled but 'custom_cmap_excel' is empty; skipping.")
        return cmap_dict
    if not os.path.exists(excel_path):
        print(f"Custom colormap Excel not found: {excel_path}; skipping.")
        return cmap_dict

    sheet_map = SET.get('custom_cmap_sheets', {})
    loaded = []
    for name, sheet in sheet_map.items():
        try:
            cmap_dict[name] = load_colormap_sheet(excel_path, sheet)
            loaded.append(name)
        except Exception as e:
            print(f"Could not load colormap '{name}' from sheet '{sheet}': {e}")

    if loaded:
        print("[INFO] Loaded custom colormaps:", ", ".join(loaded))
    return cmap_dict


# -----------------------------
# Density colors (scatter)
# -----------------------------
def _compute_density_colors(Y, nbins=150, sigma=4.0, use_log=True,
                           min_rel=0.0, max_rel=1.0, cmap=None,
                           metric="density", enrichment_eps=1e-12,
                           enrichment_use_log=True):
    if cmap is None:
        cmap = plt.get_cmap("plasma")
    if Y is None or Y.shape[0] == 0:
        return None

    metric = (metric or "density").strip().lower()
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
    Hs_lin = gaussian_filter(H, sigma=sigma, mode="nearest")

    if not np.isfinite(Hs_lin).any() or Hs_lin.max() <= 0:
        return cmap(np.zeros_like(x))

    ix = np.searchsorted(xedges, x, side="right") - 1
    iy = np.searchsorted(yedges, y, side="right") - 1
    ix = np.clip(ix, 0, nbins - 1)
    iy = np.clip(iy, 0, nbins - 1)
    dens_lin = Hs_lin[ix, iy]

    if metric == "density":
        dens_score = np.log1p(dens_lin) if use_log else dens_lin
    elif metric == "enrichment":
        base = float(np.nanmean(dens_lin)) if np.isfinite(dens_lin).any() else 0.0
        base = max(base, 0.0)
        ratio = (dens_lin + enrichment_eps) / (base + enrichment_eps)
        dens_score = np.log1p(ratio) if enrichment_use_log else ratio
    else:
        raise ValueError("metric must be 'density' or 'enrichment'")

    if not np.isfinite(dens_score).any() or dens_score.max() <= 0:
        return cmap(np.zeros_like(x))

    dmax = float(np.nanmax(dens_score))
    min_rel = max(0.0, min(1.0, float(min_rel)))
    max_rel = max(min_rel, min(1.0, float(max_rel)))

    vmin = min_rel * dmax
    vmax = max(vmin + 1e-12, max_rel * dmax)

    dens_clipped = np.clip(dens_score, vmin, vmax)
    norm = (dens_clipped - vmin) / (vmax - vmin)
    norm = np.clip(norm, 0.0, 1.0)
    return cmap(norm)


# -----------------------------
# Plotting
# -----------------------------
def plot_heatmap(SET, V, features_reorder, n_genes_full, bin_size, title_txt, custom_cmaps):
    fig, ax = plt.subplots(figsize=SET['heatmap_fig_size'], dpi=SET['figure_dpi'])

    cmap_name = SET['heatmap_colormap_name']
    cmap = get_cmap_any(cmap_name, custom_maps=custom_cmaps, fallback='plasma')

    im = ax.imshow(V, aspect='auto', origin='lower', cmap=cmap)
    if SET['heatmap_caxis_limits'] is not None:
        im.set_clim(*SET['heatmap_caxis_limits'])
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.tick_params(direction='out')
    cbar.set_label(SET['colorbar_title_string'], fontsize=SET['colorbar_title_size'])

    ax.set_xlabel('Genes ordered by codon usage similarity', fontname=SET['font_name'])
    ax.set_ylabel('Codons', fontname=SET['font_name'])
    ax.tick_params(direction='out')
    ax.tick_params(axis='x', labelsize=SET['font_size_xticks'])
    ax.tick_params(axis='y', labelsize=SET['font_size_yticks'])

    step = max(1, SET['xtick_every_genes'])
    last_multiple = (n_genes_full - 1) // step * step
    ticks_genes = np.arange(0, last_multiple + 1, step)
    tick_pos = 1 + (ticks_genes // bin_size if SET['apply_binning'] else ticks_genes)
    ax.set_xticks(tick_pos - 1)
    ax.set_xticklabels([str(x) for x in ticks_genes])

    ax.set_yticks(np.arange(len(features_reorder)))
    ax.set_yticklabels(features_reorder, fontsize=SET['font_size_yticks'])

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.set_title(title_txt, fontsize=SET['font_size_titles'], fontname=SET['font_name'])
    return fig, ax


def plot_scatter(SET, Y, usage_basis, title_cset, custom_cmaps):
    if Y is None:
        return None, None

    fig, ax = plt.subplots(figsize=SET['scatter_fig_size'], dpi=SET['figure_dpi'])

    point_size = SET.get('scatter_point_size', 2.0)
    alpha = SET.get('scatter_point_alpha', 0.6)
    edge_width = SET.get('scatter_edge_width', 0.0)
    color_mode = SET.get('scatter_color_mode', 'plain').lower().strip()

    colors = None
    if color_mode in ("density", "enrichment"):
        density_cmap_name = SET.get('density_cmap_name', 'plasma')
        density_cmap = get_cmap_any(density_cmap_name, custom_maps=custom_cmaps, fallback='plasma')

        colors = _compute_density_colors(
            Y,
            nbins=SET.get('density_nbins', 150),
            sigma=SET.get('density_sigma', 4.0),
            use_log=SET.get('density_use_log', True),
            min_rel=SET.get('density_min_rel', 0.0),
            max_rel=SET.get('density_max_rel', 1.0),
            cmap=density_cmap,
            metric=color_mode,
            enrichment_eps=SET.get("density_enrichment_eps", 1e-12),
            enrichment_use_log=SET.get("density_enrichment_use_log", True),
        )

    ax.scatter(
        Y[:, 0], Y[:, 1],
        s=point_size,
        c=colors if colors is not None else None,
        alpha=alpha,
        linewidths=edge_width
    )

    method = SET['dimred_method'].lower()
    if method == 'tsne':
        xlab, ylab = 'tSNE1', 'tSNE2'
    elif method == 'umap':
        xlab, ylab = 'UMAP1', 'UMAP2'
    elif method == 'pca':
        xlab, ylab = 'PC1', 'PC2'
    else:
        xlab, ylab = 'Dim 1', 'Dim 2'

    ax.set_xlabel(xlab, fontname=SET['font_name'], fontsize=SET['font_size_axes'])
    ax.set_ylabel(ylab, fontname=SET['font_name'], fontsize=SET['font_size_axes'])
    ax.tick_params(labelsize=max(1, SET['font_size_axes'] - 1))

    tag_usage = usage_basis.upper()
    tag_dim = SET['dimred_method'].upper()
    tag_clu = SET['cluster_method'].upper()
    title_txt = SET['scatter_title_template'].format(
        USAGE=tag_usage, CSET=title_cset, DIMRED=tag_dim, CLUSTER=tag_clu
    )

    if method == 'tsne':
        param = (f"tSNE params: perplexity={SET['tsne_perplexity']}, "
                 f"distance={SET['tsne_distance']}, "
                 f"exaggeration={SET['tsne_exaggeration']}, "
                 f"learnrate={SET['tsne_learnrate']}")
    elif method == 'umap':
        param = (f"UMAP params: n_neighbors={SET['umap_neighbors']}, "
                 f"min_dist={SET['umap_min_dist']:.3f}, "
                 f"metric={SET['umap_metric']}, "
                 f"randomize={int(SET['umap_randomize'])}")
    elif method == 'pca':
        param = (f"PCA params: nPCs={Y.shape[1]}, "
                 f"center={int(SET['pca_center'])}, "
                 f"scale={int(SET['pca_scale'])}")
    else:
        param = ""

    extra = ""
    cm = SET.get("scatter_color_mode", "plain").lower().strip()
    if cm == "enrichment":
        extra = " | color=enrichment(vs mean density)"
    elif cm == "density":
        extra = " | color=density"

    ax.set_title(f"{param}\n{title_txt}{extra}" if param else f"{title_txt}{extra}",
                 fontsize=SET['font_size_titles'], fontname=SET['font_name'])
    return fig, ax


def make_run_output_dir(base_folder, SET):
    method = SET['dimred_method'].lower()
    if method == 'tsne':
        per = SET.get('tsne_perplexity', 30)
        ex = SET.get('tsne_exaggeration', 12)
        lr = SET.get('tsne_learnrate', 200)
        name = f"tsne per{per} ex{ex} lr{lr}"
    elif method == 'umap':
        nn = SET.get('umap_neighbors', 15)
        md = SET.get('umap_min_dist', 0.1)
        md_str = f"{md:.3g}".replace('.', 'p')
        name = f"umap nn{nn} md{md_str}"
    elif method == 'pca':
        pcs = SET.get('pca_npcs', 2)
        name = f"pca npc{pcs}"
    else:
        name = "nodimred"

    out_dir = os.path.join(base_folder, name)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir, name


# =========================================================
# BIG FUNCTION MOVED OUT: run_codon_clustering
# =========================================================
def run_codon_clustering(SET):
    if SET['usage_basis'].upper() == 'ACU' and SET['codon_set'] != '61_noSTOP':
        print("Warning: for ACU, recommended codon_set is '61_noSTOP'.")
    if SET['usage_basis'].upper() == 'RCU' and SET['codon_set'] != '59_noSTOP_MW':
        print("Warning: for RCU, recommended codon_set is '59_noSTOP_MW'.")

    print("---------------------- Open input data ----------------------")

    input_mode = str(SET.get("input_mode", "excel")).strip().lower()

    # ------------------------------------------------------------------
    # Input mode A: Excel codon usage table (+ GeneIDs Excel)
    # ------------------------------------------------------------------
    if input_mode in ("excel", "xlsx", "table"):
        codon_file = choose_excel(SET.get('default_root', ''), "Select codon usage Excel (per-gene codon usage)")

        base = os.path.splitext(os.path.basename(codon_file))[0]
        prefix_hint = infer_prefix_from_codon_basename(base)
        strain_prefix = prefix_hint
        if strain_prefix and not strain_prefix.endswith("_"):
            strain_prefix += "_"

        gene_file_auto = auto_find_geneids_excel(codon_file)
        if gene_file_auto and os.path.exists(gene_file_auto):
            gene_file = gene_file_auto
        else:
            if gene_file_auto:
                print(f"Auto geneID file not found: {os.path.basename(gene_file_auto)}")
            gene_file = choose_excel(os.path.dirname(codon_file),
                                     "Select GeneIDs Excel (must contain a 'GeneSymbol' column)")

        gene_symbol_map, gene_desc_map = load_geneids_maps(gene_file)

        base_folder = os.path.dirname(codon_file)
        output_dir, output_subfolder = make_run_output_dir(base_folder, SET)
        print(f"[INFO] Output subfolder: {output_subfolder}")

        C_abs_df = pd.read_excel(codon_file, index_col=0)
        fasta_path_for_input = None

    # ------------------------------------------------------------------
    # Input mode B: CDS FASTA (codon table + GeneIDs are computed on the fly)
    # ------------------------------------------------------------------
    elif input_mode in ("fasta", "cds_fasta", "ffn", "fna"):
        from .fasta_metrics import choose_fasta, compute_codon_usage_tables_from_cds_fasta

        fasta_path_for_input = str(SET.get("fasta_path", "") or SET.get("input_fasta", "")).strip()
        if not fasta_path_for_input:
            fasta_path_for_input = choose_fasta(SET.get('default_root', ''))
        if not os.path.isfile(fasta_path_for_input):
            raise FileNotFoundError(f"FASTA not found: {fasta_path_for_input}")

        base = os.path.splitext(os.path.basename(fasta_path_for_input))[0]
        prefix_hint = infer_prefix_from_codon_basename(base)
        strain_prefix = prefix_hint
        if strain_prefix and not strain_prefix.endswith("_"):
            strain_prefix += "_"

        codon_set = str(SET.get("codon_set", "61_noSTOP")).strip().lower()
        include_stops = (codon_set == "64_all")

        freq_df, geneids_df, abs_df, gc_pct = compute_codon_usage_tables_from_cds_fasta(
            fasta_path=fasta_path_for_input,
            row_id_mode=str(SET.get("fasta_row_id_mode", "primary")),
            trim_to_multiple_of_3=bool(SET.get("fasta_trim_to_multiple_of_3", True)),
            include_stops=bool(SET.get("fasta_include_stops", include_stops)),
            keep_first_duplicate=True,
        )

        table_mode = str(SET.get("fasta_codon_table_mode", "freq")).strip().lower()
        if table_mode in ("abs", "counts", "count"):
            C_abs_df = abs_df
        else:
            C_abs_df = freq_df

        gene_symbol_map = {}
        gene_desc_map = {}
        if geneids_df is not None and (not geneids_df.empty):
            for _, r in geneids_df.iterrows():
                lt = str(r.get("LocusTag", "")).strip()
                if not lt or lt.lower() == "nan":
                    continue
                gs = str(r.get("GeneSymbol", "") or "").strip()
                ds = str(r.get("ProteinDescription", "") or "").strip()
                if gs and gs.lower() != "nan" and gs != "NA":
                    gene_symbol_map[lt] = gs
                if ds and ds.lower() != "nan" and ds != "NA":
                    gene_desc_map[lt] = ds

        codon_file = fasta_path_for_input
        gene_file = None

        base_folder = os.path.dirname(fasta_path_for_input)
        output_dir, output_subfolder = make_run_output_dir(base_folder, SET)
        print(f"[INFO] Output subfolder: {output_subfolder}")

        if bool(SET.get("fasta_write_intermediate_excels", False)):
            try:
                codon_xlsx = os.path.join(output_dir, f"{base}_codons.xlsx")
                geneids_xlsx = os.path.join(output_dir, f"{base}_geneIDs.xlsx")
                C_abs_df.to_excel(codon_xlsx, index=True)
                if geneids_df is not None:
                    geneids_df.to_excel(geneids_xlsx, index=False)
                with open(os.path.join(output_dir, "GC_content.txt"), "w", encoding="utf-8") as fh:
                    fh.write(f"GC content of concatenated CDS: {gc_pct:.2f}%\n")
                print(f"[INFO] Wrote intermediate tables in: {output_dir}")
            except Exception as e:
                print(f"[WARN] Could not write intermediate FASTA-derived tables: {e}")

    else:
        raise ValueError("SET['input_mode'] must be 'excel' or 'fasta'.")
    C_abs, _, AA_df, C_rel_df, aa_names, _, RowNames, var_names = compute_AA_ACU_RCU(C_abs_df)

    Usage, feature_labels, _ = subset_usage(
        SET, SET['usage_basis'], C_abs, C_rel_df, AA_df, aa_names, var_names
    )
    values = normalize_values(SET, Usage)
    nGenes_full = values.shape[0]

    Y = run_dimred(SET, values, feature_labels)
    if Y is None:
        print("[WARN] dimred_method=='none'; falling back to PCA(2) for coordinates.")
        SET_fallback = dict(SET)
        SET_fallback['dimred_method'] = 'pca'
        Y = run_dimred(SET_fallback, values, feature_labels)

    gene_order, labels = cluster_genes(SET, Y)
    ordered_genes = RowNames[gene_order]

    features_reorder, feat_order = reorder_features(SET, values, feature_labels)
    V = values[gene_order, :].T[feat_order, :]
    V_smooth, bin_size = smooth_and_bin(SET, V)

    custom_cmaps = load_custom_colormaps(SET)

    tag_usage = SET['usage_basis'].upper()
    tag_dim = SET['dimred_method'].upper()
    tag_clu = SET['cluster_method'].upper()
    title_cset = f"{AA_df.shape[1]} AAs" if tag_usage == 'AA' else f"{Usage.shape[1]} codons"

    heatmap_title = SET['heatmap_title_template'].format(
        USAGE=tag_usage, CSET=title_cset, DIMRED=tag_dim, CLUSTER=tag_clu
    )

    fig_h, _ = plot_heatmap(SET, V_smooth, features_reorder, nGenes_full, bin_size, heatmap_title, custom_cmaps)
    fig_s, _ = plot_scatter(SET, Y, SET['usage_basis'], title_cset, custom_cmaps)

    tag_cset = title_cset.replace(' ', '')
    tag_pipe = f"{tag_dim}_{tag_clu}"
    if SET['dimred_method'].lower() == 'pca' and Y is not None:
        tag_pipe = f"{tag_pipe}_PC{Y.shape[1]}"
    out_base = os.path.join(output_dir, f"{base}_{tag_usage}_{tag_cset}_{tag_pipe}")

    if SET.get('save_png', True):
        fig_h.savefig(out_base + "_heatmap.png", dpi=SET['figure_dpi'], bbox_inches='tight')
        if fig_s is not None:
            fig_s.savefig(out_base + "_scatter.png", dpi=SET['figure_dpi'], bbox_inches='tight')
    if SET.get('save_pdf', False):
        fig_h.savefig(out_base + "_heatmap.pdf", bbox_inches='tight')
        if fig_s is not None:
            fig_s.savefig(out_base + "_scatter.pdf", bbox_inches='tight')
    if SET.get('save_jpeg', False):
        fig_h.savefig(out_base + "_heatmap.jpeg", dpi=SET['figure_dpi'], bbox_inches='tight')
        if fig_s is not None:
            fig_s.savefig(out_base + "_scatter.jpeg", dpi=SET['figure_dpi'], bbox_inches='tight')

    print("\n[INFO] Codon-usage clustering completed (heatmap + scatter).")

    return dict(
        input_mode=input_mode,
        fasta_path=fasta_path_for_input,
        codon_file=codon_file,
        gene_file=gene_file,
        strain_prefix=strain_prefix,
        prefix_hint=prefix_hint,
        gene_symbol_map=gene_symbol_map,
        gene_desc_map=gene_desc_map,
        RowNames=RowNames,
        Usage=Usage,
        values=values,
        Y=Y,
        ordered_genes=ordered_genes,
        features_reorder=features_reorder,
        AA_df=AA_df,
        C_abs_df=C_abs_df,
        C_rel_df=C_rel_df,
        output_dir=output_dir,
        output_subfolder=output_subfolder,
    )
