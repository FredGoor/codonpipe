"""codonpipe.fasta_metrics

FASTA parsing, identifier canonicalization, and per-gene sequence metrics.
"""

# codonpipe/fasta_metrics.py
import os
import re
import math
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from tkinter import Tk
from tkinter.filedialog import askopenfilename


# -----------------------------
# FASTA parsing helpers
# -----------------------------
_LOCUS_PATTERNS = [
    re.compile(r"\[locus_tag=([^\]]+)\]"),
    re.compile(r"\blocus_tag=([^\s\]]+)"),
    re.compile(r"^>+\s*([\w.\-:]+)"),
]
_GENE_PATTERNS = [
    re.compile(r"\[gene=([^\]]+)\]"),
    re.compile(r"\bgene=([^\s\]]+)"),
]
_PRODUCT_PATTERNS = [
    re.compile(r"\[(?:product|protein)=([^\]]+)\]"),
    re.compile(r"\b(?:product|protein)=([^\[]+?)(?:\s*\[|$)"),
]

_PRIMARY_ID_RE = re.compile(r"^>\s*([^\s]+)")
_LOCUS_TAG_RE = re.compile(r"\[locus_tag=([^\]]+)\]")
_PROTEIN_ID_RE = re.compile(r"\[protein_id=([^\]]+)\]")


def _clean_text(s):
    if not isinstance(s, str):
        return s
    return s.strip().strip('"').strip("'")


def parse_fasta_header(header):
    if not header or not header.startswith(">"):
        return None, ""

    locus = None
    for pat in _LOCUS_PATTERNS:
        m = pat.search(header)
        if m:
            locus = _clean_text(m.group(1))
            break

    gene = ""
    for pat in _GENE_PATTERNS:
        m = pat.search(header)
        if m:
            gene = _clean_text(m.group(1))
            break

    if not gene:
        for pat in _PRODUCT_PATTERNS:
            m = pat.search(header)
            if m:
                gene = _clean_text(m.group(1))
                break

    return locus, gene


def fasta_records(fasta_path):
    header = None
    seq_chunks = []
    with open(fasta_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line.rstrip("\n")
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def extract_strand_from_header(header):
    loc_m = re.search(r"\[location=([^\]]+)\]", header)
    if loc_m:
        loc = loc_m.group(1)
        return (0, "-") if "complement(" in loc else (1, "+")
    if re.search(r"strand\s*=\s*-\b", header) or re.search(r"\(-\)", header):
        return 0, "-"
    if re.search(r"strand\s*=\s*\+\b", header) or re.search(r"\(\+\)", header):
        return 1, "+"
    return 1, "+"


def extract_primary_id_from_header(header: str) -> str:
    if not header:
        return ""
    m = _PRIMARY_ID_RE.search(header)
    if m:
        return _clean_text(m.group(1))
    return _clean_text(header.lstrip(">").strip().split()[0]) if header.startswith(">") else _clean_text(str(header).split()[0])


def extract_protein_id_from_header(header: str) -> str:
    if not header:
        return ""
    m = _PROTEIN_ID_RE.search(header)
    return _clean_text(m.group(1)) if m else ""


# -----------------------------
# Canonicalization utilities
# -----------------------------
def canonicalize_id(x, alias_map: dict):
    if x is None:
        return ""
    s = str(x).strip()
    if not s or s.lower() == "nan":
        return ""

    m = _LOCUS_TAG_RE.search(s)
    if m:
        return _clean_text(m.group(1))

    if s.startswith(">"):
        s = s[1:].strip()
    if " " in s:
        s = s.split()[0].strip()

    return alias_map.get(s, s)


def canonicalize_generic_map(value_map: dict, alias_map: dict) -> dict:
    if not value_map:
        return {}
    out = {}
    for k, v in value_map.items():
        kk = canonicalize_id(k, alias_map)
        if not kk:
            continue
        vv = "" if v is None else str(v).strip()
        if kk not in out or (not out[kk] and vv):
            out[kk] = vv
    return out


def canonicalize_cluster_df(cluster_df: pd.DataFrame, alias_map: dict) -> pd.DataFrame:
    df = cluster_df.copy()
    for col in df.columns:
        df[col] = df[col].replace({np.nan: ""}).astype(str).map(
            lambda v: canonicalize_id(v, alias_map) if (v and str(v).lower() != "nan") else ""
        )
    return df


def enforce_unique_after_canon(orig_ids, canon_ids):
    seen = set()
    out = []
    for o, c in zip(list(orig_ids), list(canon_ids)):
        key = c if (c and c not in seen) else str(o)
        if key in seen or not key:
            base = key if key else str(o)
            i = 2
            key2 = base
            while key2 in seen or not key2:
                key2 = f"{base}__dup{i}"
                i += 1
            key = key2
        out.append(key)
        seen.add(key)
    return np.array(out, dtype=object)




# -----------------------------
# Codon table generation from CDS FASTA
# -----------------------------

# Standard genetic code (DNA alphabet)
_CODON2AA_1 = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"STOP","TAG":"STOP",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"STOP","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

_AA1_TO_AA3 = {
    "A":"Ala","R":"Arg","N":"Asn","D":"Asp","C":"Cys",
    "Q":"Gln","E":"Glu","G":"Gly","H":"His","I":"Ile",
    "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro",
    "S":"Ser","T":"Thr","W":"Trp","Y":"Tyr","V":"Val",
    "STOP":"STOP",
}

_STOP_CODONS = {"TAA", "TAG", "TGA"}

# Lexicographic codon order: AAA, AAC, ..., TTT
_BASE_ORDER = "ACGT"
_ALL_CODONS_64 = [a + b + c for a in _BASE_ORDER for b in _BASE_ORDER for c in _BASE_ORDER]

_RE_KV = re.compile(r"\[([A-Za-z0-9_\/\-\s]+)=([^\]]+)\]")
_RE_GENEID = re.compile(r"GeneID:(\d+)")
_RE_UNIPROT = re.compile(r"UniProtKB\/Swiss-Prot:([A-Z0-9]+)")


def _fasta_iter_simple(path: str):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def _parse_header_kv(header: str) -> dict:
    out = dict(
        primary_id="",
        locus_tag="NA",
        old_locus_tag="NA",
        gene="NA",
        product="NA",
        protein_id="NA",
        entrez="NA",
        uniprot="NA",
    )
    if not header:
        return out

    out["primary_id"] = header.split()[0].strip()

    kv = {}
    for m in _RE_KV.finditer(header):
        k = m.group(1).strip()
        v = m.group(2).strip()
        kv[k] = v

    if "locus_tag" in kv:
        out["locus_tag"] = kv["locus_tag"]
    if "old_locus_tag" in kv:
        out["old_locus_tag"] = kv["old_locus_tag"]
    if "gene" in kv:
        out["gene"] = kv["gene"]

    if "protein" in kv:
        out["product"] = kv["protein"]
    elif "product" in kv:
        out["product"] = kv["product"]

    if "protein_id" in kv:
        out["protein_id"] = kv["protein_id"]

    m = _RE_GENEID.search(header)
    if m:
        out["entrez"] = m.group(1)

    m = _RE_UNIPROT.search(header)
    if m:
        out["uniprot"] = m.group(1)

    # If old_locus_tag absent but primary looks like SL1344_0001, keep it
    if out["old_locus_tag"] == "NA" and re.match(r".+_\d+$", out["primary_id"]):
        out["old_locus_tag"] = out["primary_id"]

    return out


def _choose_row_id(meta: dict, row_id_mode: str) -> str:
    mode = str(row_id_mode or "primary").lower().strip()
    if mode == "primary":
        return meta.get("primary_id", "")
    if mode == "old":
        v = meta.get("old_locus_tag", "NA")
        return v if v != "NA" else meta.get("primary_id", "")
    if mode == "locus":
        v = meta.get("locus_tag", "NA")
        return v if v != "NA" else meta.get("primary_id", "")
    raise ValueError("row_id_mode must be one of: 'primary', 'old', 'locus'")


def _clean_seq_dna(seq: str) -> str:
    s = (seq or "").upper().replace("U", "T")
    s = re.sub(r"[^ACGT]", "N", s)
    return s


def _codon_counts(seq: str, codon_list, trim_to_multiple_of_3: bool = True):
    s = _clean_seq_dna(seq)
    if trim_to_multiple_of_3:
        s = s[:len(s) - (len(s) % 3)]
    counts = {c: 0 for c in codon_list}
    total_valid = 0
    for i in range(0, len(s), 3):
        c = s[i:i+3]
        if c in counts:
            counts[c] += 1
            total_valid += 1
    return counts, total_valid


def compute_codon_usage_tables_from_cds_fasta(
    fasta_path: str,
    row_id_mode: str = "primary",
    trim_to_multiple_of_3: bool = True,
    include_stops: bool = True,
    keep_first_duplicate: bool = True,
):
    """Build per-gene codon usage tables from a CDS FASTA.

    Returns
    -------
    codon_freq_df : pandas.DataFrame
        Per-gene codon frequencies (rows=genes, cols=AA3_CODON) in AAA..TTT order.
    geneids_df : pandas.DataFrame
        Parsed identifiers/annotations from headers.
    codon_abs_df : pandas.DataFrame
        Per-gene absolute codon counts (same columns/order as codon_freq_df).
    gc_percent : float
        GC% of concatenated unique CDS sequences.
    """
    if not fasta_path or (not os.path.isfile(fasta_path)):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    codon_list = list(_ALL_CODONS_64)
    if not include_stops:
        codon_list = [c for c in codon_list if c not in _STOP_CODONS]

    freq_rows = {}
    abs_rows = {}
    meta_rows = []
    seq_by_id = {}

    for header, seq in _fasta_iter_simple(fasta_path):
        meta = _parse_header_kv(header)
        row_id = _choose_row_id(meta, row_id_mode=row_id_mode)
        if not row_id:
            continue

        if keep_first_duplicate and row_id in freq_rows:
            continue

        counts, total_valid = _codon_counts(seq, codon_list=codon_list, trim_to_multiple_of_3=trim_to_multiple_of_3)
        if total_valid == 0:
            freqs = {c: 0.0 for c in codon_list}
        else:
            freqs = {c: counts[c] / float(total_valid) for c in codon_list}

        freq_rows[row_id] = freqs
        abs_rows[row_id] = counts
        seq_by_id[row_id] = _clean_seq_dna(seq)

        meta_rows.append(dict(
            LocusTag=row_id,
            GeneSymbol=meta.get("gene", "NA"),
            EntrezGeneID=meta.get("entrez", "NA"),
            ProteinDescription=meta.get("product", "NA"),
            RefSeqProteinID=meta.get("protein_id", "NA"),
            UniProtID=meta.get("uniprot", "NA"),
            RefSeq_LocusTag_RS=meta.get("locus_tag", "NA"),
            Old_LocusTag=meta.get("old_locus_tag", "NA"),
            PrimaryID=meta.get("primary_id", ""),
        ))

    if not freq_rows:
        raise ValueError("No CDS records could be parsed from the FASTA (empty file or invalid format).")

    freq_df = pd.DataFrame.from_dict(freq_rows, orient="index")[codon_list]
    abs_df  = pd.DataFrame.from_dict(abs_rows,  orient="index")[codon_list]

    new_cols = []
    for codon in codon_list:
        aa1 = _CODON2AA_1.get(codon, "X")
        aa3 = _AA1_TO_AA3.get(aa1, "Xxx")
        new_cols.append(f"{aa3}_{codon}")
    freq_df.columns = new_cols
    abs_df.columns = new_cols

    geneids_df = pd.DataFrame(meta_rows)

    concat = "".join(seq_by_id.values()).upper()
    gc = concat.count("G") + concat.count("C")
    gc_pct = 100.0 * gc / max(1, len(concat))

    return freq_df, geneids_df, abs_df, gc_pct


# -----------------------------
# Codon group metrics
# -----------------------------
CODON_GROUPS = {
    "sensitive": {
        "ACC", "ATC", "CCG", "ATT", "CTC", "GGC", "CTA", "CCA",
        "CGC", "CGA", "CGT", "GCC", "GTC", "CTT", "TCC", "CAA", "GTG", "GGT"
    },
    "insensitive": {
        "AGA", "AGG", "CGG", "GGA", "GGG", "ATA", "TTA", "TTG", "AGC", "AGT", "TCG", "ACG"
    },
    "MiaA": {
        "TGC", "TGT", "TTA", "TTG", "TTC", "TTT", "TCC", "TCT", "TCA", "TCG",
        "TGG", "TAC", "TAT"
    },
    "MnmEG": {
        "AGG", "AGA", "CAA", "CAG", "GAA", "GAG", "GGA", "GGG", "TTA", "TTG",
        "AAA", "AAG"
    },
    "Tgt": {
        "AAC", "AAT", "GAC", "GAT", "CAC", "CAT", "TAC", "TAT"
    },
    "regulatory": {
        "AGA", "AGG", "TTA", "TTG"
    },
    "regbis": {
        "AGA", "AGG", "TTA", "TTG", "GGA", "GGG"
    },
}

THRESHOLD_FLAGS = [
    ("perc_insensitive", 17, "bin_insensitive_gt17"),
    ("perc_MiaA",       20.0, "bin_MiaA_gt20"),
    ("perc_MnmEG",      25.0, "bin_MnmEG_gt25"),
    ("perc_Tgt",        17.0, "bin_Tgt_gt17"),
    ("perc_regulatory",  5.0, "bin_regulatory_gt5"),
    ("perc_regbis",      8.0, "bin_regbis_gt8"),
]


def gc_fraction(seq):
    if not seq:
        return math.nan
    s = seq.upper().replace("U", "T")
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    denom = a + t + g + c
    if denom == 0:
        return math.nan
    return (g + c) / denom


def gc_flags(gc_pct):
    if gc_pct is None or (isinstance(gc_pct, float) and math.isnan(gc_pct)):
        return {"gc_lt_48": 0, "gc_lt_45": 0, "gc_lt_40": 0, "gc_lt_35": 0}
    return {
        "gc_lt_48": 1 if gc_pct < 48.0 else 0,
        "gc_lt_45": 1 if gc_pct < 45.0 else 0,
        "gc_lt_40": 1 if gc_pct < 40.0 else 0,
        "gc_lt_35": 1 if gc_pct < 35.0 else 0,
    }


def _valid_codon(c):
    if len(c) != 3:
        return False
    for ch in c:
        if ch not in ("A", "C", "G", "T"):
            return False
    return True


def compute_codon_percentages(seq, groups=None):
    if groups is None:
        groups = CODON_GROUPS

    if not seq:
        return 0, {f"perc_{name}": math.nan for name in groups.keys()}

    s = seq.upper().replace("U", "T")
    total = 0
    counts = {name: 0 for name in groups.keys()}

    L = len(s) - (len(s) % 3)
    for i in range(0, L, 3):
        c = s[i:i+3]
        if not _valid_codon(c):
            continue
        total += 1
        for name, codon_set in groups.items():
            if c in codon_set:
                counts[name] += 1

    if total == 0:
        pcts = {f"perc_{name}": math.nan for name in groups.keys()}
    else:
        pcts = {f"perc_{name}": 100.0 * counts[name] / total for name in groups.keys()}

    return total, pcts


def perc_threshold_flags(pcts):
    flags = {}
    for key, cutoff, outname in THRESHOLD_FLAGS:
        val = pcts.get(key, math.nan)
        if val is None or (isinstance(val, float) and math.isnan(val)):
            flags[outname] = 0
        else:
            flags[outname] = 1 if (val > cutoff) else 0
    return flags


# =========================================================
# BIG FUNCTION MOVED OUT: build_locus_index
# =========================================================
def build_locus_index(fasta_path):
    index = {}
    dup_counts = {}
    missing_locus_headers = 0

    lt_re = re.compile(r"\[locus_tag=([^\]]+)\]")
    gb_re = re.compile(r"\[gbkey=([^\]]+)\]")

    counts = Counter()
    gb_by_lt = defaultdict(Counter)

    alias_map = {}
    idmap_rows = []

    print("[INFO] Parsing FASTA and indexing by locus_tag…")
    for header, seq in fasta_records(fasta_path):
        m = lt_re.search(header)
        if m:
            lt = m.group(1)
            counts[lt] += 1
            g = gb_re.search(header)
            gb = g.group(1) if g else "NA"
            gb_by_lt[lt][gb] += 1

        primary_id = extract_primary_id_from_header(header)
        protein_id = extract_protein_id_from_header(header)
        locus, gene_name = parse_fasta_header(header)

        canonical = locus or protein_id or primary_id
        if not canonical:
            missing_locus_headers += 1
            continue

        if primary_id:
            alias_map[primary_id] = canonical
        if protein_id:
            alias_map[protein_id] = canonical
        if locus:
            alias_map[locus] = locus

        idmap_rows.append({
            "canonical_id": canonical,
            "locus_tag": locus,
            "primary_id": primary_id,
            "protein_id": protein_id,
        })

        gc = gc_fraction(seq)
        gc_pct = None if (gc is None or (isinstance(gc, float) and math.isnan(gc))) else 100.0 * gc
        strand_bin, strand_sym = extract_strand_from_header(header)
        total_codons_valid, pcts = compute_codon_percentages(seq, CODON_GROUPS)

        rec = {
            "gene_name":            gene_name,
            "gc_fraction":          gc,
            "gc_percent":           gc_pct,
            "seq_length":           len(seq) if seq else 0,
            "total_codons_valid":   total_codons_valid,
            "strand_binary":        strand_bin,
            "strand_symbol":        strand_sym,
        }
        rec.update(gc_flags(gc_pct))
        rec.update(pcts)
        rec.update(perc_threshold_flags(pcts))

        if canonical in index:
            dup_counts[canonical] = dup_counts.get(canonical, 1) + 1
            continue
        index[canonical] = rec

    id_map_df = pd.DataFrame(idmap_rows).drop_duplicates()

    print(f"[INFO] Indexed {len(index)} unique entries from FASTA (canonical IDs).")
    if missing_locus_headers:
        print(f"[WARN] {missing_locus_headers} FASTA headers without detectable locus_tag/ID (ignored).")
    if dup_counts:
        print(f"[WARN] Detected {len(dup_counts)} canonical ID(s) with duplicate entries in FASTA (first kept).")

    print("Unique locus tags:", len(counts))
    print("How many locus tags appear twice:", sum(1 for v in counts.values() if v == 2))
    for lt, n in counts.most_common(10):
        if n > 1:
            print(lt, n, dict(gb_by_lt[lt]))

    return index, alias_map, id_map_df, missing_locus_headers, dup_counts


# -----------------------------
# FASTA selection helpers (kept)
# -----------------------------
def auto_find_fasta(folder, prefix_hint=""):
    exts = (".fa", ".fna", ".fasta", ".ffn", ".faa")
    candidates = []
    ph = (prefix_hint or "").lower().strip("_")

    for fname in os.listdir(folder):
        lower = fname.lower()
        if not lower.endswith(exts):
            continue

        score = 0
        if ph and ph in lower:
            score -= 10
        if "cds_from_genomic" in lower:
            score -= 6
        if "cds" in lower:
            score -= 3
        if lower.endswith(".ffn") or lower.endswith(".fna"):
            score -= 1
        if lower.endswith(".faa") or "protein" in lower:
            score += 5

        candidates.append((score, fname))

    if not candidates:
        return None

    candidates.sort()
    return os.path.join(folder, candidates[0][1])


def choose_fasta(initialdir):
    try:
        root = Tk()
        root.withdraw()
        path = askopenfilename(
            title="Select CDS FASTA file (e.g. *.fna / *.ffn)…",
            initialdir=initialdir,
            filetypes=[("FASTA files", "*.fa *.fna *.fasta *.ffn *.faa"),
                       ("All files", "*.*")]
        )
        root.destroy()
    except Exception:
        print("[NO GUI] Please provide path to CDS FASTA file:")
        path = input("> ").strip().strip('"').strip("'")

    if not path:
        raise SystemExit("No FASTA file selected; aborting.")
    return path


# =========================================================
# BIG FUNCTIONS MOVED OUT: metrics tables
# =========================================================
def build_metrics_table(ordered_tags, locus_index, gene_symbol_map):
    rows = []
    missing = []
    zero_flags = {name: 0 for _, _, name in THRESHOLD_FLAGS}

    for lt in ordered_tags:
        rec = locus_index.get(lt)
        gene_sym = (gene_symbol_map.get(lt, "") or "") if gene_symbol_map else ""

        if rec is None:
            missing.append(lt)
            row = {
                "locus_tag": lt,
                "gene_name": gene_sym,
                "gc_fraction": math.nan,
                "gc_percent": math.nan,
                "seq_length": math.nan,
                "total_codons_valid": 0,
                "strand_binary": math.nan,
                "strand_symbol": "",
                "gc_lt_48": 0,
                "gc_lt_45": 0,
                "gc_lt_40": 0,
                "gc_lt_35": 0,
                "perc_sensitive": math.nan,
                "perc_insensitive": math.nan,
                "perc_MiaA": math.nan,
                "perc_MnmEG": math.nan,
                "perc_Tgt": math.nan,
                "perc_regulatory": math.nan,
                "perc_regbis": math.nan,
            }
            row.update(zero_flags)
        else:
            row = {"locus_tag": lt, "gene_name": gene_sym if gene_sym else rec.get("gene_name", "")}
            row.update(rec)
            if gene_sym:
                row["gene_name"] = gene_sym

        rows.append(row)

    df = pd.DataFrame(rows)

    base_cols = [
        "locus_tag", "gene_name", "gc_fraction", "gc_percent", "seq_length",
        "total_codons_valid", "strand_binary", "strand_symbol",
        "gc_lt_48", "gc_lt_45", "gc_lt_40", "gc_lt_35",
        "perc_sensitive", "perc_insensitive", "perc_MiaA", "perc_MnmEG",
        "perc_Tgt", "perc_regulatory", "perc_regbis",
    ]
    flag_cols = [flag_name for _, _, flag_name in THRESHOLD_FLAGS]

    ordered_cols = [c for c in base_cols + flag_cols if c in df.columns]
    extra_cols = [c for c in df.columns if c not in ordered_cols]
    df = df[ordered_cols + extra_cols]

    return df, missing


def build_summary(metrics_df, total_requested, missing_count):
    matched_df = metrics_df[metrics_df["gc_percent"].notna()] if "gc_percent" in metrics_df.columns else metrics_df.iloc[0:0]

    summary = {
        "total_requested": [total_requested],
        "found_in_fasta": [total_requested - missing_count],
        "missing_in_fasta": [missing_count],
        "count_gc_lt_48": [int(matched_df.get("gc_lt_48", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_45": [int(matched_df.get("gc_lt_45", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_40": [int(matched_df.get("gc_lt_40", pd.Series([0])).sum())] if not matched_df.empty else [0],
        "count_gc_lt_35": [int(matched_df.get("gc_lt_35", pd.Series([0])).sum())] if not matched_df.empty else [0],
    }

    for _, _, flag_name in THRESHOLD_FLAGS:
        summary[f"count_{flag_name}"] = [int(metrics_df[flag_name].sum())] if flag_name in metrics_df.columns else [0]

    return pd.DataFrame(summary)


def build_quantitative_bis(metrics_df):
    if metrics_df is None or metrics_df.empty:
        return pd.DataFrame()
    if "locus_tag" not in metrics_df.columns:
        raise ValueError("metrics_df must contain a 'locus_tag' column.")

    locus_tags = metrics_df["locus_tag"].astype(str).tolist()
    out_cols = {}

    for col in metrics_df.columns:
        if col in ("locus_tag", "gene_name"):
            continue

        s_num = pd.to_numeric(metrics_df[col], errors="coerce")
        uniq = pd.unique(s_num.dropna())
        if len(uniq) == 0:
            continue
        if not np.all(np.isin(uniq, [0, 1])):
            continue

        hits = [lt for lt, v in zip(locus_tags, s_num.fillna(0).astype(int).tolist()) if v == 1]
        out_cols[str(col)] = hits

    if not out_cols:
        return pd.DataFrame()

    max_len = max(len(v) for v in out_cols.values())
    data = {}
    for k, v in out_cols.items():
        vv = list(v)
        if len(vv) < max_len:
            vv += [""] * (max_len - len(vv))
        data[k] = vv

    return pd.DataFrame(data)
