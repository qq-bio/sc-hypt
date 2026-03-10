#!/usr/bin/env python
# coding: utf-8

import os, sys, argparse
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path
from scipy import sparse
import re

import popv
from popv.preprocessing import Process_Query
from popv.annotation import annotate_data

popv.settings.compute_embedding = False
popv.settings.return_probabilities = False



################################################################################

parser = argparse.ArgumentParser()
parser.add_argument("-q", "--query", required=True, help="Path to query .h5ad")
parser.add_argument("-s", "--query-species", required=True, choices=["mouse", "rat"],
                    help="Species of query data")
parser.add_argument("-r", "--reference-mode", required=True,
                    help="One of: hypomap.c2, hypomap.c3, hca.heart, kpmp.sn_sc, kpmp.sn, hca.vascular")
parser.add_argument("-m", "--mode", default="retrain", choices=["fast", "inference", "retrain"],
                    help="PopV prediction mode")
parser.add_argument("--outdir", default="/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv",
                    help="Output directory")
parser.add_argument("--models-dir", default=None,
                    help="Directory for pretrained/trained models (defaults to <outdir>/models.<reference-mode>)")
parser.add_argument("--n-samples-per-label", type=int, default=800,
                    help="Cap per-label subsampling in reference")
args = parser.parse_args()


################################################################################

OUT_DIR = args.outdir
os.makedirs(OUT_DIR, exist_ok=True)


QUERY_H5AD = args.query
q_stem = Path(QUERY_H5AD).stem

CTMAP_FILE = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/reference_ct_to_major_ct.txt"
QUERY_MAJOR_COL = "subclass_level1"

REF_MODE_FILE = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/reference.setting.txt"
if not os.path.exists(REF_MODE_FILE):
    sys.exit(f"Missing reference setting file: {REF_MODE_FILE}")

ref_df = pd.read_table(REF_MODE_FILE)
required_cols = {"ref_mod", "file_name", "ref_labels_key", "human_col"}
missing = required_cols - set(ref_df.columns)
if missing:
    sys.exit(f"reference.setting.txt missing columns: {missing}")

row = ref_df.loc[ref_df["ref_mod"] == args.reference_mode]
if row.empty:
    choices = ", ".join(sorted(ref_df["ref_mod"].unique()))
    sys.exit(f"--reference-mode must be one of: {choices}")
row = row.iloc[0]

ref_folder = "/xdisk/mliang1/qqiu/reference/single_cell_ref/"
REF_H5AD        = os.path.join(ref_folder, row["file_name"])
REF_LABELS_KEY  = row["ref_labels_key"]       # e.g. C2_named / C3_named / ...
HUMAN_COL       = row["human_col"]   


if not os.path.exists(REF_H5AD):
    sys.exit(f"Reference file not found: {REF_H5AD}")


# Ortholog mapping
if args.query_species == "mouse":
    ORTHO_CSV = "/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2human.out.txt"
elif args.query_species == "rat":
    ORTHO_CSV = "/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.rat2human.out.txt"
elif args.query_species == "human":
    ORTHO_CSV = None  # no mapping needed
else:
    sys.exit("--query-species must be mouse, rat, or human")

if ORTHO_CSV and not os.path.exists(ORTHO_CSV):
    sys.exit(f"Ortholog CSV not found: {ORTHO_CSV}")


PRED_MODE = args.mode          # 'fast' | 'inference' | 'retrain'
HVG = 2000 if PRED_MODE == "retrain" else None


MODEL_DIR = args.models_dir or os.path.join(OUT_DIR, f"models.{args.reference_mode}")
os.makedirs(MODEL_DIR, exist_ok=True)


# PopV settings
REF_BATCH_KEY  = None
QUERY_BATCH_KEY = None

USE_ONTOLOGY = False 
CL_OBO_FOLDER = False 


METHODS = ["CELLTYPIST","Support_Vector","Random_Forest","KNN_BBKNN"]


# output
OUT_H5AD = os.path.join(OUT_DIR, f"{q_stem}.{args.reference_mode}.{PRED_MODE}.popv.h5ad")
OUT_META = os.path.join(OUT_DIR, f"{q_stem}.{args.reference_mode}.{PRED_MODE}.popv_labels.tsv")



print("\n[PopV config]")
print(f" query:            {QUERY_H5AD}")
print(f" query_species:    {args.query_species}")
print(f" reference_mode:   {args.reference_mode}")
print(f" reference_file:   {REF_H5AD}")
print(f" ref_labels_key:   {REF_LABELS_KEY}")
print(f" human_col:        {HUMAN_COL}")
print(f" ortholog_csv:     {ORTHO_CSV if ORTHO_CSV else '(none, human query)'}")
print(f" mode:             {PRED_MODE}")
print(f" HVG:              {HVG}")
print(f" n_samples/label:  {args.n_samples_per_label}")
print(f" methods:  {METHODS}")
print(f" models_dir:       {MODEL_DIR}")
print(f" out_h5ad:         {OUT_H5AD}")
print(f" out_labels_csv:   {OUT_META}\n")


################################################################################
def choose_hvg_by_size(n_cells: int) -> int:
    if n_cells < 500:
        return 1000
    elif n_cells < 1000:
        return 1500
    elif n_cells < 5000:
        return 2000
    elif n_cells < 10000:
        return 3000
    else:
        return 4000

def _norm_symbols(v):
    return pd.Index([str(g).strip() for g in v]).astype(str)

def load_adata(path):
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()
    adata.var_names = _norm_symbols(adata.var_names)
    return adata

def map_genes_to_target(adata, mapping_dict):
    import anndata as ad
    import numpy as np
    import pandas as pd

    def _norm(v): return pd.Index([str(x).strip() for x in v])

    v = _norm(adata.var_names)
    mask = v.isin(mapping_dict.keys())
    sub = adata[:, mask].copy()
    sub.var_names = _norm(pd.Index([mapping_dict[g] for g in sub.var_names]))

    # Collapse duplicates by summing
    if sub.var_names.has_duplicates:
        uniq = pd.unique(sub.var_names)
        if hasattr(sub.X, "tocsr"):
            X = sub.X.tocsr()
            from scipy.sparse import vstack
            cols = {g: np.where(sub.var_names == g)[0] for g in uniq}
            new_cols = []
            for g in uniq:
                Xg = X[:, cols[g]]
                new_cols.append(Xg.sum(axis=1))
            from scipy.sparse import hstack
            Xnew = hstack(new_cols)
            sub = ad.AnnData(Xnew, obs=sub.obs.copy(), var=pd.DataFrame(index=uniq))
        else:
            df = pd.DataFrame(sub.X, columns=sub.var_names, index=sub.obs_names)
            df = df.groupby(axis=1, level=0).sum()
            sub = ad.AnnData(df.values, obs=sub.obs.copy(), var=pd.DataFrame(index=df.columns))
    return sub

def load_biomart_orthologs(
    path,
    mouse_col="Gene name",
    human_col="Human gene name",
    confidence_col="Human orthology confidence [0 low, 1 high]",
    mouse_to_human_identity_col="Human % identity",
    human_to_mouse_identity_col="Mouse % identity",
    sep=None,  # auto-detect; set "\t" if you know it's TSV
    strategy="strict",  # "strict" | "best_identity" | "first"
):
    """
    Returns a dict mapping mouse_symbol -> human_symbol according to BioMart export.
    Filters to confidence==1, drops missing, resolves one-to-many per `strategy`.
    """

    # Auto-detect delimiter if not provided
    if sep is None:
        # Infer: if file has tabs in first line, assume TSV; else CSV
        with open(path, "r") as fh:
            head = fh.readline()
        sep = "\t" if "\t" in head else ","

    df = pd.read_csv(path, sep=sep, dtype=str)

    df.columns = [c.strip() for c in df.columns]

    required = {mouse_col, human_col, confidence_col}
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"Missing columns in BioMart file: {missing}")

    def _norm(s):
        return s.astype(str).str.strip()

    df[mouse_col] = _norm(df[mouse_col])
    df[human_col] = _norm(df[human_col])

    # Filter: confidence == 1
    conf = df[confidence_col].astype(str).str.strip()
    df = df[conf == "1"].copy()

    # Drop rows with empty symbols
    df = df[(df[mouse_col] != "") & (df[human_col] != "")]
    df = df.dropna(subset=[mouse_col, human_col])

    # Optional identity columns (for tie-breaks)
    for col in (mouse_to_human_identity_col, human_to_mouse_identity_col):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Remove exact duplicate pairs
    df = df.drop_duplicates(subset=[mouse_col, human_col])

    # Resolve one-to-many / many-to-one
    if strategy == "strict":
        # Keep only mouse symbols that map to exactly one human symbol
        counts = df.groupby(mouse_col)[human_col].nunique()
        keep_mouse = counts[counts == 1].index
        df = df[df[mouse_col].isin(keep_mouse)].copy()

        # Also ensure human side is unique (pure one-to-one)
        counts_h = df.groupby(human_col)[mouse_col].nunique()
        keep_human = counts_h[counts_h == 1].index
        df = df[df[human_col].isin(keep_human)].copy()

    elif strategy == "best_identity":
        # Prefer the highest % identity; fall back from human% to mouse% if needed
        score = None
        if mouse_to_human_identity_col in df.columns and df[mouse_to_human_identity_col].notna().any():
            score = mouse_to_human_identity_col
        elif human_to_mouse_identity_col in df.columns and df[human_to_mouse_identity_col].notna().any():
            score = human_to_mouse_identity_col

        if score is not None:
            df = (
                df.sort_values([mouse_col, score], ascending=[True, False])
                  .drop_duplicates(subset=[mouse_col], keep="first")
            )
        else:
            df = (
                df.sort_values([mouse_col, human_col])
                  .drop_duplicates(subset=[mouse_col], keep="first")
            )

        df = df.sort_values([human_col, mouse_col]).drop_duplicates(subset=[human_col], keep="first")

    elif strategy == "first":
        df = (
            df.sort_values([mouse_col, human_col])
              .drop_duplicates(subset=[mouse_col], keep="first")
        )
        df = df.sort_values([human_col, mouse_col]).drop_duplicates(subset=[human_col], keep="first")
    else:
        raise ValueError("strategy must be one of {'strict','best_identity','first'}")

    # Final mapping
    mapping = dict(zip(df[mouse_col], df[human_col]))

    # ---- Reporting ----
    print("Ortholog mapping summary:")
    print(f"  File: {Path(path).name}")
    print(f"  Strategy: {strategy}")
    print(f"  Rows after confidence==1 & non-empty: {len(df)}")
    print(f"  Unique mouse symbols mapped: {len(set(mapping.keys()))}")
    print(f"  Unique human symbols mapped: {len(set(mapping.values()))}")

    return mapping


def to_layer_compatible(arr, dtype=np.float32):
    # Accepts dense or any scipy.sparse; returns ndarray or CSR
    if sparse.issparse(arr):
        return arr.tocsr().astype(dtype)
    else:
        return np.asarray(arr, dtype=dtype)

def ensure_raw_counts_layer(adata, prefer_layers=("counts","raw","raw_counts")):
    # If already present but wrong format, fix it
    if "raw_counts" in adata.layers:
        adata.layers["raw_counts"] = to_layer_compatible(adata.layers["raw_counts"])
        return adata

    # Try common layer names
    for name in prefer_layers:
        if name in adata.layers:
            adata.layers["raw_counts"] = to_layer_compatible(adata.layers[name])
            return adata

    # Try .raw
    if adata.raw is not None:
        adata.layers["raw_counts"] = to_layer_compatible(adata.raw.X)
        return adata

    # Last resort (only if X is truly counts)
    adata.layers["raw_counts"] = to_layer_compatible(adata.X)
    return adata

def ensure_csr_X_and_counts(adata, counts_layer="raw_counts", use_counts_for_X=True, dtype=np.float32):
    # Ensure counts layer exists and is CSR/CSC/ndarray
    if counts_layer not in adata.layers:
        raise KeyError(f"'{counts_layer}' not found in adata.layers")

    Xc = adata.layers[counts_layer]
    if sparse.issparse(Xc):
        if not (sparse.isspmatrix_csr(Xc) or sparse.isspmatrix_csc(Xc)):
            Xc = Xc.tocsr()
        adata.layers[counts_layer] = Xc.astype(dtype)
    else:
        adata.layers[counts_layer] = np.asarray(Xc, dtype=dtype)

    # Ensure .X is also CSR/CSC/ndarray; optionally set it to counts
    X = adata.layers[counts_layer] if use_counts_for_X else adata.X
    if sparse.issparse(X):
        if not (sparse.isspmatrix_csr(X) or sparse.isspmatrix_csc(X)):
            X = X.tocsr()
        adata.X = X.astype(dtype)
    else:
        adata.X = np.asarray(X, dtype=dtype)
    return adata

def layer_int_check(A, layer="raw_counts", who="adata"):
    assert layer in A.layers, f"{who}: missing layer '{layer}'"
    X = A.layers[layer]
    is_ok_type = (sparse.isspmatrix_csr(X) or sparse.isspmatrix_csc(X) or isinstance(X, np.ndarray))
    print(f"{who}: type={type(X)}, ok_type={is_ok_type}")
    if sparse.issparse(X):
        data = X.data
    else:
        data = np.asarray(X).ravel()
    print(f"{who}: dtype={data.dtype}, min={data.min()}, max={data.max()}")
    nonneg = np.all(data >= 0)
    # “integer” means exactly whole numbers AND integer dtype (usually)
    whole = np.all(np.isfinite(data)) and np.all(np.equal(data, np.floor(data)))
    print(f"{who}: nonneg={nonneg}, whole_values={whole}")

def coerce_raw_counts_integer(A, layer="raw_counts"):
    X = A.layers[layer]
    if sparse.issparse(X):
        X = X.tocsr()
        data = X.data
    else:
        X = np.asarray(X)
        X = sparse.csr_matrix(X)
        data = X.data

    data = np.rint(data)
    data[data < 0] = 0

    data = data.astype(np.int64)
    X.data = data

    A.layers[layer] = X
    A.X = X.copy()
    return A


################################################################################

mapping = load_biomart_orthologs(
    ORTHO_CSV,
    mouse_col="Gene name",
    human_col=HUMAN_COL,
    confidence_col="Human orthology confidence [0 low, 1 high]",
    mouse_to_human_identity_col="Human % identity", 
    human_to_mouse_identity_col="Mouse % identity", 
    sep="\t", 
    strategy="strict"
)


q = load_adata(QUERY_H5AD)
r = load_adata(REF_H5AD)

q = ensure_raw_counts_layer(q, prefer_layers="raw_counts")
r = ensure_raw_counts_layer(r, prefer_layers="raw_counts")

q_mapped = map_genes_to_target(q, mapping)

common = q_mapped.var_names.intersection(r.var_names)
if len(common) < 500:
    raise ValueError(f"Too few shared genes after mapping: {len(common)}. Check mapping and species direction.")
q_mapped = q_mapped[:, common].copy()
r_sub = r[:, common].copy()

del q, r

q_mapped  = ensure_raw_counts_layer(q_mapped)
r_sub  = ensure_raw_counts_layer(r_sub)

q_mapped = coerce_raw_counts_integer(q_mapped, "raw_counts")
r_sub   = coerce_raw_counts_integer(r_sub,   "raw_counts")




# ---------------------------
# Celltype-to-celltype mapping
# ---------------------------
if QUERY_MAJOR_COL not in q_mapped.obs.columns:
    sys.exit(f"'{QUERY_MAJOR_COL}' not found in query obs.")

if not os.path.exists(CTMAP_FILE):
    sys.exit(f"Mapping table not found: {CTMAP_FILE}")

ctmap = pd.read_csv(CTMAP_FILE, sep="\t")
need_cols = {"ref_mod", "ref_labels", "query_labels"}
if not need_cols.issubset(ctmap.columns):
    sys.exit(f"Mapping table must have columns: {need_cols}")

annot_mask_global = pd.Series(False, index=q_mapped.obs_names)

# Iterate mapping rows
mask_rm = ctmap["ref_mod"].astype(str).str.strip().eq(args.reference_mode)
for q_label in sorted(ctmap.loc[mask_rm, "query_labels"].astype(str).str.strip().unique()):
    sub = ctmap[mask_rm & ctmap["query_labels"].astype(str).str.strip().eq(q_label)]
    if len(sub) <= 1:
        continue  # only proceed if >1 rows

    refs = []
    for s in sub["ref_labels"].astype(str):
        refs.extend([x.strip() for x in re.split(r"[;,]", s) if x.strip()])
    ref_list = sorted(set(refs))
    if len(ref_list) < 2:
        continue  

    if q_label not in set(q_mapped.obs[QUERY_MAJOR_COL].astype(str)):
        continue
    q_mask = q_mapped.obs[QUERY_MAJOR_COL].astype(str).eq(q_label)
    if q_mask.sum() == 0:
        continue
    
    if REF_LABELS_KEY not in r_sub.obs.columns:
        sys.exit(f"Reference missing labels column '{REF_LABELS_KEY}'")
    r_mask = r_sub.obs[REF_LABELS_KEY].astype(str).isin(ref_list)
    if r_mask.sum() < 50:
        print(f"[ct2ct] Skipping '{q_label}': only {r_mask.sum()} reference cells for list {ref_list}")
        continue

    q_ct = q_mapped[q_mask].copy() 
    r_ct = r_sub[r_mask].copy() 
    
    if PRED_MODE == "retrain":
        hvg_this = choose_hvg_by_size(q_ct.n_obs)
        hvg_this = min(hvg_this, q_ct.n_vars, r_ct.n_vars)
    else:
        hvg_this = None

    n_labels_present = r_ct.obs[REF_LABELS_KEY].astype(str).nunique()
    if n_labels_present < 2:
        only_lbl = r_ct.obs[REF_LABELS_KEY].astype(str).mode()[0]
        print(f"[ct2ct] '{q_label}': only one reference class present -> fallback assign '{only_lbl}'")
        continue

    pq_ct = Process_Query(
        query_adata = q_ct.copy(),
        ref_adata   = r_ct.copy(),
        ref_labels_key = REF_LABELS_KEY,
        ref_batch_key  = REF_BATCH_KEY,
        cl_obo_folder  = CL_OBO_FOLDER if USE_ONTOLOGY else False,
        query_batch_key = QUERY_BATCH_KEY,
        query_layer_key = "raw_counts",
        ref_layer_key   = "raw_counts",
        prediction_mode = PRED_MODE,
        unknown_celltype_label = "unknown",
        n_samples_per_label = args.n_samples_per_label,
        save_path_trained_models = MODEL_DIR,
        pretrained_scvi_path = None,
        relabel_reference_cells = False,
        hvg = hvg_this
    )
    
    A = pq_ct.adata
    lbl = REF_LABELS_KEY
    
    if "_ref_subsample" in A.obs.columns:
        ref_train_mask = A.obs["_ref_subsample"].astype(bool)
    else:
        ref_train_mask = pd.Series(False, index=A.obs_names)
    
    if ref_train_mask.any():
        n_classes_train = A.obs.loc[ref_train_mask, REF_LABELS_KEY].astype(str).nunique()
    else:
        n_classes_train = 0
    
    if n_classes_train < 2:
        print(f"[ct2ct] After Process_Query, only one training class -> fallback '{only_lbl}'")
        continue
    
    annotate_data(pq_ct.adata, methods=METHODS, save_path=OUT_DIR, methods_kwargs=None)

    obs_ct = pq_ct.adata.obs.copy()
    popv_cols = [c for c in obs_ct.columns if ("popv" in c.lower() or "pred" in c.lower() or "vote" in c.lower())]

    if len(popv_cols) == 0:
        print(f"[ct2ct] Warning: no PopV columns found to copy for '{q_label}'")
    else:
        target_idx = q_mapped.obs.index[q_mask]
    
        for c in popv_cols:
            if c not in q_mapped.obs.columns:
                q_mapped.obs[c] = pd.Series(pd.NA, index=q_mapped.obs.index, dtype="object")
    
            vals = obs_ct.reindex(target_idx)[c].values
            q_mapped.obs.loc[target_idx, c] = vals


# Write combined outputs
def _is_bad_obj_col(s: pd.Series) -> bool:
    if s.dtype != "object": return False
    types = {type(x) for x in s.dropna().tolist()}
    return not (len(types) <= 1 and (str in types or len(types) == 0))

bad_cols = []
for col in q_mapped.obs.columns:
    s = q_mapped.obs[col]
    if _is_bad_obj_col(s):
        bad_cols.append(col)
        # Try numeric coercion first
        ss = s.apply(lambda x: (np.array(x).ravel()[0] if isinstance(x,(list,tuple,np.ndarray)) else x))
        num = pd.to_numeric(ss, errors="ignore")
        if pd.api.types.is_numeric_dtype(num):
            q_mapped.obs[col] = num.astype("float32")
        else:
            # last resort: stringify
            q_mapped.obs[col] = s.astype(str)

if bad_cols:
    print(f"[sanitize] Fixed object columns in obs: {bad_cols}")

q_mapped.write(OUT_H5AD)
meta = q_mapped.obs.copy()
meta["cell"] = meta.index
meta.to_csv(OUT_META, sep="\t", index=False)



