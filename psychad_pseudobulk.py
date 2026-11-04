#!/usr/bin/env python
"""
ARI GeroScreen — PsychAD lifespan DLPFC atlas: excitatory-neuron pseudobulk
===========================================================================
Sums raw counts of excitatory neurons per donor from the PsychAD lifespan
atlas (Yang et al., Nature 2026; CELLxGENE dataset 8c3ae5ad-...-7c5f84bb0664).
The per-donor matrix is then used by psychad_EN_age_signature.R to derive
age-associated genes in excitatory neurons.

Two ways to get the cells:
  1. CELLxGENE Census (default)     -- streams only the cells needed
       python psychad_pseudobulk.py --outdir PsychAD_EN
  2. Local .h5ad downloaded from CELLxGENE Discover (use if the dataset is
     not in your Census version, or if the Census connection is blocked)
       python psychad_pseudobulk.py --h5ad psychad_lifespan.h5ad --outdir PsychAD_EN

First run with --inspect to print the available cell types and columns
without aggregating anything:
       python psychad_pseudobulk.py --inspect

Outputs (in --outdir):
  psychad_EN_pseudobulk_counts.csv.gz  genes x donors, summed raw counts
  psychad_EN_donor_meta.csv            donor_id, age_years, sex, n_cells, ...
  psychad_celltype_counts.csv          cells per cell type (for checking)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

DATASET_ID = "8c3ae5ad-d4c1-4cbe-a8b8-7c5f84bb0664"
# CELLxGENE standard label for excitatory neurons; the regex also catches
# author labels such as "EN_L2_3_IT" if --celltype-col points to them
DEFAULT_EN_REGEX = r"glutamatergic|excitatory|^EN[_ ]"
OBS_COLS = ["donor_id", "development_stage", "sex", "disease", "cell_type",
            "is_primary_data"]
CHUNK = 50_000


# -----------------------------------------------------------------------------
def parse_age_years(stage):
    """'65-year-old human stage' -> 65; '9-month-old ...' -> 0.75; else NaN."""
    s = str(stage)
    m = re.search(r"(\d+)-year-old", s)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+)-month-old", s)
    if m:
        return float(m.group(1)) / 12
    if re.search(r"week|embryonic|fetal|newborn", s, re.I):
        return 0.0
    return np.nan


def select_en(obs, celltype_col, en_regex):
    if celltype_col not in obs.columns:
        sys.exit(f"Column '{celltype_col}' not in obs. Available: {list(obs.columns)}")
    keep = obs[celltype_col].astype(str).str.contains(en_regex, case=False, regex=True)
    if "disease" in obs.columns:
        keep &= obs["disease"].astype(str).str.lower().eq("normal")
    if keep.sum() == 0:
        sys.exit(f"No cells matched '{en_regex}' in '{celltype_col}'. "
                 "Run with --inspect and set --en-regex / --celltype-col.")
    return keep.values


def donor_table(obs_en):
    meta = (obs_en.groupby("donor_id", observed=True)
            .agg(n_cells=("donor_id", "size"),
                 development_stage=("development_stage", "first"),
                 sex=("sex", "first"))
            .reset_index())
    meta["age_years"] = meta["development_stage"].map(parse_age_years)
    return meta


def write_outputs(agg, donors, genes, meta, obs_all, celltype_col, outdir):
    os.makedirs(outdir, exist_ok=True)
    counts = pd.DataFrame(np.asarray(agg).T, index=genes, columns=donors)
    counts = counts.loc[counts.sum(axis=1) > 0]
    counts.round(0).astype(np.int64).to_csv(
        os.path.join(outdir, "psychad_EN_pseudobulk_counts.csv.gz"),
        compression="gzip")
    meta.set_index("donor_id").loc[donors].reset_index().to_csv(
        os.path.join(outdir, "psychad_EN_donor_meta.csv"), index=False)
    obs_all[celltype_col].value_counts().rename_axis("cell_type").rename("n_cells") \
        .to_csv(os.path.join(outdir, "psychad_celltype_counts.csv"))
    print(f"Wrote {counts.shape[0]} genes x {counts.shape[1]} donors to {outdir}")


# -----------------------------------------------------------------------------
# Route 1: CELLxGENE Census
# -----------------------------------------------------------------------------
def run_census(args):
    import cellxgene_census
    import tiledbsoma as soma

    with cellxgene_census.open_soma(census_version=args.census_version) as census:
        exp = census["census_data"]["homo_sapiens"]
        value_filter = f"dataset_id == '{DATASET_ID}'"
        obs_cols = ["soma_joinid"] + OBS_COLS
        if args.celltype_col not in obs_cols:
            obs_cols.append(args.celltype_col)
        with exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=value_filter)) as q:
            obs = q.obs(column_names=obs_cols).concat().to_pandas()
            if obs.empty:
                sys.exit("Dataset not found in this Census version. Download the .h5ad "
                         "from CELLxGENE Discover and rerun with --h5ad.")
            print(f"{len(obs):,} cells in dataset")
            print(obs[args.celltype_col].value_counts().head(40).to_string())
            if args.inspect:
                return

        keep = select_en(obs, args.celltype_col, args.en_regex)
        obs_en = obs.loc[keep].copy()
        meta = donor_table(obs_en)
        donors = meta["donor_id"].tolist()
        print(f"{len(obs_en):,} excitatory neurons from {len(donors)} donors")

        # Query only the excitatory neurons, by their Census row ids
        with exp.axis_query(
                "RNA",
                obs_query=soma.AxisQuery(coords=(obs_en["soma_joinid"].values,))) as q:
            var = q.var(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
            var_pos = pd.Series(np.arange(len(var)), index=var["soma_joinid"].values)
            donor_pos = pd.Series(pd.Categorical(obs_en["donor_id"], categories=donors).codes,
                                  index=obs_en["soma_joinid"].values)
            agg = np.zeros((len(donors), len(var)), dtype=np.float64)
            n_done = 0
            for tbl in q.X("raw").tables():
                rows = donor_pos.loc[tbl["soma_dim_0"].to_numpy()].to_numpy()
                cols = var_pos.loc[tbl["soma_dim_1"].to_numpy()].to_numpy()
                vals = tbl["soma_data"].to_numpy()
                agg += sp.coo_matrix((vals, (rows, cols)), shape=agg.shape).toarray()
                n_done += len(vals)
                print(f"  {n_done:,} nonzero entries aggregated", end="\r")
        print()
        write_outputs(agg, donors, var["feature_name"].values, meta, obs,
                      args.celltype_col, args.outdir)


# -----------------------------------------------------------------------------
# Route 2: local .h5ad (read in chunks; never loads the full matrix)
# -----------------------------------------------------------------------------
def run_h5ad(args):
    import h5py
    try:
        from anndata.io import read_elem, sparse_dataset          # anndata >= 0.11
    except ImportError:
        from anndata.experimental import read_elem, sparse_dataset  # anndata 0.10

    with h5py.File(args.h5ad, "r") as f:
        obs = read_elem(f["obs"])
        print(f"{len(obs):,} cells; obs columns: {list(obs.columns)}")
        print(obs[args.celltype_col].value_counts().head(40).to_string())
        if args.inspect:
            return

        # CELLxGENE files keep raw counts in raw/X (X is normalized)
        if "raw" in f and "X" in f["raw"]:
            X = sparse_dataset(f["raw/X"])
            var = read_elem(f["raw/var"])
            print("Using raw/X (raw counts)")
        else:
            X = sparse_dataset(f["X"])
            var = read_elem(f["var"])
            print("WARNING: no raw/X found; using X. Check it contains raw counts.")
        genes = var["feature_name"].astype(str).values if "feature_name" in var.columns \
            else var.index.astype(str).values

        keep = select_en(obs, args.celltype_col, args.en_regex)
        obs_en = obs.loc[keep]
        meta = donor_table(obs_en)
        donors = meta["donor_id"].tolist()
        codes = pd.Categorical(obs_en["donor_id"], categories=donors).codes
        idx = np.flatnonzero(keep)
        print(f"{len(idx):,} excitatory neurons from {len(donors)} donors")

        agg = np.zeros((len(donors), X.shape[1]), dtype=np.float64)
        for start in range(0, len(idx), CHUNK):
            sl = idx[start:start + CHUNK]
            Xc = X[sl]
            D = sp.csr_matrix((np.ones(len(sl)), (codes[start:start + CHUNK], np.arange(len(sl)))),
                              shape=(len(donors), len(sl)))
            agg += (D @ Xc).toarray()
            print(f"  {min(start + CHUNK, len(idx)):,} / {len(idx):,} cells", end="\r")
        print()
    write_outputs(agg, donors, genes, meta, obs, args.celltype_col, args.outdir)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", help="local .h5ad (otherwise uses CELLxGENE Census)")
    ap.add_argument("--census-version", default="stable")
    ap.add_argument("--celltype-col", default="cell_type")
    ap.add_argument("--en-regex", default=DEFAULT_EN_REGEX)
    ap.add_argument("--outdir", default="PsychAD_EN")
    ap.add_argument("--inspect", action="store_true",
                    help="print cell types and columns, then stop")
    args = ap.parse_args()
    run_h5ad(args) if args.h5ad else run_census(args)
