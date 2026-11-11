#!/usr/bin/env python3
"""
Fetch a subset of ABIDE Preprocessed (PCP) resting-state connectomes and save them in a
format the R pipeline (Codes/core/abide_pipeline.R) can read directly.

ABIDE PCP hosts precomputed ROI time series on S3 (no data-use agreement, no nilearn
dependency). For each selected subject we download the ROI time series for a given atlas,
form a Fisher-z transformed Pearson correlation connectome, zero the diagonal, and stack
all connectomes into a single CSV.

Outputs (under data/abide/):
  connectomes_full.csv  : (N*V) x V, no header. Subject i is rows [i*V : (i+1)*V] (a full
                          symmetric V x V connectome). Stacked-block layout keeps edge
                          ordering unambiguous for R.
  subjects.csv          : N rows aligned to the block order, with phenotype columns.
  meta.csv              : N, V, atlas, pipeline, strategy.

Usage:
  python3 Codes/fetch/abide_fetch.py --n 200 --atlas cc200
"""
import argparse, io, os, sys, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

S3 = "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data", "abide"))

# ROI count per atlas (for a sanity check on downloaded files)
ATLAS_ROIS = {"cc200": 200, "aal": 116, "ho": 111, "tt": 97, "ez": 116,
              "dosenbach160": 161, "cc400": 392}


def roi_url(file_id, pipeline, strategy, atlas):
    return f"{S3}/Outputs/{pipeline}/{strategy}/rois_{atlas}/{file_id}_rois_{atlas}.1D"


def fetch_one(file_id, pipeline, strategy, atlas):
    """Download one subject's ROI time series and return a Fisher-z connectome, or None.

    Uses curl via subprocess: macOS Python's urllib hits SSL cert-verification failures
    against this S3 bucket, while the system curl works reliably.
    """
    try:
        out = subprocess.run(
            ["curl", "-s", "--max-time", "90", roi_url(file_id, pipeline, strategy, atlas)],
            capture_output=True, check=True)
        raw = out.stdout
        if not raw or b"<Error>" in raw[:200] or b"NoSuchKey" in raw[:400]:
            return file_id, None
        ts = np.loadtxt(io.BytesIO(raw))            # timepoints x ROIs
        if ts.ndim != 2 or ts.shape[0] < 10:
            return file_id, None
        # drop ROIs with zero variance (occasionally an ROI is empty for a subject)
        good = ts.std(axis=0) > 0
        C = np.corrcoef(ts.T)
        C = np.nan_to_num(C, nan=0.0)
        C = np.clip(C, -0.999999, 0.999999)
        Z = np.arctanh(C)                            # Fisher z-transform (standard for FC)
        np.fill_diagonal(Z, 0.0)
        Z[~good, :] = 0.0
        Z[:, ~good] = 0.0
        return file_id, Z
    except Exception as e:
        return file_id, None


def main():
    global DATA
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="max subjects to download")
    ap.add_argument("--atlas", default="cc200", choices=list(ATLAS_ROIS))
    ap.add_argument("--pipeline", default="cpac")
    ap.add_argument("--strategy", default="filt_global")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--outdir", default=None, help="output directory (default data/abide)")
    args = ap.parse_args()
    if args.outdir:
        DATA = args.outdir

    os.makedirs(DATA, exist_ok=True)
    pheno_path = os.path.join(DATA, "Phenotypic_V1_0b_preprocessed1.csv")
    if not os.path.exists(pheno_path):
        subprocess.run(["curl", "-s", "--max-time", "120", "-o", pheno_path,
                        f"{S3}/Phenotypic_V1_0b_preprocessed1.csv"], check=True)
    pheno = pd.read_csv(pheno_path)
    pheno = pheno[pheno.FILE_ID != "no_filename"].copy()
    pheno = pheno.dropna(subset=["AGE_AT_SCAN"])
    # prefer the largest sites first so we get many subjects from few sites (less batch noise),
    # but still span several sites
    site_order = pheno.SITE_ID.value_counts().index.tolist()
    pheno["__site_rank"] = pheno.SITE_ID.map({s: i for i, s in enumerate(site_order)})
    pheno = pheno.sort_values(["__site_rank", "FILE_ID"]).head(args.n)
    file_ids = pheno.FILE_ID.tolist()
    print(f"Selected {len(file_ids)} subjects across {pheno.SITE_ID.nunique()} sites; "
          f"atlas={args.atlas} pipeline={args.pipeline} strategy={args.strategy}")

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_one, fid, args.pipeline, args.strategy, args.atlas): fid
                for fid in file_ids}
        done = 0
        for fut in as_completed(futs):
            fid, Z = fut.result()
            done += 1
            if Z is not None:
                results[fid] = Z
            if done % 25 == 0:
                print(f"  downloaded {done}/{len(file_ids)} ({len(results)} usable)")

    # keep phenotype order, drop subjects with no usable connectome / inconsistent V
    Vs = {Z.shape[0] for Z in results.values()}
    if not Vs:
        sys.exit("No usable connectomes downloaded.")
    V = max(Vs, key=lambda v: sum(Z.shape[0] == v for Z in results.values()))
    kept = [fid for fid in file_ids if fid in results and results[fid].shape[0] == V]
    print(f"Usable: {len(kept)} connectomes, V={V} nodes")

    blocks = [results[fid] for fid in kept]
    stacked = np.vstack(blocks)                      # (len(kept)*V) x V
    os.makedirs(DATA, exist_ok=True)
    np.savetxt(os.path.join(DATA, "connectomes_full.csv"), stacked, delimiter=",", fmt="%.6g")

    sub = pheno.set_index("FILE_ID").loc[kept].reset_index()
    cols = ["FILE_ID", "SITE_ID", "DX_GROUP", "AGE_AT_SCAN", "SEX", "FIQ", "VIQ", "PIQ"]
    sub[cols].to_csv(os.path.join(DATA, "subjects.csv"), index=False)

    pd.DataFrame([{"N": len(kept), "V": V, "atlas": args.atlas,
                   "pipeline": args.pipeline, "strategy": args.strategy}]).to_csv(
        os.path.join(DATA, "meta.csv"), index=False)

    print(f"Saved connectomes_full.csv ({stacked.shape}), subjects.csv ({len(kept)} rows), meta.csv")
    print(f"  age range: {sub.AGE_AT_SCAN.min():.1f}-{sub.AGE_AT_SCAN.max():.1f}, "
          f"mean {sub.AGE_AT_SCAN.mean():.1f}; "
          f"DX autism={int((sub.DX_GROUP==1).sum())} control={int((sub.DX_GROUP==2).sum())}")


if __name__ == "__main__":
    main()
