#!/usr/bin/env python3
"""
Fetch a subset of ACPI (Addiction Connectome Preprocessed Initiative) resting-state
connectomes from the Preprocessed Connectomes Project S3 bucket and save them in the SAME
format as Codes/fetch/abide_fetch.py, so Codes/core/abide_pipeline.R can read them directly
(load_abide(dir="data/acpi")).

ACPI is a second fully-open functional dataset (no data-use agreement) on the fcp-indi S3
bucket. Time series live at:
  data/Projects/ACPI/OutputsTimeSeries/{mta_1,nyu_1}/{ID}-session_1/ts_{atlas}_rois.csv
and the aggregated phenotype at:
  data/Projects/ACPI/PhenotypicData/acpi_aggregated_phenotypic_data.csv
(columns: ID, DATASET, SESSION, AGE_AT_SCAN_1, SEX, HANDEDNESS, BMI, ...).

NOTE on layout: ACPI ts files are ROIs x (1 metadata col + timepoints) -- the TRANSPOSE of
ABIDE's timepoints x ROIs. We drop the first column and transpose before correlating.

The cannabis/marijuana-use group label is NOT in the aggregated file; it is in the
site-specific files (mta_1_phenotypic_data.csv etc.). This script uses the aggregated file
(age/sex/BMI). To predict cannabis group, merge a site phenotype file on ID.

Usage:
  python3 Codes/fetch/acpi_fetch.py --n 150 --atlas cc200
"""
import argparse, io, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

S3 = "https://s3.amazonaws.com/fcp-indi/data/Projects/ACPI"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data", "acpi"))
SITES = ["mta_1", "nyu_1"]


def curl(url):
    out = subprocess.run(["curl", "-s", "--max-time", "120", url],
                         capture_output=True, check=True)
    return out.stdout


def list_subject_ids(site):
    """List subject session prefixes for a site via the S3 REST API."""
    url = (f"https://s3.amazonaws.com/fcp-indi?list-type=2"
           f"&prefix=data/Projects/ACPI/OutputsTimeSeries/{site}/&delimiter=/&max-keys=1000")
    xml = curl(url).decode("utf-8", "ignore")
    import re
    return re.findall(r"<Prefix>data/Projects/ACPI/OutputsTimeSeries/" + site + r"/([^<]+)/</Prefix>", xml)


def fetch_one(site, sess, atlas):
    url = f"{S3}/OutputsTimeSeries/{site}/{sess}/ts_{atlas}_rois.csv"
    try:
        raw = curl(url)
        if not raw or b"<Error>" in raw[:200]:
            return sess, None
        df = pd.read_csv(io.BytesIO(raw))
        # rows = ROIs, col 0 = node/volume metadata, remaining cols = timepoints -> transpose
        ts = df.iloc[:, 1:].to_numpy(dtype=float).T          # timepoints x ROIs
        if ts.shape[0] < 10:
            return sess, None
        good = ts.std(axis=0) > 0
        C = np.corrcoef(ts.T)
        C = np.nan_to_num(C, nan=0.0)
        C = np.clip(C, -0.999999, 0.999999)
        Z = np.arctanh(C)
        np.fill_diagonal(Z, 0.0)
        Z[~good, :] = 0.0; Z[:, ~good] = 0.0
        return sess, Z
    except Exception:
        return sess, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--atlas", default="cc200", choices=["cc200", "aal", "ho"])
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)

    pheno = pd.read_csv(io.BytesIO(curl(f"{S3}/PhenotypicData/acpi_aggregated_phenotypic_data.csv")))
    # phenotype ID is an int (e.g. 28046); folder IDs are zero-padded ("0028046").
    # Join on the integer value.
    pheno["ID"] = pheno["ID"].astype(int)

    # collect (site, session_dir, subject_id) tuples
    tasks = []
    for site in SITES:
        for sess in list_subject_ids(site):
            sid = sess.split("-")[0]
            tasks.append((site, sess, sid))
    tasks = tasks[: args.n]
    print(f"ACPI: attempting {len(tasks)} subject-sessions, atlas={args.atlas}")

    results = {}
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_one, site, sess, args.atlas): (sess, sid)
                for (site, sess, sid) in tasks}
        done = 0
        for fut in as_completed(futs):
            sess, Z = fut.result()
            sid = futs[fut][1]
            done += 1
            if Z is not None:
                results[sess] = (sid, Z)
            if done % 25 == 0:
                print(f"  {done}/{len(tasks)} ({len(results)} usable)")

    if not results:
        sys.exit("No usable ACPI connectomes downloaded.")
    Vs = [Z.shape[0] for _, (sid, Z) in [(k, v) for k, v in results.items()]]
    V = max(set(Vs), key=Vs.count)
    kept = [(sess, sid, Z) for sess, (sid, Z) in results.items() if Z.shape[0] == V]
    print(f"Usable: {len(kept)} connectomes, V={V}")

    stacked = np.vstack([Z for _, _, Z in kept])
    np.savetxt(os.path.join(DATA, "connectomes_full.csv"), stacked, delimiter=",", fmt="%.6g")

    ids = [sid for _, sid, _ in kept]
    sub = pheno.set_index("ID").reindex([int(s) for s in ids]).reset_index()
    # rename to the shared schema so abide_pipeline.R works unchanged for age/sex
    out = pd.DataFrame({
        "FILE_ID": ids,
        "SITE_ID": "ACPI",
        "AGE_AT_SCAN": sub["AGE_AT_SCAN_1"].values,
        "SEX": sub["SEX"].values,
        "BMI": sub.get("BMI", pd.Series([np.nan] * len(ids))).values,
    })
    out.to_csv(os.path.join(DATA, "subjects.csv"), index=False)
    pd.DataFrame([{"N": len(kept), "V": V, "atlas": args.atlas,
                   "pipeline": "acpi", "strategy": "OutputsTimeSeries"}]).to_csv(
        os.path.join(DATA, "meta.csv"), index=False)
    print(f"Saved to {DATA}: connectomes_full.csv {stacked.shape}, subjects.csv ({len(kept)}), meta.csv")
    print(f"  age range {out.AGE_AT_SCAN.min():.1f}-{out.AGE_AT_SCAN.max():.1f}, "
          f"mean {out.AGE_AT_SCAN.mean():.1f}")


if __name__ == "__main__":
    main()
