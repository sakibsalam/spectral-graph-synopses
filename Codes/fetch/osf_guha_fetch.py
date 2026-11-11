#!/usr/bin/env python3
"""
Fetch OSF "Power/Schaefer" functional-connectivity datasets into the shared
connectomes_full.csv / subjects.csv / meta.csv format.

DEFAULT DATASET (unchanged from before this file was generalized): OSF project 5zyws -- the
Guha et al. cognitive-control task-fMRI data the PI flagged (https://osf.io/5zyws), Power
229-node atlas, filenames "Power229_10_sub-<ID>_<cond>.txt", 1007 raw listing entries / 1003
unique files, 144 subjects, 7 conditions (rest, nbk_0back, nbk_1back, nbk_2back,
gng_initiation, gng_inhibition, tsw_shifting). Each file is a VxV Pearson correlation matrix.
API: https://api.osf.io/v2/nodes/<node>/files/osfstorage/  (download links in
data[i].links.download)

GENERALIZED (2026-09-11) to also fetch the two sibling parcellations of the SAME 144
subjects x 7 conditions, added via --node/--prefix (all other behavior, including the
explicit-page-number enumeration + de-dup-by-file-ID fix from 2026-09-10, is unchanged and
applies identically to every node):
  --node hnj7s --prefix Schaefer100_17   (component "Schaefer 100 parcel 17 networks", V=100)
  --node e4nph --prefix Schaefer200_17   (component "Schaefer 200 parcel 17 networks", V=200)
Every existing call site that does not pass --node/--prefix/--outdir keeps fetching exactly
the original 5zyws/Power229_10/data/osf_guha dataset, unchanged.

RESUMABLE: every matrix is cached as <outdir>/cache/<name>.npy and skipped if present, so the
download survives interruptions -- just re-run until "missing: 0". Each run (re)assembles
whatever is cached into <outdir>/connectomes_full.csv + subjects.csv (FILE_ID, SUBJECT, COND).

Usage:
  python3 Codes/fetch/osf_guha_fetch.py                                   # default: 5zyws / Power229_10 / data/osf_guha, all conditions
  python3 Codes/fetch/osf_guha_fetch.py --conds rest nbk_2back
  python3 Codes/fetch/osf_guha_fetch.py --node hnj7s --prefix Schaefer100_17 --outdir data/osf_schaefer100
  python3 Codes/fetch/osf_guha_fetch.py --node e4nph --prefix Schaefer200_17 --outdir data/osf_schaefer200
"""
import argparse, json, os, re, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

DEFAULT_NODE = "5zyws"
DEFAULT_PREFIX = "Power229_10"
HERE = os.path.dirname(os.path.abspath(__file__))
# HERE = <project_root>/Codes/fetch -- two levels up (fetch -> Codes -> project root) to reach
# the canonical <project_root>/data/... directory. (Previously only one level up was used for
# the default dataset, which resolved to a stray <project_root>/Codes/data/osf_guha shadow
# directory instead -- fixed 2026-09-10; preserved here for every dataset/outdir.)
DEFAULT_OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "data", "osf_guha"))


def curl_json(url, tries=8):
    for k in range(tries):
        out = subprocess.run(["curl", "-sL", "--max-time", "50", url], capture_output=True).stdout
        try:
            return json.loads(out)
        except Exception:
            time.sleep(2 + 2 * k)
    return None


def list_files(node):
    """Enumerate every file on the given OSF node by explicit page number (NOT by following the
    `links.next` cursor -- that under-enumerated by ~45 files for node 5zyws in an earlier
    version of this script, silently dropping entries across some page boundary). Walks page=
    1,2,3,... until a page comes back with zero entries, and de-duplicates by file ID (a handful
    of raw entries repeat the same file ID across pages on 5zyws; OSF's paging is not perfectly
    stable under concurrent edits to a node, so duplicates must be expected and collapsed rather
    than treated as new files -- the same de-dup logic is applied here for every node, including
    the two Schaefer siblings, without assuming their listing is clean).

    Returns dict: name -> (file_id, download_url), plus (raw_seen, n_unique_ids, n_unique_names)
    for the printed self-check. If two different file IDs ever map to the same name, the last one
    seen wins for download purposes and BOTH ids still count in the unique-id total.
    """
    by_id = {}          # file_id -> (name, url)
    raw_seen = 0
    page = 1
    while True:
        url = (f"https://api.osf.io/v2/nodes/{node}/files/osfstorage/"
               f"?page%5Bsize%5D=100&page={page}")
        d = curl_json(url)
        if d is None:
            print(f"  WARNING: page {page} failed to fetch after retries -- stopping enumeration"
                  f" early (results may be incomplete)")
            break
        data = d.get("data", [])
        if not data:
            break
        raw_seen += len(data)
        for f in data:
            by_id[f["id"]] = (f["attributes"]["name"], f["links"]["download"])
        page += 1
    names_to_id = {}
    for fid, (name, url) in by_id.items():
        names_to_id.setdefault(name, []).append(fid)
    files = {name: by_id[ids[-1]][1] for name, ids in names_to_id.items()}
    return files, raw_seen, len(by_id), len(names_to_id)


def fetch_one(name, url, cache_dir):
    cpath = os.path.join(cache_dir, name + ".npy")
    if os.path.exists(cpath):
        return name, True
    for k in range(4):
        out = subprocess.run(["curl", "-sL", "--max-time", "60", url], capture_output=True).stdout
        try:
            m = np.loadtxt(out.decode("utf-8", "ignore").splitlines())
            if m.shape[0] == m.shape[1] and m.shape[0] > 50:
                m = (m + m.T) / 2.0
                m = np.clip(m, -0.999999, 0.999999)
                z = np.arctanh(m); np.fill_diagonal(z, 0.0)
                np.save(cpath, z.astype(np.float32))
                return name, True
        except Exception:
            time.sleep(1.5 * (k + 1))
    return name, False


def assemble(outdir, cache_dir, pat, atlas, pipeline, conds_filter=None):
    """Build connectomes_full.csv / subjects.csv from everything in the cache."""
    cached = sorted(f[:-4] for f in os.listdir(cache_dir) if f.endswith(".npy"))
    rows = []
    for name in cached:
        m = pat.match(name)
        if not m:
            continue
        if conds_filter and m.group(2) not in conds_filter:
            continue
        rows.append((name, m.group(1), m.group(2)))
    if not rows:
        return 0
    V = int(np.load(os.path.join(cache_dir, rows[0][0] + ".npy")).shape[0])
    rows = [r for r in rows if np.load(os.path.join(cache_dir, r[0] + ".npy")).shape[0] == V]
    stacked = np.vstack([np.load(os.path.join(cache_dir, n + ".npy")) for (n, s, c) in rows])
    np.savetxt(os.path.join(outdir, "connectomes_full.csv"), stacked, delimiter=",", fmt="%.6g")
    import csv
    with open(os.path.join(outdir, "subjects.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["FILE_ID", "SUBJECT", "COND"])
        for (n, s, c) in rows:
            w.writerow([n, s, c])
    with open(os.path.join(outdir, "meta.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["N", "V", "atlas", "pipeline", "strategy"])
        w.writerow([len(rows), V, atlas, pipeline, "all"])
    import collections
    print(f"Assembled {len(rows)} connectomes (V={V}); conditions:",
          dict(collections.Counter(c for _, _, c in rows)))
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--node", default=DEFAULT_NODE,
                    help=f"OSF node id to fetch from (default: {DEFAULT_NODE}, the Power-229 dataset)")
    ap.add_argument("--prefix", default=DEFAULT_PREFIX,
                    help=f"filename prefix before '_sub-<id>_<cond>.txt' (default: {DEFAULT_PREFIX})")
    ap.add_argument("--outdir", default=None,
                    help="output directory (default: data/osf_guha for the default node; "
                         "required to differ for any other node)")
    ap.add_argument("--conds", nargs="+", default=None, help="conditions to download (default: all)")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    outdir = args.outdir
    if outdir is None:
        outdir = DEFAULT_OUTDIR if args.node == DEFAULT_NODE and args.prefix == DEFAULT_PREFIX \
            else os.path.normpath(os.path.join(HERE, "..", "..", "data", "osf_" + args.node))
    outdir = os.path.normpath(outdir) if os.path.isabs(outdir) else \
        os.path.normpath(os.path.join(HERE, "..", "..", outdir))
    cache_dir = os.path.join(outdir, "cache")
    pat = re.compile(re.escape(args.prefix) + r"_sub-(\d+)_(.+)\.txt")
    ## atlas/pipeline labels for meta.csv: for the untouched default call these come out
    ## EXACTLY as before ("Power229", "osf_guha") -- atlas = prefix up to the first "_", and
    ## pipeline = the output directory's basename ("osf_guha" for the default outdir).
    atlas = args.prefix.split("_")[0]
    pipeline = os.path.basename(outdir.rstrip("/"))

    os.makedirs(cache_dir, exist_ok=True)

    print(f"Node={args.node}  prefix={args.prefix}  outdir={outdir}")
    print("Listing OSF files (explicit page-number enumeration, deduped by file ID) ...")
    allf, raw_seen, n_unique_ids, n_unique_names = list_files(args.node)
    print(f"  self-check: {raw_seen} raw entries seen across all pages; "
          f"{n_unique_ids} unique file IDs; {n_unique_names} unique filenames")
    if n_unique_ids != n_unique_names:
        print(f"  NOTE: unique-ID count != unique-filename count "
              f"({n_unique_ids} vs {n_unique_names}) -- some name collides across distinct IDs")

    want = [n for n in allf if pat.match(n) and (args.conds is None or pat.match(n).group(2) in args.conds)]
    cached_names = set(n for n in want if os.path.exists(os.path.join(cache_dir, n + ".npy")))
    missing = sorted(n for n in want if n not in cached_names)
    print(f"  self-check: target {len(want)} files on OSF matching prefix '{args.prefix}'; "
          f"{len(cached_names)} cached; {len(missing)} missing")

    todo = missing
    ok = 0
    failed = []
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(fetch_one, n, allf[n], cache_dir): n for n in todo}
        done = 0
        for fut in as_completed(futs):
            name, good = fut.result(); done += 1
            if good:
                ok += 1
            else:
                failed.append(name)
            if done % 25 == 0 or done == len(todo):
                print(f"  {done}/{len(todo)} ({ok} new ok, {len(failed)} failed so far)")

    cached_now = len([f for f in os.listdir(cache_dir) if f.endswith(".npy")])
    still_missing = sorted(n for n in want if not os.path.exists(os.path.join(cache_dir, n + ".npy")))
    print(f"  self-check: cache now holds {cached_now} matrices; "
          f"{len(want) - len(still_missing)}/{len(want)} of target cached; "
          f"still MISSING: {len(still_missing)}")
    if still_missing:
        print("  MISSING files (of target):")
        for n in still_missing:
            print(f"    {n}")
    assemble(outdir, cache_dir, pat, atlas, pipeline, conds_filter=args.conds)


if __name__ == "__main__":
    main()
