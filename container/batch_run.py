"""Batch trial over the VUMC BIDS cohort inside the pet-fdg container.

Picks one whole-body AC series per subject, runs the full runner() chain
against /templates, and writes one aggregated summary CSV. Per-series
failures are logged, never fatal.

    pet-fdg batch /bids --pattern '*3DWB*AC*PT.nii.gz' -n 12
"""

import argparse
import csv
import re
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, "/opt/src/pet-py/tests")
from run import runner  # noqa: E402

import pet_py  # noqa: E402


def pick_series(bids_root: Path, pattern: str, limit: int, all_subjects: bool):
    seen = set()
    out = []
    for nii in sorted(bids_root.glob(f"sub-*/PT/{pattern}")):
        mrn = re.search(r"sub-(\d{7,})_", nii.name).group(1)
        if mrn in seen and not all_subjects:
            continue
        seen.add(mrn)
        out.append(nii)
        if len(out) >= limit:
            break
    return out


def main() -> None:
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("bids_root", type=Path)
    ap.add_argument("--pattern", default="*AC*3DWB*PT.nii.gz")
    ap.add_argument("-n", type=int, default=10)
    ap.add_argument("--all-subjects", action="store_true")
    ap.add_argument(
        "--line",
        type=int,
        help="run only the 1-based K-th selected series (SLURM array task mode)",
    )
    ap.add_argument(
        "--head-at",
        choices=("high", "low"),
        help="force the head-end for robustfov on ambiguous-orientation "
        "series (triage override; default: content heuristic)",
    )
    ap.add_argument("--work", type=Path, default=Path("/output/batch"))
    ap.add_argument("--templates", type=Path, default=Path("/templates"))
    ap.add_argument("--summary", type=Path, default=Path("/output/batch_summary.csv"))
    args = ap.parse_args()

    series = pick_series(args.bids_root, args.pattern, args.n, args.all_subjects)
    if args.line:
        if not 1 <= args.line <= len(series):
            raise SystemExit(f"--line {args.line} outside 1..{len(series)}")
        series = series[args.line - 1:args.line]
    print(f"{len(series)} series selected", flush=True)

    rows = []
    for src in series:
        mrn = re.search(r"sub-(\d{7,})_", src.name).group(1)
        stage = args.work / src.parent.parent.name / "PT"
        stage.mkdir(parents=True, exist_ok=True)
        t0 = time.time()
        row = {"mrn": mrn, "series": src.name, "status": "ok", "secs": "", "csv": ""}
        try:
            for src_file, dst in (
                (src, stage / src.name),
                (src.with_name(src.name.replace(".nii.gz", ".json")),
                 stage / src.name.replace(".nii.gz", ".json")),
            ):
                if not dst.is_file():
                    shutil.copy(src_file, dst)
            js = pet_py.analysis.get_suvbw_scale_factor(stage / src.name.replace(".nii.gz", ".json"))
            runner(stage / src.name, args.work, args.templates,
                   args.templates / "param27.txt", head_at=args.head_at)
            row["scale"] = round(js, 6)
            outs = sorted(args.work.glob(f"derivatives/sub-{mrn}/**/*SUVR.csv"))
            row["csv"] = str(outs[-1]) if outs else ""
            if not outs:
                row["status"] = "no-csv"
        except Exception as e:
            row["status"] = f"FAIL: {type(e).__name__}: {str(e)[:160]}"
        row["secs"] = round(time.time() - t0, 1)
        print(mrn, row["status"], flush=True)
        rows.append(row)

    cols = ["mrn", "series", "status", "secs", "scale", "csv"]
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    with open(args.summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("summary:", args.summary)


if __name__ == "__main__":
    main()
