#!/usr/bin/env python3
"""
Cluster 5′-end coverage (start sites) from an indexed BAM
and write BED files of peaks via DBSCAN, HDBSCAN, scipy.find_peaks,
MACS3, or DBSCAN+sliding-window.
"""
import os
import sys
import argparse
import warnings
import subprocess

import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import hdbscan
from scipy.signal import find_peaks


def parse_args():
    p = argparse.ArgumentParser(
        description="Cluster strand-aware 5′-end coverage from a BAM."
    )
    p.add_argument(
        "bam",
        help="Indexed BAM file (e.g. sample.sorted.bam)"
    )
    p.add_argument(
        "out_dir",
        help="Directory to write output BED files (created if needed)"
    )
    p.add_argument(
        "--chrom",
        help="Restrict to a single chromosome (e.g. 'chr1')",
        default=None
    )

    # methods
    p.add_argument("--run_dbscan",       action="store_true", help="Enable DBSCAN")
    p.add_argument("--dbscan_eps",       nargs="+", type=float, help="ε values for DBSCAN")
    p.add_argument("--dbscan_min_samples", nargs="+", type=int, help="min_samples for DBSCAN")

    p.add_argument("--run_hdbscan",      action="store_true", help="Enable HDBSCAN")
    p.add_argument("--hdbscan_min_cluster_size", nargs="+", type=int, help="min_cluster_size for HDBSCAN")
    p.add_argument("--hdbscan_min_samples", type=int, default=1, help="min_samples for HDBSCAN")

    p.add_argument("--run_find_peaks",   action="store_true", help="Enable scipy.find_peaks")
    p.add_argument("--peak_height",      nargs="+", type=int, help="peak height(s)")
    p.add_argument("--peak_distance",    nargs="+", type=int, help="peak distance(s)")

    p.add_argument("--run_macs3",        action="store_true", help="Enable MACS3")
    p.add_argument("--macs3_gsize",      default="hs", help="Genome size for MACS3")
    p.add_argument("--macs3_qvalue",     type=float, default=0.01, help="q-value for MACS3")

    p.add_argument("--run_dbscan_sw",    action="store_true", help="Enable DBSCAN+sliding-window")
    p.add_argument("--sw_bin",           nargs="+", type=int, help="sliding-window bin sizes")
    p.add_argument("--sw_threshold",     nargs="+", type=int, help="read-count thresholds per bin")

    return p.parse_args()


def build_df_from_bam(bam_path, chrom_filter=None):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    recs = []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        chrom = bam.get_reference_name(read.reference_id)
        if chrom_filter and chrom != chrom_filter:
            continue
        if read.is_reverse:
            # 5′ end on reverse‐strand is reference_end - 1
            pos = read.reference_end - 1
            strand = "-"
        else:
            pos = read.reference_start
            strand = "+"
        recs.append((chrom, pos, strand))
    bam.close()
    df = (
        pd.DataFrame(recs, columns=["chrom", "pos", "strand"])
          .groupby(["chrom", "pos", "strand"], sort=False)
          .size()
          .reset_index(name="freq")
    )
    return df


def write_bed(clusters, out_path):
    if not clusters:
        warnings.warn(f"No clusters to write for {out_path}")
        return
    with open(out_path, "w") as f:
        for cl in clusters:
            name = f"{cl['method']}_{cl['param']}"
            f.write(
                f"{cl['chrom']}\t{cl['start']}\t{cl['end']}\t"
                f"{name}\t{cl['support']}\t{cl['strand']}\n"
            )


def select_dbscan_params(df):
    data = df["pos"].values.reshape(-1, 1)
    if data.shape[0] < 2:
        return 5, 2
    best_score, best_eps, best_ms = -1, 5, 2
    for eps in range(1, 21):
        for ms in range(2, 6):
            db = DBSCAN(eps=eps, min_samples=ms, metric="chebyshev").fit(data)
            lbl = db.labels_
            mask = lbl != -1
            if mask.sum() < 2 or len(np.unique(lbl[mask])) < 2:
                continue
            try:
                score = silhouette_score(data[mask], lbl[mask], metric="chebyshev")
            except:
                continue
            if score > best_score:
                best_score, best_eps, best_ms = score, eps, ms
    return best_eps, best_ms


def run_dbscan(df, eps, ms):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        pts = grp["pos"].values.reshape(-1, 1)
        db = DBSCAN(eps=eps, min_samples=ms, metric="chebyshev").fit(pts)
        grp2 = grp.assign(cluster=db.labels_)
        for cid, sub in grp2[grp2.cluster != -1].groupby("cluster"):
            start = int(sub.pos.min())
            end = int(sub.pos.max()) + 1
            support = int(sub.freq.sum())
            clusters.append({
                "chrom": chrom, "start": start, "end": end,
                "support": support, "strand": strand,
                "method": "DBSCAN", "param": f"eps{eps}_ms{ms}"
            })
    return clusters


def run_hdbscan(df, mcs, ms):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        pts = grp["pos"].values.reshape(-1, 1)
        cl = hdbscan.HDBSCAN(
            min_cluster_size=mcs,
            min_samples=ms,
            metric="chebyshev"
        ).fit(pts)
        grp2 = grp.assign(cluster=cl.labels_)
        for cid, sub in grp2[grp2.cluster != -1].groupby("cluster"):
            start = int(sub.pos.min())
            end = int(sub.pos.max()) + 1
            support = int(sub.freq.sum())
            clusters.append({
                "chrom": chrom, "start": start, "end": end,
                "support": support, "strand": strand,
                "method": "HDBSCAN", "param": f"mcs{mcs}_ms{ms}"
            })
    return clusters


def run_find_peaks(df, height, distance, chrom_lengths):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        length = chrom_lengths.get(chrom, grp.pos.max() + 1)
        cov = np.zeros(length, dtype=int)
        for _, r in grp.iterrows():
            cov[int(r.pos)] = int(r.freq)
        cov_sm = np.convolve(cov, np.ones(3)/3.0, mode="same")
        peaks, _ = find_peaks(cov_sm, height=height, distance=distance)
        for p in peaks:
            left, right = p, p
            while left > 0 and cov_sm[left] >= height/2:
                left -= 1
            while right < length - 1 and cov_sm[right] >= height/2:
                right += 1
            support = int(cov[left:right+1].sum())
            clusters.append({
                "chrom": chrom, "start": left, "end": right + 1,
                "support": support, "strand": strand,
                "method": "findPeaks", "param": f"h{height}_d{distance}"
            })
    return clusters


def run_macs3(df, gsize, qvalue, out_dir):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        tmp = os.path.join(out_dir, f"tmp_{chrom}_{strand}.bed")
        with open(tmp, "w") as f:
            for _, r in grp.iterrows():
                f.write(f"{r.chrom}\t{r.pos}\t{r.pos+1}\t.\t{r.freq}\t{r.strand}\n")

        for mode in ("narrow", "broad"):
            prefix = f"MACS3_{chrom}_{strand}_{mode}"
            cmd = [
                "macs3", "callpeak",
                "-f", "BED", "-t", tmp,
                "-n", prefix,
                "--gsize", str(gsize),
                "--qvalue", str(qvalue),
                "--nomodel", "--keep-dup", "all",
                "--outdir", out_dir
            ]
            if mode == "broad":
                cmd += ["--broad", "--broad-cutoff", str(qvalue)]
            subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL)

            peakfile = os.path.join(out_dir, f"{prefix}_peaks.{mode}Peak")
            if os.path.isfile(peakfile):
                with open(peakfile) as pf:
                    for line in pf:
                        cols = line.strip().split("\t")
                        c, s, e, _, sc, st = cols[:6]
                        clusters.append({
                            "chrom": c, "start": int(s), "end": int(e),
                            "support": int(float(sc)), "strand": st,
                            "method": f"MACS3_{mode}", "param": f"g{gsize}_q{qvalue}"
                        })
        os.remove(tmp)
    return clusters


def run_dbscan_sw(df, eps, ms, bin_size, threshold, chrom_lengths):
    clusters = run_dbscan(df, eps, ms)
    covered = {
        (cl["chrom"], cl["strand"], p)
        for cl in clusters
        for p in range(cl["start"], cl["end"])
    }
    residual = df[~df.apply(lambda r: (r.chrom, r.strand, r.pos) in covered, axis=1)]

    for (chrom, strand), grp in residual.groupby(["chrom", "strand"]):
        length = chrom_lengths.get(chrom, grp.pos.max() + 1)
        nbin = length // bin_size + 1
        hist = np.zeros(nbin, dtype=int)
        for _, r in grp.iterrows():
            hist[int(r.pos) // bin_size] += int(r.freq)
        hot = hist >= threshold

        i = 0
        while i < nbin:
            if hot[i]:
                j, supp = i, 0
                while j < nbin and hot[j]:
                    supp += hist[j]
                    j += 1
                clusters.append({
                    "chrom": chrom,
                    "start": i * bin_size,
                    "end": min(j * bin_size, length),
                    "support": supp,
                    "strand": strand,
                    "method": "DBSCAN_SW",
                    "param": f"eps{eps}_bin{bin_size}_th{threshold}"
                })
                i = j
            else:
                i += 1

    return clusters


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # build coverage df of 5′-ends
    df = build_df_from_bam(args.bam, args.chrom)
    if df.empty:
        sys.exit("No 5′-end coverage found. Exiting.")

    # precompute chromosome lengths
    chrom_lengths = {c: df[df.chrom == c].pos.max() + 1 for c in df.chrom.unique()}

    # 1) DBSCAN
    if args.run_dbscan:
        eps_list = args.dbscan_eps or [select_dbscan_params(df)[0]]
        ms_list = args.dbscan_min_samples or [select_dbscan_params(df)[1]]
        for eps in eps_list:
            for ms in ms_list:
                cls = run_dbscan(df, eps, ms)
                write_bed(cls, os.path.join(args.out_dir, f"DBSCAN_eps{eps}_ms{ms}.bed"))

    # 2) HDBSCAN
    if args.run_hdbscan and args.hdbscan_min_cluster_size:
        for mcs in args.hdbscan_min_cluster_size:
            cls = run_hdbscan(df, mcs, args.hdbscan_min_samples)
            write_bed(cls, os.path.join(
                args.out_dir,
                f"HDBSCAN_mcs{mcs}_ms{args.hdbscan_min_samples}.bed"
            ))

    # 3) find_peaks
    if args.run_find_peaks and args.peak_height and args.peak_distance:
        for h in args.peak_height:
            for d in args.peak_distance:
                cls = run_find_peaks(df, h, d, chrom_lengths)
                write_bed(cls, os.path.join(
                    args.out_dir,
                    f"findPeaks_h{h}_d{d}.bed"
                ))

    # 4) MACS3
    if args.run_macs3:
        cls = run_macs3(df, args.macs3_gsize, args.macs3_qvalue, args.out_dir)
        write_bed(cls, os.path.join(args.out_dir, "MACS3_all.bed"))

    # 5) DBSCAN + sliding window
    if args.run_dbscan_sw and args.sw_bin and args.sw_threshold:
        eps, ms = (args.dbscan_eps or [5])[0], (args.dbscan_min_samples or [2])[0]
        for b in args.sw_bin:
            for thr in args.sw_threshold:
                cls = run_dbscan_sw(df, eps, ms, b, thr, chrom_lengths)
                write_bed(cls, os.path.join(
                    args.out_dir,
                    f"DBSCAN_SW_eps{eps}_bin{b}_th{thr}.bed"
                ))


if __name__ == "__main__":
    main()
