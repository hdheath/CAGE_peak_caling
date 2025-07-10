#!/usr/bin/env python3
"""
QC routines for TSS peak‐calling pipeline, with caching, conda‐env isolation,
parallelism, detailed logging of bedtools calls and timing, and extended metrics:

1) length of each peak window
2) number of peak windows
3) number of reads in each window (from reads TSS)
4) distances for reads not in windows (strand aware)
5) precision & recall for reads
6) precision & recall for annotation TSS
"""
import os
import subprocess
import time
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import pandas as pd
import numpy as np
import math

# Maximum number of threads for recall calculation
MAX_RECALL_THREADS = 4

# cache for total‐read counts
_peak_read_count_cache = {}


def region_length_dist_file(bed_file, out_dir):
    base = os.path.basename(bed_file).replace(".bed", "")
    out_tsv = os.path.join(out_dir, f"{base}_lengths.tsv")
    if os.path.exists(out_tsv):
        print(f"[QC] Length TSV exists, skipping: {out_tsv}")
        return out_tsv
    print(f"[QC] Computing region_length_dist_file for {bed_file}")
    lengths = []
    with open(bed_file) as fh:
        for line in fh:
            _, start_pos, end_pos, *_ = line.rstrip("\n").split("\t")
            lengths.append(int(end_pos) - int(start_pos))
    with open(out_tsv, "w") as o:
        o.write("length\n")
        for L in lengths:
            o.write(f"{L}\n")
    print(f"[QC] region_length_dist_file wrote {out_tsv}")
    return out_tsv


def count_peak_windows(bed_file):
    with open(bed_file) as f:
        return sum(1 for _ in f)


def extract_tss_from_bam(bam_path, out_bed):
    if os.path.exists(out_bed):
        print(f"[QC] TSS BED exists, skipping: {out_bed}")
        return out_bed
    print(f"[QC] extract_tss_from_bam -> {out_bed}")
    import pysam
    bamf = pysam.AlignmentFile(bam_path, "rb")
    with open(out_bed, "w") as o:
        for r in bamf.fetch(until_eof=True):
            if r.is_unmapped:
                continue
            chrom = bamf.get_reference_name(r.reference_id)
            pos = (r.reference_end - 1) if r.is_reverse else r.reference_start
            strand = '-' if r.is_reverse else '+'
            o.write(f"{chrom}\t{pos}\t{pos+1}\t.\t0\t{strand}\n")
    bamf.close()
    print(f"[QC] Completed extract_tss_from_bam: {out_bed}")
    return out_bed


def extract_annotation_tss(annotation_bed, out_bed):
    if os.path.exists(out_bed):
        print(f"[QC] Annotation TSS BED exists, skipping: {out_bed}")
        return out_bed
    print(f"[QC] extract_annotation_tss -> {out_bed}")
    df = pd.read_csv(annotation_bed, sep="\t", header=None, comment="#").iloc[:, :6]
    rows = []
    for _, row in df.iterrows():
        chrom, start_pos, end_pos, name, score, strand = row
        t0 = int(start_pos) if strand=='+' else int(end_pos)-1
        rows.append((chrom, t0, t0+1, name, score, strand))
    pd.DataFrame(rows).to_csv(out_bed, sep="\t", header=False, index=False)
    print(f"[QC] Completed extract_annotation_tss: {out_bed}")
    return out_bed


@lru_cache(maxsize=None)
def _ensure_sorted(bed, env):
    sorted_bed = bed + ".sorted"
    need_sort = (
        not os.path.exists(sorted_bed)
        or os.path.getsize(sorted_bed) == 0
        or sum(1 for _ in open(sorted_bed)) == 0
    )
    if need_sort:
        print(f"[QC] Sorting BED: {bed} -> {sorted_bed}")
        with open(sorted_bed, 'w') as out:
            cmd = ["conda", "run", "-n", env, "bedtools", "sort", "-i", bed]
            subprocess.run(cmd, stdout=out, check=True)
        print(f"[QC] Sorted BED written: {sorted_bed}")
    else:
        print(f"[QC] Using existing sorted BED: {sorted_bed}")
    return sorted_bed


def _run_closest(a_bed, b_bed, env):
    sa = _ensure_sorted(a_bed, env)
    sb = _ensure_sorted(b_bed, env)
    cmd = ["conda", "run", "-n", env, "bedtools", "closest", "-a", sa, "-b", sb, "-s", "-d"]
    print(f"[QC] Running bedtools closest: {' '.join(cmd)}")
    t0 = time.time()
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
    print(f"[QC] bedtools closest finished in {time.time()-t0:.2f}s, lines={len(proc.stdout.splitlines())}")
    return proc.stdout.splitlines()


def compute_read_recall(reads_bed, peaks_bed, env):
    """Fast unique-read recall via bedtools intersect -u."""
    sr = _ensure_sorted(reads_bed, env)
    sp = _ensure_sorted(peaks_bed, env)

    total = _peak_read_count_cache.get(sr, sum(1 for _ in open(sr)))
    cmd = [
        "conda", "run", "-n", env,
        "bedtools", "intersect", "-a", sr,
        "-b", sp, "-s", "-u"
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True)
    inside = len(proc.stdout.splitlines())

    recall = inside / total if total else 0.0
    print(f"[QC] Read recall = {inside}/{total} = {recall:.3f}")
    return recall


def compute_precision_recall(truth_bed, peak_bed, window, env):
    print(f"[QC] PR for truth={truth_bed}, peaks={peak_bed}, window={window}")
    # Precision via closest -d
    lines = _run_closest(peak_bed, truth_bed, env)
    total_p = len(lines)
    matched_p = sum(1 for ln in lines if abs(int(ln.split('\t')[-1] or 0)) <= window)
    precision = matched_p / total_p if total_p else 0.0

    # Recall via fast intersect -u
    recall = compute_read_recall(truth_bed, peak_bed, env)

    print(f"[QC] Precision={precision:.3f} Recall={recall:.3f}")
    return {'precision': precision, 'recall': recall}


def write_metrics_tsv(metrics, sample, method, out_dir):
    fn = f"{sample}_{method}_metrics.tsv"
    out = os.path.join(out_dir, fn)
    if os.path.exists(out):
        print(f"[QC] Metrics TSV exists, skipping: {out}")
        return out
    pd.DataFrame([{**{'sample': sample, 'method': method}, **metrics}]) \
      .to_csv(out, sep="\t", index=False)
    print(f"[QC] Wrote metrics TSV: {out}")
    return out


def write_peak_window_counts(peak_bed, reads_bed, out_dir, env):
    sorted_peak  = _ensure_sorted(peak_bed, env)
    sorted_reads = _ensure_sorted(reads_bed, env)
    out_file     = os.path.join(
        out_dir, f"{os.path.basename(peak_bed).replace('.bed','')}_read_counts_per_window.tsv"
    )
    if os.path.exists(out_file):
        print(f"[QC] Read‐counts TSV exists, skipping: {out_file}")
        return out_file

    peak_lines = sum(1 for _ in open(sorted_peak))
    if sorted_reads not in _peak_read_count_cache:
        _peak_read_count_cache[sorted_reads] = sum(1 for _ in open(sorted_reads))
    reads_lines = _peak_read_count_cache[sorted_reads]

    print(f"[QC][CHECK] {sorted_peak} → {peak_lines} lines")
    print(f"[QC][CHECK] {sorted_reads} → {reads_lines} lines")

    cmd = [
        "conda", "run", "-n", env, "bedtools", "intersect",
        "-a", sorted_peak, "-b", sorted_reads, "-s", "-c"
    ]
    print(f"[QC] Counting reads per peak window: {' '.join(cmd)} -> {out_file}")
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)

    with open(out_file, "w") as o:
        o.write("Chrom\tStart\tEnd\tName\tScore\tStrand\tReadCount\n")
        for i, line in enumerate(proc.stdout.splitlines(), start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or not fields[-1].isdigit():
                print(f"[WARN] Line {i} malformed or missing count: {line}")
                continue
            o.write(line + "\n")

    print(f"[QC] Wrote read counts per window to: {out_file}")
    return out_file


def write_distances_for_unmatched(reads_bed, peak_bed, out_dir, env, prefix, suffix):
    sorted_reads = _ensure_sorted(reads_bed, env)
    sorted_peak  = _ensure_sorted(peak_bed,  env)
    out_file     = os.path.join(out_dir, f"{prefix}_{suffix}_distances.tsv")
    if os.path.exists(out_file):
        print(f"[QC] Distances TSV exists, skipping: {out_file}")
        return out_file

    cmd = [
        "conda", "run", "-n", env, "bedtools", "closest",
        "-a", sorted_reads, "-b", sorted_peak, "-s", "-d"
    ]
    print(f"[QC] Calculating distances for unmatched reads: {' '.join(cmd)} -> {out_file}")
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)

    header = [
        "Chrom","Start","End","Name","Score","Strand",
        "PeakChrom","PeakStart","PeakEnd","PeakName","PeakScore","PeakStrand","Distance"
    ]
    with open(out_file, "w") as o:
        o.write("\t".join(header) + "\n")
        for i, line in enumerate(proc.stdout.splitlines(), start=1):
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            try:
                dist = int(cols[-1])
            except (ValueError, IndexError):
                print(f"[WARN] Line {i} malformed or missing distance, skipping: {line!r}")
                continue
            if dist > 0:
                o.write(line + "\n")

    print(f"[QC] Wrote unmatched‐reads distances to: {out_file}")
    return out_file


def _qc_task(args):
    sample_dir, bed_name, prefix, lr_tss, annot_tss, window, qc_base, env = args
    peak_bed = os.path.join(sample_dir, bed_name)
    method   = bed_name.replace(".bed", "")

    print(f"[QC][{prefix}][{method}] Starting QC task")
    # 1) region lengths
    region_length_dist_file(peak_bed, qc_base)
    # 2) read counts
    write_peak_window_counts(peak_bed, lr_tss, qc_base, env)
    # 3) distances
    write_distances_for_unmatched(lr_tss, peak_bed, qc_base, env, prefix, f"{method}_readsNotInWindows")
    write_distances_for_unmatched(annot_tss, peak_bed, qc_base, env, prefix, f"{method}_annotNotInWindows")

    # 4) precision & recall
    pr_tsv = os.path.join(qc_base, f"{prefix}_{method}_PR_metrics.tsv")
    if os.path.exists(pr_tsv):
        print(f"[QC] PR TSV exists, skipping precision/recall: {pr_tsv}")
    else:
        print(f"[QC][{prefix}][{method}] Starting precision-recall computation...")
        pr_read  = compute_precision_recall(lr_tss,    peak_bed, window, env)
        pr_annot = compute_precision_recall(annot_tss, peak_bed, window, env)

        metrics = {
            'read_precision':       pr_read ['precision'],
            'read_recall':          pr_read ['recall'],
            'annotation_precision': pr_annot['precision'],
            'annotation_recall':    pr_annot['recall'],
        }
        print(f"[QC][{prefix}][{method}] PR metrics → {metrics}")
        write_metrics_tsv(metrics, prefix, f"{method}_PR", qc_base)

    print(f"[QC][{prefix}][{method}] Finished QC task")


def run_qc_parallel(sample_dir, beds_to_qc, prefix, lr_tss, annot_tss, window, qc_base, env, max_workers=6):
    print(f"[QC] Running parallel QC with up to {max_workers} workers")
    tasks = [(sample_dir, b, prefix, lr_tss, annot_tss, window, qc_base, env) for b in beds_to_qc]
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        exe.map(_qc_task, tasks)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="QC pipeline with extended metrics")
    parser.add_argument("sample_dir")
    parser.add_argument("beds_to_qc", nargs='+')
    parser.add_argument("prefix")
    parser.add_argument("lr_tss")
    parser.add_argument("annot_tss")
    parser.add_argument("window", type=int)
    parser.add_argument("qc_base")
    parser.add_argument("env")
    parser.add_argument("--max_workers", type=int, default=6)
    args = parser.parse_args()
    run_qc_parallel(
        args.sample_dir, args.beds_to_qc, args.prefix,
        args.lr_tss, args.annot_tss, args.window,
        args.qc_base, args.env, args.max_workers
    )
