#!/usr/bin/env python3
"""
Utility functions: argument parsing, BAM-based DataFrame builder,
BED writer, and DBSCAN parameter selector.
"""
import argparse
import warnings
import pandas as pd
import pysam
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score

def parse_args():
    p = argparse.ArgumentParser(
        description="Cluster strand-aware 5′-end coverage from a BAM."
    )
    p.add_argument(
        "--config",
        help="Path to a TOML config file with all inputs & flags",
        default=None
    )
    p.add_argument(
        "bam",
        nargs="*",
        help="One or more indexed BAMs (e.g. sample1.bam sample2.bam), or supply via config.toml as a TOML array"
    )
    p.add_argument(
        "out_dir",
        nargs="?",
        help="Directory to write output BED files — or supply via config.toml"
    )
    p.add_argument(
        "--chrom",
        help="Restrict to a single chromosome (e.g. 'chr1')",
        default=None
    )
    p.add_argument("--cager_threshold", nargs="+", type=int, default=[], help="Minimum tag count(s) for CAGEr")
    p.add_argument("--cager_maxdist",   nargs="+", type=int, default=[], help="Maximum cluster distance(s) for CAGEr")
    
    p.add_argument(
        "--dry_run",
        action="store_true",
        help="Print what would be run/written, but don’t actually execute anything"
    )

    # methods
    p.add_argument("--run_dbscan", action="store_true", help="Enable DBSCAN")
    p.add_argument(
        "--dbscan_eps",
        nargs="+",
        type=float,
        help="ε values for DBSCAN"
    )
    p.add_argument(
        "--dbscan_min_samples",
        nargs="+",
        type=int,
        help="min_samples for DBSCAN"
    )

    p.add_argument("--run_hdbscan", action="store_true", help="Enable HDBSCAN")
    p.add_argument(
        "--hdbscan_min_cluster_size",
        nargs="+",
        type=int,
        help="min_cluster_size for HDBSCAN"
    )
    p.add_argument(
        "--hdbscan_min_samples",
        type=int,
        default=1,
        help="min_samples for HDBSCAN"
    )

    p.add_argument("--run_find_peaks", action="store_true", help="Enable scipy.find_peaks")
    p.add_argument(
        "--peak_height",
        nargs="+",
        type=int,
        help="peak height(s)"
    )
    p.add_argument(
        "--peak_distance",
        nargs="+",
        type=int,
        help="peak distance(s)"
    )



    p.add_argument("--run_dbscan_sw", action="store_true", help="Enable DBSCAN+sliding-window")
    p.add_argument(
        "--sw_bin",
        nargs="+",
        type=int,
        help="sliding-window bin sizes"
    )
    p.add_argument(
        "--sw_threshold",
        nargs="+",
        type=int,
        help="read-count thresholds per bin"
    )
    p.add_argument(
        "--macs3_gsize",
        nargs="+",
        default=[],
        help="Genome size(s) for MACS3 (e.g. hs mm dm)"
    )
    p.add_argument(
        "--macs3_qvalue",
        nargs="+",
        type=float,
        default=[],
        help="q-value(s) for MACS3 (e.g. 0.01 0.05)"
    )
    # after your existing macs3 args
    p.add_argument(
        "--macs3_nomodel", nargs="*", type=lambda x: x.lower() == "true",
        help="List of booleans for --nomodel"
    )
    p.add_argument(
        "--macs3_shift", nargs="*", type=int,
        help="List of shift values for MACS3"
    )
    p.add_argument(
        "--macs3_extsize", nargs="*", type=int,
        help="List of extension sizes for MACS3"
    )
    p.add_argument(
        "--macs3_call_summits", nargs="*", type=lambda x: x.lower() == "true",
        help="List of booleans for --call-summits"
    )

    # and for CAGEr
    p.add_argument(
        "--cager_correctFirstG", nargs="*", type=lambda x: x.lower() == "true",
        help="List of booleans for correctFirstG"
    )


    # QC options
    p.add_argument("--run_qc",           action="store_true", help="Enable QC at end")
    p.add_argument("--qc_out_dir",       help="Directory for QC outputs")
    p.add_argument("--long_reads",       nargs="+",             help="Long-read BAM(s)")
    p.add_argument("--annotation_bed",   help="Annotation BED file")
    p.add_argument("--bedtools_env",     default="bedtools_env", help="Conda env with bedtools")
    p.add_argument("--qc_dist_cutoff",   type=int, default=50,  help="Distance cutoff for bedtools closest")

    return p.parse_args()


def build_df_from_bam(bam_paths, chrom_filter=None):
    """Count 5′ ends per (chrom,pos,strand) across one or more BAMs."""
    recs = []
    for bam_path in bam_paths:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            if chrom_filter and chrom != chrom_filter:
                continue
            if read.is_reverse:
                pos, strand = read.reference_end - 1, "-"
            else:
                pos, strand = read.reference_start, "+"
            recs.append((chrom, pos, strand))
        bam.close()
    # aggregate counts across all BAMs
    df = (
        pd.DataFrame(recs, columns=["chrom","pos","strand"])
          .groupby(["chrom","pos","strand"], sort=False)
          .size()
          .reset_index(name="freq")
    )
    return df



def write_bed(clusters, out_path):
    """Write either a list of cluster dicts *or* a DataFrame to a BED6 file."""
    # If it's a DataFrame, assume it already has the right columns
    if isinstance(clusters, pd.DataFrame):
        df = clusters
        if df.empty:
            warnings.warn(f"No clusters to write for {out_path}")
            return
        with open(out_path, "w") as f:
            # Expecting columns: chrom/seqnames, start, end, name, score, strand
            # Adjust `seqnames`→`chrom` if necessary
            chrom_col = "chrom" if "chrom" in df.columns else "seqnames"
            for _, row in df.iterrows():
                f.write(
                    f"{row[chrom_col]}\t{row['start']}\t{row['end']}\t"
                    f"{row['name']}\t{row['score']}\t{row['strand']}\n"
                )
        return

    # Otherwise expect a list of dicts with keys chrom/start/end/name/score/strand
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
    """Grid-search eps & min_samples via silhouette_score."""
    data = df["pos"].values.reshape(-1, 1)
    if data.shape[0] < 2:
        return 5, 2
    best_score, best_eps, best_ms = -1, 5, 2
    for eps in range(1, 21):
        for ms in range(2, 6):
            db = DBSCAN(eps=eps, min_samples=ms, metric="chebyshev").fit(data)
            lbl = db.labels_
            mask = lbl != -1
            if mask.sum() < 2 or len(pd.unique(lbl[mask])) < 2:
                continue
            try:
                score = silhouette_score(data[mask], lbl[mask], metric="chebyshev")
            except:
                continue
            if score > best_score:
                best_score, best_eps, best_ms = score, eps, ms
    return best_eps, best_ms
