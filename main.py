#!/usr/bin/env python3
"""
Driver script: imports all modules and orchestrates execution

Synopsis:
  python main.py --config config.toml --dry_run
  python main.py sample1.bam sample2.bam output_dir [flags…]
"""
import os
import sys
import warnings
import pandas as pd

from utils       import parse_args, build_df_from_bam, write_bed, select_dbscan_params
from clustering  import run_dbscan, run_hdbscan, run_find_peaks, run_dbscan_sw
from macs3       import run_macs3
from cager       import run_cager
from config      import load_config


def _filter_bed_by_chrom(in_bed, out_bed, chrom):
    """Write only lines from in_bed that start with chrom to out_bed."""
    with open(in_bed) as inf, open(out_bed, "w") as outf:
        for line in inf:
            if line.startswith(f"{chrom}\t"):
                outf.write(line)
    return out_bed


def print_dry_run_plan(args, sample_dir, prefix):
    import os
    print(f"[DRY RUN][{prefix}] Planned actions:")

    # 1) DBSCAN
    if args.run_dbscan:
        for eps in (args.dbscan_eps or ["<auto_eps>"]):
            for ms in (args.dbscan_min_samples or ["<auto_ms>"]):
                fn = f"{prefix}_DBSCAN_eps{eps}_ms{ms}.bed"
                print(f"  [DBSCAN] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [DBSCAN] skipped")

    # 2) HDBSCAN
    if args.run_hdbscan:
        if not args.hdbscan_min_cluster_size:
            print("  [HDBSCAN] skipped (no min_cluster_size provided)")
        else:
            for mcs in args.hdbscan_min_cluster_size:
                fn = f"{prefix}_HDBSCAN_mcs{mcs}_ms{args.hdbscan_min_samples}.bed"
                print(f"  [HDBSCAN] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [HDBSCAN] skipped")

    # 3) find_peaks
    if args.run_find_peaks:
        if not (args.peak_height and args.peak_distance):
            print("  [find_peaks] skipped (missing height or distance)")
        else:
            for h in args.peak_height:
                for d in args.peak_distance:
                    fn = f"{prefix}_findPeaks_h{h}_d{d}.bed"
                    print(f"  [find_peaks] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [find_peaks] skipped")

    # 4) MACS3
    if args.run_macs3:
        for g in args.macs3_gsize or ["hs"]:
            for q in args.macs3_qvalue or [0.01]:
                fn = f"{prefix}_MACS3_g{g}_q{q}.bed"
                print(f" [MACS3] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [MACS3] skipped")

    # 5) DBSCAN + sliding-window
    if args.run_dbscan_sw:
        if not (args.sw_bin and args.sw_threshold):
            print("  [DBSCAN_SW] skipped (missing bin or threshold)")
        else:
            base_eps = (args.dbscan_eps or [5])[0]
            for b in args.sw_bin:
                for thr in args.sw_threshold:
                    fn = f"{prefix}_DBSCAN_SW_eps{base_eps}_bin{b}_th{thr}.bed"
                    print(f"  [DBSCAN_SW] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [DBSCAN_SW] skipped")

    # 6) CAGEr
    if args.run_cager:
        for thr in (args.cager_threshold or [1]):
            for md in (args.cager_maxdist or [20]):
                fn = f"{prefix}_CAGEr_thr{thr}_md{md}.bed"
                print(f"  [CAGEr] would write {os.path.join(sample_dir, fn)}")
    else:
        print("  [CAGEr] skipped")
            
    # 7) QC
    if args.run_qc:
        qc_base = args.qc_out_dir or os.path.join(args.out_dir, "QC")
        print(f"  [QC] would create QC directory → {qc_base}")

        # a) TSS extraction steps (one per sample, not per method)
        print(f"  [QC] would extract long-read TSS → {os.path.join(qc_base, f'{prefix}_longReads_TSS.bed')}")
        print(f"  [QC] would extract annotation TSS → {os.path.join(qc_base, 'annotation_TSS.bed')}")

        # b) per‐method region‐length & closest runs
        # reconstruct the exact BED filenames you planned above
        planned = []
        if args.run_dbscan:
            for eps in (args.dbscan_eps or ["<auto_eps>"]):
                for ms in (args.dbscan_min_samples or ["<auto_ms>"]):
                    planned.append(f"{prefix}_DBSCAN_eps{eps}_ms{ms}.bed")
        if args.run_hdbscan:
            for mcs in (args.hdbscan_min_cluster_size or []):
                planned.append(f"{prefix}_HDBSCAN_mcs{mcs}_ms{args.hdbscan_min_samples}.bed")
        if args.run_find_peaks:
            for h in (args.peak_height or []):
                for d in (args.peak_distance or []):
                    planned.append(f"{prefix}_findPeaks_h{h}_d{d}.bed")
        if args.run_macs3:
            planned.append(f"{prefix}_MACS3_all.bed")
        if args.run_dbscan_sw:
            base_eps = (args.dbscan_eps or [5])[0]
            for b in (args.sw_bin or []):
                for thr in (args.sw_threshold or []):
                    planned.append(f"{prefix}_DBSCAN_SW_eps{base_eps}_bin{b}_th{thr}.bed")
        if args.run_cager:
            planned.append(f"{prefix}_CAGEr_clusters.bed")

        for bed in planned:
            # region lengths
            lengths_tsv = f"{bed.replace('.bed','')}_lengths.tsv"
            print(f"  [QC] would compute region lengths → {os.path.join(qc_base, lengths_tsv)}")
            # long-read vs peaks
            lr_dist = f"{bed.replace('.bed','')}_lr_dist.tsv"
            print(f"  [QC] would run long-read vs {bed} → {os.path.join(qc_base, lr_dist)}")
            # annotation vs peaks
            annot_dist = f"{bed.replace('.bed','')}_annot_dist.tsv"
            print(f"  [QC] would run annotation vs {bed} → {os.path.join(qc_base, annot_dist)}")
    else:
        print("  [QC] skipped")

    print(f"[DRY RUN][{prefix}] End of plan.\n")



def main():
    cli = parse_args()

    # load TOML config
    cfg = {}
    if cli.config:
        if not os.path.isfile(cli.config):
            sys.exit(f"Config file not found: {cli.config}")
        cfg = load_config(cli.config)
        print(f"[CONFIG] loaded {cli.config}")

    # overlay config then CLI
    args = cli
    # config may use 'bams' or 'bam'
    if 'bams' in cfg:
        args.bam = cfg['bams']
    elif 'bam' in cfg:
        args.bam = cfg['bam'] if isinstance(cfg['bam'], list) else [cfg['bam']]
    # override out_dir if provided in config and not on CLI
    if 'out_dir' in cfg and not args.out_dir:
        args.out_dir = cfg['out_dir']
    # apply any other flags from config if CLI didn't set them
    for k, v in cfg.items():
        if k in ('bams','bam','out_dir'):
            continue
        if getattr(args, k, None) in (None, False, "", []):
            setattr(args, k, v)

    # ensure inputs
    if not getattr(args, 'bam', None) or not args.out_dir:
        sys.exit("Error: must specify one or more BAM(s) and an out_dir via CLI or config")

    # normalize to list
    bam_list = args.bam if isinstance(args.bam, list) else [args.bam]
    os.makedirs(args.out_dir, exist_ok=True)

    # dry-run
    if args.dry_run:
        for bam_path in bam_list:
            prefix     = os.path.basename(bam_path).split('.',1)[0]
            sample_dir = os.path.join(args.out_dir, prefix)
            print_dry_run_plan(args, sample_dir, prefix)
        sys.exit(0)

    # process each BAM
    for bam_path in bam_list:
        prefix     = os.path.basename(bam_path).split('.',1)[0]
        sample_dir = os.path.join(args.out_dir, prefix)
        os.makedirs(sample_dir, exist_ok=True)
        print(f"[RUNNING] Sample: {prefix}")

        # coverage cache
        cov_cache = os.path.join(sample_dir, f"{prefix}_coverage.tsv")
        if os.path.exists(cov_cache):
            df = pd.read_csv(cov_cache, sep='\t')
            print(f"[CACHE][{prefix}] Loaded coverage from {cov_cache}")
        else:
            df = build_df_from_bam([bam_path], args.chrom)
            if df.empty:
                print(f"[WARNING] No coverage in {bam_path}; skipping.")
                continue
            df.to_csv(cov_cache, sep='\t', index=False)
            print(f"[CACHE][{prefix}] Wrote coverage cache → {cov_cache}")

        # chromosome lengths map
        chrom_lengths = {c: df[df.chrom==c].pos.max()+1 for c in df.chrom.unique()}

        # DBSCAN
        if args.run_dbscan:
            for eps in (args.dbscan_eps or [select_dbscan_params(df)[0]]):
                for ms in (args.dbscan_min_samples or [select_dbscan_params(df)[1]]):
                    fn      = f"{prefix}_DBSCAN_eps{eps}_ms{ms}.bed"
                    out_bed = os.path.join(sample_dir, fn)
                    if os.path.exists(out_bed):
                        print(f"[DBSCAN][{prefix}] {fn} exists; skipping.")
                    else:
                        clusters = run_dbscan(df, eps, ms)
                        write_bed(clusters, out_bed)
                        print(f"[DBSCAN][{prefix}] Wrote {len(clusters)} clusters to {fn}")
        else:
            print(f"[DBSCAN][{prefix}] Skipped")

        # HDBSCAN
        if args.run_hdbscan:
            if not args.hdbscan_min_cluster_size:
                warnings.warn(f"[HDBSCAN][{prefix}] missing --hdbscan_min_cluster_size; skipping.")
            else:
                for mcs in args.hdbscan_min_cluster_size:
                    fn      = f"{prefix}_HDBSCAN_mcs{mcs}_ms{args.hdbscan_min_samples}.bed"
                    out_bed = os.path.join(sample_dir, fn)
                    if os.path.exists(out_bed):
                        print(f"[HDBSCAN][{prefix}] {fn} exists; skipping.")
                    else:
                        clusters = run_hdbscan(df, mcs, args.hdbscan_min_samples)
                        write_bed(clusters, out_bed)
                        print(f"[HDBSCAN][{prefix}] Wrote {len(clusters)} clusters to {fn}")
        else:
            print(f"[HDBSCAN][{prefix}] Skipped")

        # find_peaks
        if args.run_find_peaks:
            if not (args.peak_height and args.peak_distance):
                warnings.warn(f"[find_peaks][{prefix}] missing height/distance; skipping.")
            else:
                for h in args.peak_height:
                    for d in args.peak_distance:
                        fn      = f"{prefix}_findPeaks_h{h}_d{d}.bed"
                        out_bed = os.path.join(sample_dir, fn)
                        if os.path.exists(out_bed):
                            print(f"[find_peaks][{prefix}] {fn} exists; skipping.")
                        else:
                            clusters = run_find_peaks(df, h, d, chrom_lengths)
                            write_bed(clusters, out_bed)
                            print(f"[find_peaks][{prefix}] Wrote {len(clusters)} peaks to {fn}")
        else:
            print(f"[find_peaks][{prefix}] Skipped")

        # MACS3
        if args.run_macs3:
            for g in args.macs3_gsize or ["hs"]:
                for q in args.macs3_qvalue or [0.01]:
                    fn = f"{prefix}_MACS3_g{g}_q{q}.bed"
                    out_bed = os.path.join(sample_dir, fn)
                    if os.path.exists(out_bed):
                        print(f"[MACS3][{prefix}] {fn} exists; skipping.")
                    else:
                        clusters = run_macs3(df, g, q, sample_dir)
                        write_bed(clusters, out_bed)
                        print(f"[MACS3][{prefix}] Wrote {len(clusters)} peaks to {fn}")

        else:
            print(f"[MACS3][{prefix}] Skipped")

        # DBSCAN + sliding-window
        if args.run_dbscan_sw:
            if not (args.sw_bin and args.sw_threshold):
                warnings.warn(f"[DBSCAN_SW][{prefix}] missing bin/threshold; skipping.")
            else:
                eps0, ms0 = (args.dbscan_eps or [5])[0], (args.dbscan_min_samples or [2])[0]
                for b in args.sw_bin:
                    for thr in args.sw_threshold:
                        fn      = f"{prefix}_DBSCAN_SW_eps{eps0}_bin{b}_th{thr}.bed"
                        out_bed = os.path.join(sample_dir, fn)
                        if os.path.exists(out_bed):
                            print(f"[DBSCAN_SW][{prefix}] {fn} exists; skipping.")
                        else:
                            clusters = run_dbscan_sw(df, eps0, ms0, b, thr, chrom_lengths)
                            write_bed(clusters, out_bed)
                            print(f"[DBSCAN_SW][{prefix}] Wrote {len(clusters)} regions to {fn}")
        else:
            print(f"[DBSCAN_SW][{prefix}] Skipped")

        # CAGEr
        if args.run_cager:
            fn      = f"{prefix}_CAGEr_clusters.bed"
            out_bed = os.path.join(sample_dir, fn)
            if os.path.exists(out_bed):
                print(f"[CAGEr][{prefix}] {fn} exists; skipping.")
            else:
                try:
                    bed_path = run_cager(
                     bam_files=[bam_path],
                     out_dir=sample_dir,
                     threshold=args.cager_threshold,
                     maxdist=args.cager_maxdist,
                     env=args.cager_env
                    )
                    print(f"[CAGEr][{prefix}] Clusters written to {os.path.basename(bed_path)}")
                except Exception as e:
                    warnings.warn(f"[CAGEr][{prefix}] Failed: {e}")
        else:
            print(f"[CAGEr][{prefix}] Skipped")

        # ------------------------------------------------------
        # QC
        if args.run_qc:
            import qc
            qc_base = args.qc_out_dir or os.path.join(args.out_dir, 'QC')
            os.makedirs(qc_base, exist_ok=True)

            # 1) extract TSS beds once
            lr_tss_raw    = os.path.join(qc_base, f"{prefix}_longReads_TSS.bed")
            annot_tss_raw = os.path.join(qc_base, "annotation_TSS.bed")
            qc.extract_tss_from_bam(bam_path, lr_tss_raw)
            qc.extract_annotation_tss(args.annotation_bed, annot_tss_raw)

            # 2) filter both truth beds by args.chrom if set
            if getattr(args, 'chrom', None):
                filtered_lr    = os.path.join(qc_base, f"{prefix}_longReads_TSS.{args.chrom}.bed")
                filtered_annot = os.path.join(qc_base, f"annotation_TSS.{args.chrom}.bed")
                lr_tss    = _filter_bed_by_chrom(lr_tss_raw,    filtered_lr,    args.chrom)
                annot_tss = _filter_bed_by_chrom(annot_tss_raw, filtered_annot, args.chrom)
            else:
                lr_tss, annot_tss = lr_tss_raw, annot_tss_raw

            # 3) collect all method‐specific bed filenames
            beds_to_qc = []
            if args.run_dbscan:
                for eps in (args.dbscan_eps or [select_dbscan_params(df)[0]]):
                    for ms in (args.dbscan_min_samples or [select_dbscan_params(df)[1]]):
                        beds_to_qc.append(f"{prefix}_DBSCAN_eps{eps}_ms{ms}.bed")
            if args.run_hdbscan:
                for mcs in (args.hdbscan_min_cluster_size or []):
                    beds_to_qc.append(f"{prefix}_HDBSCAN_mcs{mcs}_ms{args.hdbscan_min_samples}.bed")
            if args.run_find_peaks:
                for h in (args.peak_height or []):
                    for d in (args.peak_distance or []):
                        beds_to_qc.append(f"{prefix}_findPeaks_h{h}_d{d}.bed")
            if args.run_macs3:
                beds_to_qc.append(f"{prefix}_MACS3_all.bed")
            if args.run_dbscan_sw:
                base_eps = (args.dbscan_eps or [5])[0]
                for b in (args.sw_bin or []):
                    for thr in (args.sw_threshold or []):
                        beds_to_qc.append(f"{prefix}_DBSCAN_SW_eps{base_eps}_bin{b}_th{thr}.bed")
            if args.run_cager:
                beds_to_qc.append(f"{prefix}_CAGEr_clusters.bed")

            # DEBUG: print out the beds to QC and whether they exist
            print(f"[QC][{prefix}] beds_to_qc = {beds_to_qc!r}")
            for bed in beds_to_qc:
                path = os.path.join(sample_dir, bed)
                print(f"[QC][{prefix}]   - {bed}: {'FOUND' if os.path.exists(path) else 'MISSING'}")

            # 4) run QC in parallel
            qc.run_qc_parallel(
                sample_dir=sample_dir,
                beds_to_qc=beds_to_qc,
                prefix=prefix,
                lr_tss=lr_tss,
                annot_tss=annot_tss,
                window=args.qc_dist_cutoff,
                qc_base=qc_base,
                env=args.bedtools_env,
                max_workers=4
            )
            # optionally make plots
            if args.make_plots:
                script = os.path.join(os.path.dirname(__file__), "plot_qc.py")
                print(f"[PLOTS][{prefix}] Generating plots under {qc_base}/plots …")
                subprocess.run(
                    [sys.executable, script, qc_base],
                    check=True
                )
                print(f"[PLOTS][{prefix}] Done.")
        else:
            print(f"[QC][{prefix}] Skipped")


if __name__ == '__main__':
    main()

