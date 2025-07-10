#!/usr/bin/env python3
"""
plot_qc.py

Plot QC results from the TSS peak-calling QC pipeline.
Generates and saves the following plots under <qc_dir>/plots/<sample>/:

For each sample:
1) Peak-length distributions (box + violin overlaid)
2) Bar plot of number of peaks per method
3) Multi-panel histograms of read counts per peak window (one panel per method)
4) Distance-to-read-TSS distributions for reads not in windows (box + violin)
5) Distance-to-annot-TSS distributions for annotation TSSs not in windows (box + violin)
6) Grouped bar plot of annotation TSS precision & recall per method
7) Bar plot of read TSS recall per method
"""
import os
import glob
import argparse
import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    return path


def load_lengths(qc_dir):
    print(f"[LOAD] Looking for length files in {qc_dir}")
    pattern = os.path.join(qc_dir, '*_lengths.tsv')
    files = sorted(glob.glob(pattern))
    data = {}
    for f in files:
        base = os.path.basename(f).replace('_lengths.tsv', '')
        sample, method = base.split('_', 1)
        df = pd.read_csv(f, sep='\t')
        data.setdefault(sample, {})[method] = df['length'].values
        print(f"[LOAD] {len(df)} lengths loaded for {sample}/{method}")
    return data


def load_read_counts(qc_dir):
    print(f"[LOAD] Looking for read-count files in {qc_dir}")
    pattern = os.path.join(qc_dir, '*_read_counts_per_window.tsv')
    files = sorted(glob.glob(pattern))
    data = {}
    for f in files:
        base = os.path.basename(f).replace('_read_counts_per_window.tsv', '')
        sample, method = base.split('_', 1)
        df = pd.read_csv(f, sep='\t')
        if 'ReadCount' not in df.columns:
            print(f"[WARN] ReadCount missing in {f}, skipping")
            continue
        data.setdefault(sample, {})[method] = df['ReadCount'].values
        print(f"[LOAD] {len(df)} read counts loaded for {sample}/{method}")
    return data


def load_distances(qc_dir, suffix):
    print(f"[LOAD] Looking for distance files ({suffix}) in {qc_dir}")
    pattern = os.path.join(qc_dir, f'*_{suffix}_distances.tsv')
    files = sorted(glob.glob(pattern))
    data = {}
    for f in files:
        base = os.path.basename(f).replace(f'_{suffix}_distances.tsv', '')
        sample, method = base.split('_', 1)
        df = pd.read_csv(f, sep='\t')
        if 'Distance' not in df.columns:
            print(f"[WARN] Distance missing in {f}, skipping")
            continue
        data.setdefault(sample, {})[method] = df['Distance'].values
        print(f"[LOAD] {len(df)} distances loaded for {sample}/{method}")
    return data


def load_metrics(qc_dir):
    print(f"[LOAD] Looking for PR metrics in {qc_dir}")
    pattern = os.path.join(qc_dir, '*_PR_metrics.tsv')
    files = sorted(glob.glob(pattern))
    records = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        records.append(df)
        print(f"[LOAD] {len(df)} PR records loaded from {f}")
    if records:
        combined = pd.concat(records, ignore_index=True)
        print(f"[LOAD] Combined PR metrics: {len(combined)} records total")
        return combined
    else:
        print("[LOAD] No PR metrics files found")
        return pd.DataFrame()


def plot_length_distribution(methods_dict, out_dir):
    methods = sorted(methods_dict.keys())
    print(f"[PLOT] Peak length distribution for methods: {methods}")
    values = [methods_dict[m] for m in methods]
    pos = np.arange(len(methods)) + 1

    fig, ax = plt.subplots(figsize=(max(6, len(methods)), 6))
    ax.violinplot(values, positions=pos, showextrema=False)
    ax.boxplot(values, positions=pos, widths=0.1)
    ax.set_xticks(pos)
    ax.set_xticklabels(methods, rotation=45, ha='right')
    ax.set_ylabel('Peak length')
    ax.set_title('Peak Length Distribution')
    fig.tight_layout()
    path = os.path.join(out_dir, 'peak_length_distribution.png')
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def plot_peak_counts(methods_dict, out_dir):
    methods = sorted(methods_dict.keys())
    print(f"[PLOT] Number of peaks per method: {methods}")
    counts = [len(methods_dict[m]) for m in methods]
    fig, ax = plt.subplots(figsize=(max(6, len(methods)), 4))
    ax.bar(methods, counts)
    ax.set_ylabel('Number of peaks')
    ax.set_title('Number of Peaks per Method')
    ax.set_xticklabels(methods, rotation=45, ha='right')
    fig.tight_layout()
    path = os.path.join(out_dir, 'peak_counts_bar.png')
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def plot_multi_panel_histograms(methods_dict, out_dir, cols=2, bins=50):
    methods = sorted(methods_dict.keys())
    n = len(methods)
    rows = math.ceil(n / cols)
    print(f"[PLOT] Multi-panel histograms for methods: {methods} (rows={rows}, cols={cols})")
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), squeeze=False)

    for idx, method in enumerate(methods):
        r, c = divmod(idx, cols)
        ax = axes[r][c]
        values = np.clip(methods_dict[method], 0, 500)
        ax.hist(values, bins=bins)
        ax.set_title(method, fontsize='small')
        ax.set_xlabel('Read counts (capped at 500)')
        ax.set_ylabel('Frequency')

    for idx in range(n, rows * cols):
        r, c = divmod(idx, cols)
        axes[r][c].axis('off')

    fig.tight_layout()
    path = os.path.join(out_dir, 'read_counts_histograms.png')
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def plot_distances(methods_dict, out_dir, title, filename):
    methods = sorted(methods_dict.keys())
    print(f"[PLOT] {title} for methods: {methods}")
    values = [methods_dict[m] for m in methods]
    pos = np.arange(len(methods)) + 1

    fig, ax = plt.subplots(figsize=(max(6, len(methods)), 6))
    ax.violinplot(values, positions=pos, showextrema=False)
    ax.boxplot(values, positions=pos, widths=0.1)
    ax.set_xticks(pos)
    ax.set_xticklabels(methods, rotation=45, ha='right')
    ax.set_ylabel('Distance to nearest window')
    ax.set_title(title)
    fig.tight_layout()
    path = os.path.join(out_dir, filename)
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def plot_annotation_pr(metrics_df, out_dir):
    df = metrics_df.copy()
    df['method'] = df['method'].str.replace('_PR$', '', regex=True)
    methods = sorted(df['method'])
    print(f"[PLOT] Annotation precision & recall for methods: {methods}")
    precision = df['annotation_precision'].values
    recall = df['annotation_recall'].values

    x = np.arange(len(methods))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, len(methods)), 5))
    ax.bar(x - width/2, precision, width, label='Precision')
    ax.bar(x + width/2, recall, width, label='Recall')
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha='right')
    ax.set_ylabel('Value')
    ax.set_title('Annotation TSS Precision & Recall per Method')
    ax.legend()
    fig.tight_layout()
    path = os.path.join(out_dir, 'annotation_pr_bar.png')
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def plot_read_recall(metrics_df, out_dir):
    df = metrics_df.copy()
    df['method'] = df['method'].str.replace('_PR$', '', regex=True)
    methods = sorted(df['method'])
    print(f"[PLOT] Read recall for methods: {methods}")
    recall = df['read_recall'].values

    fig, ax = plt.subplots(figsize=(max(6, len(methods)), 4))
    ax.bar(methods, recall)
    ax.set_ylabel('Read Recall')
    ax.set_title('Read TSS Recall per Method')
    ax.set_xticklabels(methods, rotation=45, ha='right')
    fig.tight_layout()
    path = os.path.join(out_dir, 'read_recall_bar.png')
    fig.savefig(path)
    plt.close(fig)
    print(f"[SAVE] -> {path}")


def main():
    parser = argparse.ArgumentParser(description='Plot QC results per sample')
    parser.add_argument('qc_dir', help='Directory containing QC TSV outputs')
    parser.add_argument('--plots_dir', help='Base directory to save plots (default: <qc_dir>/plots)')
    parser.add_argument('--hist_cols', type=int, default=2, help='Number of columns for multi-panel histograms')
    args = parser.parse_args()

    qc_dir = args.qc_dir
    base_plots = args.plots_dir or os.path.join(qc_dir, 'plots')
    ensure_dir(base_plots)

    print(f"[START] plot_qc.py invoked with qc_dir={qc_dir}")

    lengths_all = load_lengths(qc_dir)
    counts_all = load_read_counts(qc_dir)
    reads_dist_all = load_distances(qc_dir, 'readsNotInWindows')
    annot_dist_all = load_distances(qc_dir, 'annotNotInWindows')
    metrics_all = load_metrics(qc_dir)

    samples = sorted(set(lengths_all) | set(counts_all) | set(reads_dist_all) | set(annot_dist_all) | set(metrics_all.get('sample', [])))

    for sample in samples:
        print(f"[SAMPLE] Generating plots for {sample}")
        sample_dir = os.path.join(base_plots, sample)
        ensure_dir(sample_dir)

        if sample in lengths_all:
            plot_length_distribution(lengths_all[sample], sample_dir)
            plot_peak_counts(lengths_all[sample], sample_dir)
        if sample in counts_all:
            plot_multi_panel_histograms(counts_all[sample], sample_dir, cols=args.hist_cols)
        if sample in reads_dist_all:
            plot_distances(reads_dist_all[sample], sample_dir,
                           'Distance Distribution: Reads not in Windows',
                           'reads_notinwindows_dist.png')
        if sample in annot_dist_all:
            plot_distances(annot_dist_all[sample], sample_dir,
                           'Distance Distribution: Annot. TSS not in Windows',
                           'annot_notinwindows_dist.png')
        if not metrics_all.empty:
            df_sample = metrics_all[metrics_all['sample'] == sample]
            if not df_sample.empty:
                plot_annotation_pr(df_sample, sample_dir)
                plot_read_recall(df_sample, sample_dir)

    print(f"[DONE] All plots generated under {base_plots}")

if __name__ == '__main__':
    main()
