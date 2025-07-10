# === File: clustering.py ===
#!/usr/bin/env python3
"""
Clustering routines: DBSCAN, HDBSCAN, find_peaks, and sliding-window.
"""
import numpy as np
import hdbscan
from sklearn.cluster import DBSCAN
from scipy.signal import find_peaks


def run_dbscan(df, eps, ms):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        pts = grp["pos"].values.reshape(-1,1)
        db = DBSCAN(eps=eps, min_samples=ms, metric="chebyshev").fit(pts)
        grp2 = grp.assign(cluster=db.labels_)
        for cid, sub in grp2[grp2.cluster != -1].groupby("cluster"):
            start, end = int(sub.pos.min()), int(sub.pos.max())+1
            support = int(sub.freq.sum())
            clusters.append({
                "chrom":chrom, "start":start, "end":end,
                "support":support, "strand":strand,
                "method":"DBSCAN", "param":f"eps{eps}_ms{ms}"
            })
    return clusters


def run_hdbscan(df, mcs, ms):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        pts = grp["pos"].values.reshape(-1,1)
        cl = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms, metric="chebyshev").fit(pts)
        grp2 = grp.assign(cluster=cl.labels_)
        for cid, sub in grp2[grp2.cluster != -1].groupby("cluster"):
            start, end = int(sub.pos.min()), int(sub.pos.max())+1
            support = int(sub.freq.sum())
            clusters.append({
                "chrom":chrom, "start":start, "end":end,
                "support":support, "strand":strand,
                "method":"HDBSCAN", "param":f"mcs{mcs}_ms{ms}"
            })
    return clusters


def run_find_peaks(df, height, distance, chrom_lengths):
    clusters = []
    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        length = chrom_lengths.get(chrom, grp.pos.max()+1)
        cov = np.zeros(length, dtype=int)
        for _,r in grp.iterrows(): cov[int(r.pos)] = int(r.freq)
        cov_sm = np.convolve(cov, np.ones(3)/3.0, mode="same")
        peaks, _ = find_peaks(cov_sm, height=height, distance=distance)
        for p in peaks:
            l, r = p, p
            while l>0 and cov_sm[l]>=height/2: l-=1
            while r<length-1 and cov_sm[r]>=height/2: r+=1
            support=int(cov[l:r+1].sum())
            clusters.append({
                "chrom":chrom, "start":l, "end":r+1,
                "support":support, "strand":strand,
                "method":"findPeaks", "param":f"h{height}_d{distance}"
            })
    return clusters


def run_dbscan_sw(df, eps, ms, bin_size, threshold, chrom_lengths):
    from clustering import run_dbscan
    clusters = run_dbscan(df, eps, ms)
    covered = {(cl['chrom'], cl['strand'], p)
               for cl in clusters for p in range(cl['start'], cl['end'])}
    residual = df[~df.apply(lambda r: (r.chrom,r.strand,r.pos) in covered, axis=1)]
    for (chrom, strand), grp in residual.groupby(["chrom","strand"]):
        length = chrom_lengths.get(chrom, grp.pos.max()+1)
        nbin = length//bin_size+1
        hist = np.zeros(nbin, dtype=int)
        for _,r in grp.iterrows(): hist[int(r.pos)//bin_size]+=int(r.freq)
        hot = hist>=threshold
        i=0
        while i<nbin:
            if hot[i]:
                j,supp=i,0
                while j<nbin and hot[j]: supp+=hist[j]; j+=1
                clusters.append({
                    "chrom":chrom, "start":i*bin_size, "end":min(j*bin_size,length),
                    "support":supp, "strand":strand,
                    "method":"DBSCAN_SW", "param":f"eps{eps}_bin{bin_size}_th{threshold}"
                })
                i=j
            else: i+=1
    return clusters