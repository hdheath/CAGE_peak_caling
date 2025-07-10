# === File: macs3.py ===
#!/usr/bin/env python3
"""
MACS3 wrapper: calls peak calling per strand using conda environment,
and parses output for both narrow and broad peaks.
"""

import os
import subprocess

def run_macs3(df, gsize, qvalue, out_dir):
    clusters = []

    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        tmp_bed = os.path.join(out_dir, f"tmp_{chrom}_{strand}.bed")
        with open(tmp_bed, "w") as f:
            for _, r in grp.iterrows():
                f.write(f"{r.chrom}\t{r.pos}\t{r.pos + 1}\t.\t{r.freq}\t{r.strand}\n")

        param_str = f"{chrom}_{strand}_g{gsize}_q{qvalue}"
        for mode in ["narrow", "broad"]:
            prefix = f"MACS3_{param_str}_{mode}"
            cmd = [
                "conda", "run", "-n", "macs3-env", "macs3", "callpeak",
                "-f", "BED", "-t", tmp_bed, "-n", prefix,
                "--gsize", str(gsize), "--qvalue", str(qvalue),
                "--nomodel", "--keep-dup", "all", "--outdir", out_dir
            ]
            if mode == "broad":
                cmd += ["--broad", "--broad-cutoff", str(qvalue)]

            print(f"[MACS3][{mode}][{strand}] Running: {' '.join(cmd)}")
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
            except subprocess.CalledProcessError as e:
                print(f"  MACS3 {mode} run failed for {strand} strand:\n{e.stderr.strip()}")
                continue

            peak_file = os.path.join(out_dir, f"{prefix}_peaks.{mode}Peak")
            if os.path.isfile(peak_file):
                with open(peak_file) as pf:
                    for line in pf:
                        cols = line.strip().split("\t")
                        clusters.append({
                            "chrom": cols[0],
                            "start": int(cols[1]),
                            "end": int(cols[2]),
                            "support": int(float(cols[4])),
                            "strand": strand,
                            "method": f"MACS3_{mode}",
                            "param": param_str
                        })

        os.remove(tmp_bed)

    return clusters
