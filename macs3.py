# === File: macs3.py ===
#!/usr/bin/env python3
"""
MACS3 wrapper: calls peak calling per strand using a Conda environment,
and parses output for both narrow and broad peaks.

Parameters
----------
df        : pandas.DataFrame
    DataFrame with columns ["chrom","pos","freq","strand"].
gsize     : str or float
    Genome size string (e.g. "hs" or 1.2e7).
qvalue    : float
    q-value cutoff.
nomodel   : bool, optional
    If True, force a fixed-length model (skip model building).
shift     : int, optional
    If nonzero, pass --shift <shift> (only when building model).
extsize   : int, optional
    If provided, pass --extsize <extsize> when building model.
call_summits : bool, optional
    If True, pass --call-summits (only for narrow peaks with model).
env       : str, optional
    Name of the Conda env where MACS3 is installed.
out_dir   : str, optional
    Directory to write MACS3 outputs (defaults to cwd).

Behavior
--------
- For each chrom/strand, if `nomodel` is False and there are ≥100 tags, build a model
  (narrow mode only), else use `--nomodel --extsize 147` fallback.
- Only narrow mode supports summit calling; broad mode uses fixed extsize.

Returns
-------
clusters : list of dict
    One dict per peak with keys:
    chrom, start, end, support, strand, method, param
"""
import os
import subprocess

# Minimum tag count to attempt model building
MODEL_THRESHOLD = 100

def run_macs3(df,
              gsize,
              qvalue,
              nomodel=False,
              shift=0,
              extsize=None,
              call_summits=False,
              env="macs3-env",
              out_dir=None):
    if out_dir is None:
        out_dir = os.getcwd()
    clusters = []

    def invoke(cmd):
        print(f"[MACS3] Running: {' '.join(cmd)}")
        return subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )

    for (chrom, strand), grp in df.groupby(["chrom", "strand"]):
        tmp_bed = os.path.join(out_dir, f"tmp_{chrom}_{strand}.bed")
        with open(tmp_bed, "w") as f:
            for _, r in grp.iterrows():
                f.write(f"{r.chrom}\t{r.pos}\t{r.pos+1}\t.\t{r.freq}\t{r.strand}\n")

        param_str = f"{chrom}_{strand}_g{gsize}_q{qvalue}"
        for mode in ("narrow", "broad"):
            prefix = f"MACS3_{param_str}_{mode}"
            base_cmd = [
                "conda", "run", "-n", env,
                "macs3", "callpeak",
                "-t", tmp_bed,
                "-f", "BED",
                "-n", prefix,
                "--gsize", str(gsize),
                "--qvalue", str(qvalue),
                "--keep-dup", "all",
                "--outdir", out_dir
            ]
            # broad-specific options
            if mode == "broad":
                base_cmd += ["--broad", "--broad-cutoff", str(qvalue)]

            # decide flags based on model threshold
            do_model = not nomodel and mode == "narrow" and len(grp) >= MODEL_THRESHOLD
            if do_model:
                flags = []
                if shift:
                    flags += ["--shift", str(shift)]
                if extsize:
                    flags += ["--extsize", str(extsize)]
                if call_summits:
                    flags.append("--call-summits")
            else:
                # fallback: fixed model skip
                flags = ["--nomodel", "--extsize", str(extsize or 147)]

            cmd = base_cmd + flags
            try:
                invoke(cmd)
            except subprocess.CalledProcessError as e:
                stderr = e.stderr or ""
                print(f"  [MACS3][{mode}][{strand}] Failed: {stderr.strip()}")
                continue

            # parse output
            peak_file = os.path.join(out_dir, f"{prefix}_peaks.{mode}Peak")
            if os.path.isfile(peak_file):
                with open(peak_file) as pf:
                    for line in pf:
                        cols = line.strip().split("\t")
                        clusters.append({
                            "chrom": cols[0],
                            "start": int(cols[1]),
                            "end":   int(cols[2]),
                            "support": int(float(cols[4])),
                            "strand": strand,
                            "method": f"MACS3_{mode}",
                            "param":   param_str
                        })

        try:
            os.remove(tmp_bed)
        except OSError:
            pass

    return clusters
