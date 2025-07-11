# === File: cager.py ===
#!/usr/bin/env python3
"""
CAGEr wrapper: runs an Rscript in a specified Conda environment,
and returns the path to the generated clusters BED.

Parameters
----------
bam_files       : list of str
    One or more BAM file paths.
out_dir         : str
    Directory to write CAGEr output.
threshold       : int
    Minimum tag count for clusters.
maxdist         : int
    Maximum clustering distance.
correctFirstG   : bool, optional
    If True, pass --correctFirstG TRUE to the R script.
env             : str
    Name of the Conda environment where CAGEr is installed.
rscript_path    : str, optional
    Full path to your run_cager.R wrapper (defaults to sibling file).
Returns
-------
bed_path        : str
    Path to the <prefix>_CAGEr_clusters.bed produced in out_dir.
"""
import os
import subprocess

def run_cager(bam_files,
              out_dir,
              threshold,
              maxdist,
              correctFirstG=False,
              env="cager-env",
              rscript_path=None):
    # ensure list
    if isinstance(bam_files, str):
        bam_files = [bam_files]
    bam_arg = ",".join(bam_files)
    prefix = os.path.basename(bam_files[0]).split(".")[0]

    # locate R wrapper
    if rscript_path is None:
        module_dir = os.path.dirname(__file__)
        script = os.path.join(module_dir, "run_cager.R")
    else:
        script = rscript_path
    if not os.path.isfile(script):
        raise FileNotFoundError(f"CAGEr wrapper not found at {script}")

    # build command
    cmd = [
        "conda", "run", "-n", env,
        "Rscript", script,
        bam_arg, prefix,
        str(threshold), str(maxdist)
    ]
    if correctFirstG:
        cmd += ["--correctFirstG", "TRUE"]

    # run inside out_dir
    os.makedirs(out_dir, exist_ok=True)
    print(f"[CAGEr] Running: {' '.join(cmd)} (cwd={out_dir})")
    subprocess.run(cmd, cwd=out_dir, check=True)

    # verify output
    bed = os.path.join(out_dir, f"{prefix}_CAGEr_clusters.bed")
    if not os.path.exists(bed):
        raise FileNotFoundError(f"CAGEr did not produce expected BED: {bed}")
    return bed


