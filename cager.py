import os
import subprocess
"""
Run CAGEr via Rscript in the given conda env on one or more BAMs
and write out <sample>_CAGEr_clusters.bed in out_dir.

Parameters
----------
bam_files : str
    Path to a single BAM, or comma-separated BAM paths.
out_dir : str
    Directory to write CAGEr output.
threshold : int
    Minimum tag count for CAGEr (passed to R script).
maxdist : int
    Maximum cluster distance for CAGEr (passed to R script).
rscript_path : str, optional
    Path to the run_cager.R wrapper (defaults to ./run_cager.R next to this file).
env : str
    Name of the Conda environment where CAGEr is installed.

Returns
-------
bed_path : str
    Path to the BED of clusters, e.g. <out_dir>/<sample>_CAGEr_clusters.bed
"""
import subprocess
import os

def run_cager(bam_files, out_dir, threshold, maxdist, env, rscript_path=None):
    """
    ... add rscript_path: optional path to the R wrapper script ...
    """
    # 1) join BAMs
    bam_arg = ",".join(bam_files)
    # 2) derive prefix
    prefix = os.path.basename(bam_files[0]).split(".")[0]

    # 3) locate the R script
    if rscript_path is None:
        # assume we shipped run_cager.R next to this Python file
        module_dir = os.path.dirname(__file__)
        script = os.path.join(module_dir, "run_cager.R")
    else:
        script = rscript_path

    if not os.path.isfile(script):
        raise FileNotFoundError(f"Cannot find R wrapper at {script}")

    # 4) build the command with the full script path
    cmd = [
        "conda", "run", "-n", env,
        "Rscript", script,
        bam_arg, prefix,
        str(threshold), str(maxdist)
    ]

    # 5) make sure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # 6) run *inside* out_dir so your R script's getwd() → out_dir
    subprocess.run(cmd, cwd=out_dir, check=True)

    # 7) verify output
    bed = os.path.join(out_dir, f"{prefix}_CAGEr_clusters.bed")
    if not os.path.exists(bed):
        raise FileNotFoundError(f"CAGEr failed to produce {bed}")
    return bed

