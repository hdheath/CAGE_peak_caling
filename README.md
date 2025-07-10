# Peak Calling on CAGE Data

This repository provides a streamlined pipeline for strand-aware peak calling on CAGE (Cap Analysis Gene Expression) data and subsequent QC. Intended for use by the Brook's Lab.

---

## Project Layout

```plaintext
CAGE_peak_calling/
├── config.toml              # Pipeline configuration (inputs, methods, params)
├── main.py                  # Driver script (dry-run & real run)
├── utils.py                 # Argument parsing and core utilities
├── macs3.py                 # MACS3 wrapper
├── cager.py                 # CAGEr wrapper
├── clustering/              # DBSCAN, HDBSCAN, find_peaks, sliding-window modules
├── qc/                      # QC functions (TSS extraction, bedtools calls)
├── plot_qc.py               # Script to generate QC plots
└── output/                  # Results (one subfolder per sample)
```

Data files (BAMs, etc.) live outside the repo; see `config.toml` for their paths.

---

## Installation

1. Clone the repo:

   ```bash
   git clone git@github.com:YOUR_USERNAME/CAGE_peak_calling.git
   cd CAGE_peak_calling
   ```

2. Create and activate environments:

   ```bash
   # for peak calling
   conda env create -f peak_calling_env.yaml  # defines pysam, sklearn, scipy, etc.
   conda activate peak_calling_env

   # for CAGEr (R-based)
   conda create -n cager-env r-base r-cager

   # for bedtools (QC)
   conda create -n bedtools_env bedtools
   ```

---

## Configuration (`config.toml`)

All inputs, methods, and parameter settings are driven by `config.toml`.

Key sections:

* **Required inputs**: `bams`, `out_dir`
* **Peak-calling methods**: flags `run_dbscan`, `run_hdbscan`, `run_find_peaks`, `run_macs3`, `run_cager`, `run_dbscan_sw` and their parameter arrays (e.g. `dbscan_eps = [5,10]`).
* **QC**: `long_reads`, `annotation_bed`, `qc_dist_cutoff`, `make_plots`.

Example:

```toml
run_macs3    = true
macs3_gsize  = ["hs"]
macs3_qvalue = [0.01, 0.05]
```

---

## Usage

### Dry-Run Mode (preview actions)

```bash
python main.py --config config.toml --dry_run
```

### Real Run

```bash
python main.py --config config.toml
```

Outputs are written to `output/<sample>/` as BED files.

---

## Methods Overview

### 1. DBSCAN

* **What it does**: Groups TSS 5′-end coverage points into clusters based on density.
* **Key parameters**:

  * `dbscan_eps`: maximum distance between two points to be considered in the same neighborhood (e.g. `[5, 10]` bp).
  * `dbscan_min_samples`: minimum number of points to form a cluster (e.g. `[4, 6]`).
* **Config keys**: `run_dbscan`, `dbscan_eps`, `dbscan_min_samples`.

### 2. DBSCAN + Sliding Window

* **What it does**: Aggregates reads in sliding bins, then calls DBSCAN on bin centers exceeding thresholds.
* **Key parameters**:

  * `sw_bin`: bin sizes in bp (e.g. `[5, 10]`).
  * `sw_threshold`: minimum read count per bin (e.g. `[4, 6]`).
* **Config keys**: `run_dbscan_sw`, `sw_bin`, `sw_threshold`, `dbscan_eps`, `dbscan_min_samples`.

### 3. SciPy `find_peaks`

* **What it does**: Identifies local maxima in per-position coverage profiles.
* **Key parameters**:

  * `peak_height`: minimum read count at a peak (e.g. `[4, 6]`).
  * `peak_distance`: minimum distance between peaks (e.g. `[5, 10]`).
* **Config keys**: `run_find_peaks`, `peak_height`, `peak_distance`.

### 4. MACS3

* **What it does**: Model-based peak calling widely used for ChIP-Seq, repurposed here for CAGE.
* **Key parameters**:

  * `macs3_gsize`: genome size code (`"hs"`, `"mm"`, etc.).
  * `macs3_qvalue`: false discovery rate cutoff (e.g. `[0.01, 0.05]`).
* **Config keys**: `run_macs3`, `macs3_gsize`, `macs3_qvalue`.

### 5. CAGEr (R)

* **What it does**: Clusters tag start sites into CAGE-defined peaks.
* **Key parameters**:

  * `cager_threshold`: minimum tag count per cluster (e.g. `[4]`).
  * `cager_maxdist`: maximum distance to merge tags (e.g. `[10]` bp).
* **Config keys**: `run_cager`, `cager_threshold`, `cager_maxdist`, `cager_env`.

---

## Quality Control (QC)

1. **Extract TSS**

   * From long-read BAM(s) and annotation BED.

2. **Bedtools `closest`**

   * Compare CAGE peaks to long-read TSS and annotation TSS.
   * Metrics:

     * Precision/Recall of long reads within peaks.
     * Distance distributions (% within `qc_dist_cutoff` bp).

3. **Region Length Distribution**

   * Compute lengths of each called region across methods.

4. **Summarization**

   * Boxplots of peak lengths by method.
   * Bar plots for precision, recall, and distance metrics.


Plots will appear under `output/QC/plots/`.

---

## Placeholders for Figures

* **Peak-length distribution**
  ![Peak-length boxplot placeholder](path/to/lengths_boxplot.png)

* **Precision/Recall barplot**
  ![Precision/Recall placeholder](path/to/pr_barplot.png)

