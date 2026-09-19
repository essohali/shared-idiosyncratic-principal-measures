# Allen Neuropixels application

This directory contains the frozen derived data and code used for Section 8 of the manuscript.

## Statistical reproduction

From the repository root:

```bash
Rscript application/R/run_application.R
```

No NWB download is required. Results are written to `results/application/`.

## Full data reconstruction

To reconstruct the derived dataset from DANDI:

```bash
python application/python/stream_allen_to_excel_final_corrected.py
```

Dandiset 000022, version 0.251116.2247, session 794812542. The analysis uses 75 presentations at orientation 0 degrees, contrast 0.8, temporal frequency 2 Hz, and spatial frequency 0.04 cycles/degree, with a [0,2] s analysis window.
