# hic-boundary-persistence

Persistence-weighted boundary hierarchy scoring for Hi-C topologically associating domains (TADs). Developed as part of the Fatima Institute AI Research Fellowship, under guidance from Dr. Muhammad Fatima.

## What this is

Standard TAD boundary calling (insulation score) gives a flat yes/no per boundary, with no formal confidence measure, and cannot represent nested or hierarchical domain structure. This project adds a persistent homology-based refinement layer on top of standard breakpoint calling:

1. Stage 1 -- standard insulation-score breakpoint detection (unchanged, using the existing method as implemented in cooltools/cooler).
2. Stage 2 -- windowed persistent homology within each candidate domain, turning boundary strength into a formal persistence value, with a depth-calibrated confidence bound derived from the persistence stability theorem.

See hic_boundary_persistence.pdf for the full methods writeup, including the falsifiable claims being tested and the mathematical derivation of the confidence bound.

## Project status

Early-stage / exploratory. Currently validating the Stage 2 method on synthetic data with known ground truth before applying it to real Hi-C data at scale. See src/synthetic_validation.py for the validation approach and a documented bug fix (diagonal masking) that was necessary to get correct results.

## Repository structure

- src/io_utils.py -- load .hic files via straw, pull contact matrices, plot heatmaps
- src/synthetic_validation.py -- validates Stage 2 persistent homology against known planted boundaries
- src/real_data_persistence.py -- applies the validated method to real GM12878 Hi-C data
- src/insulation.py -- Stage 1: insulation-score breakpoint calling (in progress)
- src/persistence.py -- Stage 2: windowed persistent homology (in progress)
- src/benchmark.py -- benchmark metrics: CTCF/cohesin enrichment, reproducibility (in progress)
- figures/ -- generated plots
- requirements.txt

## Setup

This project requires a Linux environment (native Linux, or WSL2 on Windows) because two core dependencies (cooltools, hic-straw) are C/C++ extensions that fail to build on Windows/MSVC.

### System prerequisites (Ubuntu/WSL2)

    sudo apt update
    sudo apt install -y python3-pip python3-venv build-essential libcurl4-openssl-dev

libcurl4-openssl-dev is required specifically for hic-straw to build (it streams .hic files over HTTP).

### Python environment

    python3 -m venv venv
    source venv/bin/activate
    pip install -r requirements.txt

## Running

Pull a real contact matrix and view it:

    python3 src/io_utils.py

Validate the persistent homology method against synthetic data with known planted boundaries:

    python3 src/synthetic_validation.py

Run the validated method against real GM12878 Hi-C data:

    python3 src/real_data_persistence.py

Output plots land in figures/.

## Data

Uses public Hi-C data (Rao et al. 2014, GM12878, GEO GSE63525) via direct HTTP streaming -- no local download or data redistribution. ENCODE intact Hi-C data is also being explored for higher-resolution analysis (see project notes).

## Compute

Genome-wide runs are intended to run on Modal (https://modal.com), a serverless compute platform, once the method is validated locally. This is not yet implemented.
