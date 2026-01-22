# SeekRNA GA Project

This is a modularized version of the SeekRNA Genetic Algorithm for RNA sequence evolution.

## Project Structure

- `main.py`: Entry point for the GA.
- `config.json`: Configuration for sequences, target sites, and GA parameters.
- `config_loader.py`: Shorthand utility to load configurations.
- `core/`:
    - `ga_engine.py`: The main GA loop and parallel processing.
    - `operators.py`: Mutation and crossover logic.
    - `scoring.py`: Tier 1 and Tier 2 scoring functions.
    - `sampling.py`: Weighted selection methods.
- `utils/`:
    - `sequence.py`: Basic RNA sequence utilities.
    - `coords.py`: Coordinate mapping and update logic.
    - `folding.py`: External folding tool wrappers (RNAfold, mxfold2, centroid_fold).

## Prerequisites

Ensure the following tools are installed and available in your `PATH`:
- `RNAfold` (ViennaRNA Package)
- `mxfold2`
- `centroid_fold`

## How to Run

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Place your initial candidates CSV file (`initial_candidates_ISEc11_ISPa11.csv`) in the project root.

3. Run the GA:
   ```bash
   python main.py
   ```

