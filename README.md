# paleo-polarity-standardizer
An automated, iterative Python tool for standardizing paleomagnetic Virtual Geomagnetic Pole (VGP) coordinates to a single canonical hemisphere.

This code provides a systematic data-cleaning step before building Apparent Polar Wander Paths (APWPs), computing running means, or exporting plate motion files (.rot) for plate reconstruction software such as GPlates.

# What This Code Does
In paleomagnetic literature and global databases (e.g., PINT, PDB, or GPMDB), pole positions are recorded as reported by the original authors—sometimes in the Northern Hemisphere, sometimes in the Southern. Before executing time-series analysis or constructing spatial means, these directions must be mapped into a unified hemisphere convention (typically Geocentric Axial Dipole normal polarity in the Northern Hemisphere).

This script:

Takes raw VGP coordinates (Plat, Plon), age constraints (Age in Ma), and tectonic descriptors (Craton).
Standardizes polarity using a two-phase iterative algorithm anchored by local temporal-spatial peer groups.
Outputs antipodally converted coordinates (Plat, Plon) with explicit status classifications (PTest).

# What This Code Does NOT Do

This script is NOT a statistical reversal test.
The McFadden & McElhinny (1990) reversal test is a hypothesis test evaluating whether pre-separated Normal and Reversed directional populations are statistically antipodal by testing concentrations and critical angles. This code does not calculate angular precision parameters (kappa), confidence ellipses (alpha95), or reversal grades. Instead, it performs the spatial-temporal coordinate unification that precedes such testing.

# How the Algorithm Works
# Phase 1: Canonical-Anchored Global Seed
Standard iterative mean calculations can suffer from majority-vote bias. For example, if a dataset contains 80% reversed poles, a simple iterative Fisher mean will converge toward the reversed hemisphere rather than the canonical convention.

Fisher Mean Vector: The script calculates the unweighted unit vector sum across all valid poles.

Canonical Anchor: If the initial seed mean latitude falls outside the selected canonical hemisphere (e.g., latitude < 0° when CANONICAL_HEMISPHERE = 'N'), the seed mean is inverted before checking distances.

Iterative Inversion: Any pole at an angular distance > 90° from the seed mean is inverted. The mean is recalculated until the flip-set stabilizes.

# Phase 2: Local Coeval Window Refinement
Because tectonic plates undergo significant polar wander over hundreds of millions of years, a single global mean becomes invalid for long-duration datasets. Phase 2 uses a dynamic, localized moving-window reference:

Coeval Peer Grouping: For each pole, the algorithm selects coeval peer poles within a temporal window (+/- 15 Ma) belonging to the same cratonic family.

Local Reference Mean: If 3 or more peers exist, the local reference direction is calculated from the peers. Otherwise, it falls back to the global mean.

RAW Coordinate Comparison: The script evaluates the pole's RAW (unmodified) coordinates against the reference mean. This prevents double-flipping artifacts across iterations.

State-History Cycle Detection: The script tracks the full state history of flip-sets. If the algorithm enters an oscillation cycle, iteration breaks and oscillating poles are flagged as PTest = 'U' (Unstable).

Vector Singularity Guard: If the resultant vector length R < 0.01 * N (a 50/50 N/R split), the script defaults to the mean of the canonical sub-cluster.

# Input & Output Formats
Input CSV Requirements
Input files must be UTF-8 formatted CSVs containing at least these column headers:
Plat: Pole latitude in degrees (-90 to +90)
Plon: Pole longitude in degrees (0 to 360 or -180 to +180)
Age: Mean age of the rock unit in Ma
Craton: Tectonic identifier (e.g., Baltica, Laurentia)

Output Classifications (PTest Column)
The script appends or updates a PTest column in the output file:
PTest = 'N': Normal polarity. Coordinates remain identical to raw input.
PTest = 'R': Reversed polarity. Coordinates flipped to canonical hemisphere.
PTest = 'U': Unstable. The pole oscillated between states across iterations (flagged for manual review).

# How to Run
Place PolarityStandardizer_v4.py in the same folder as your input CSV files.

Open terminal or command prompt and run:
python PolarityStandardizer_v4.py

Output files will be generated with the prefix PolarityStd_ in the same folder.

# Scientific References
If you use this script in your workflow, please cite the foundational methods below:

Fisher Spherical Statistics:
Fisher, R. A. (1953). Dispersion on a sphere. Proc. R. Soc. Lond. A, 217, 295-305.

Reversal Test Concepts:
McFadden, P. L., & McElhinny, M. W. (1990). Classification of the reversal test in palaeomagnetism. Geophys. J. Int., 103, 725-729.

Apparent Polar Wander & Coeval Windowing:
Torsvik, T. H., et al. (2012). Phanerozoic polar wander, palaeogeography and dynamics. Earth-Science Reviews, 114, 325-368.

# License
This project is open-source and distributed under the MIT License. You are free to modify, distribute, and integrate this software into academic workflows provided original attribution is retained.
