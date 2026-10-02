"""
Polarity Standardizer v4
========================
Hemisphere-unification of palaeomagnetic VGPs (virtual geomagnetic poles).

Reads CSV files containing paleopole coordinates (Plat, Plon) and assigns
each pole a polarity test label (PTest):

    N  — Normal   : raw direction consistent with coeval peers; coords unchanged.
    R  — Reversed : raw direction antipodal to coeval peers; Plat/Plon flipped.
    U  — Unstable : pole oscillates between N and R across iterations.
    '' — Skipped  : Age or coordinates could not be parsed.

Output is written as  PolarityStd_<input_filename>.csv  in the same directory
as the input file.

Algorithm
---------
Phase 1 (global seed):
    Iterative Fisher-mean classification anchored to the canonical hemisphere.
    Eliminates majority-vote bias by forcing the seed mean into the canonical
    hemisphere before iteration begins.

Phase 2 (local refinement):
    Each pole is reclassified against the Fisher mean of coeval peers from the
    same craton family (±WINDOW_MA Ma, ≥MIN_PEERS peers).  Full state-history
    cycle detection (frozenset) prevents infinite loops.  Oscillating poles are
    marked PTest='U'.

Usage
-----
    python PolarityStandardizer_v4.py                        # scan current directory
    python PolarityStandardizer_v4.py path/to/folder         # scan a specific folder
    python PolarityStandardizer_v4.py file1.csv file2.csv    # explicit files

Required CSV columns: Plat, Plon, Age, Craton

Author : pateljitendran@gmail.com
License: MIT
Version: 4.0.0
"""

import os
import sys
import math
import glob
import csv
import warnings
from collections import defaultdict


# ---------------------------------------------------------------------------
# CONFIGURATION
# Adjust these constants to tune the algorithm for your dataset.
# ---------------------------------------------------------------------------

# Canonical hemisphere for output: 'N' (Plat > 0) or 'S' (Plat < 0).
# All output poles are placed in this hemisphere.
# 'N' is the standard palaeomagnetic convention (GAD normal polarity = northern).
CANONICAL_HEMISPHERE = 'N'

WINDOW_MA     = 15    # ±Ma for local coeval peer window
MIN_PEERS     = 3     # minimum coeval peers required to use local mean
MAX_ITER      = 50    # safety cap on iteration count (Phases 1 and 2)
SINGULARITY_R = 0.01  # resultant vector length below which mean is indeterminate


# ---------------------------------------------------------------------------
# SPHERICAL TRIGONOMETRY & FISHER STATISTICS HELPERS
# ---------------------------------------------------------------------------

def _rad(deg):
    return deg * math.pi / 180.0


def _deg(rad):
    return rad * 180.0 / math.pi


def pole_to_xyz(lat_deg, lon_deg):
    """Convert (lat, lon) in degrees to a unit Cartesian vector."""
    la = _rad(lat_deg)
    lo = _rad(lon_deg)
    return (
        math.cos(la) * math.cos(lo),
        math.cos(la) * math.sin(lo),
        math.sin(la),
    )


def xyz_to_pole(x, y, z):
    """
    Convert a Cartesian vector to (lat, lon) in degrees.

    Returns (None, None) if the vector has near-zero length.
    """
    r = math.sqrt(x * x + y * y + z * z)
    if r < 1e-12:
        return None, None
    x /= r
    y /= r
    z /= r
    return _deg(math.asin(max(-1.0, min(1.0, z)))), _deg(math.atan2(y, x))


def vector_mean(lat_list, lon_list):
    """
    Compute the Fisher vector mean of a list of poles.

    Parameters
    ----------
    lat_list : list[float]
    lon_list : list[float]

    Returns
    -------
    mean_lat : float or None
    mean_lon : float or None
    R        : float — resultant vector length (0 … N)
    """
    sx = sy = sz = 0.0
    for la, lo in zip(lat_list, lon_list):
        x, y, z = pole_to_xyz(la, lo)
        sx += x
        sy += y
        sz += z
    R = math.sqrt(sx * sx + sy * sy + sz * sz)
    la, lo = xyz_to_pole(sx, sy, sz)
    return la, lo, R


def angular_distance(lat1, lon1, lat2, lon2):
    """Return the great-circle angular distance in degrees between two poles."""
    p1, l1 = _rad(lat1), _rad(lon1)
    p2, l2 = _rad(lat2), _rad(lon2)
    cos_g = (
        math.sin(p1) * math.sin(p2)
        + math.cos(p1) * math.cos(p2) * math.cos(l1 - l2)
    )
    return _deg(math.acos(max(-1.0, min(1.0, cos_g))))


def flip(lat, lon):
    """Return the antipodal pole: (-lat, lon+180°)."""
    return -lat, (lon + 180.0) % 360.0


def is_canonical(lat, hemisphere):
    """Return True if lat already lies in the canonical hemisphere."""
    if hemisphere == 'N':
        return lat >= 0.0
    return lat <= 0.0


def force_canonical(lat, lon, hemisphere):
    """Flip (lat, lon) into the canonical hemisphere if necessary."""
    if not is_canonical(lat, hemisphere):
        lat, lon = flip(lat, lon)
    return lat, lon


def nearest_sub_cluster(lat_list, lon_list, canonical_hemisphere):
    """
    Fallback mean for degenerate (near-zero R) vector sums.

    Splits poles into two antipodal sub-clusters by angular proximity to a
    canonical reference direction, then returns the mean of the sub-cluster
    that lies in the canonical hemisphere.

    Parameters
    ----------
    lat_list            : list[float]
    lon_list            : list[float]
    canonical_hemisphere: str  — 'N' or 'S'

    Returns
    -------
    (mean_lat, mean_lon) : tuple[float, float]
        Falls back to (±45°, 0°) if both sub-clusters are empty.
    """
    canonical_ref = (45.0, 0.0) if canonical_hemisphere == 'N' else (-45.0, 0.0)
    near, far_p = [], []
    for la, lo in zip(lat_list, lon_list):
        ang = angular_distance(la, lo, *canonical_ref)
        (near if ang <= 90.0 else far_p).append((la, lo))

    for group in [near, far_p]:
        if group:
            lats = [p[0] for p in group]
            lons = [p[1] for p in group]
            mla, mlo, R = vector_mean(lats, lons)
            if mla is not None:
                return mla, mlo

    return (45.0, 0.0) if canonical_hemisphere == 'N' else (-45.0, 0.0)


# ---------------------------------------------------------------------------
# DATA INGESTION & FORMATTING HELPERS
# ---------------------------------------------------------------------------

REQUIRED_COLS = {"Plat", "Plon", "Age", "Craton"}


def read_csv_file(filepath):
    """
    Read a CSV file and return (fieldnames, rows).

    Raises
    ------
    ValueError
        If any column from REQUIRED_COLS is absent.
    """
    with open(filepath, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        fieldnames = [f.strip() for f in (reader.fieldnames or [])]
        reader.fieldnames = fieldnames
        missing = REQUIRED_COLS - set(fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns: {sorted(missing)}  |  Found: {fieldnames}"
            )
        return fieldnames, list(reader)


def write_csv_file(filepath, fieldnames, rows):
    """Write rows to a CSV file with the given field order."""
    with open(filepath, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def craton_family(s):
    """
    Extract the craton family name from a sub-craton string.

    Example: 'Baltica-Karelia' -> 'Baltica'
    """
    if not s or str(s).strip().lower() in ('', 'nan'):
        return None
    return str(s).split('-')[0].strip()


def parse_coords(row):
    """Return (plat, plon) as floats, or raise ValueError."""
    return float(row['Plat']), float(row['Plon'])


def parse_age(row):
    """Return Age as float, or None if unparseable."""
    try:
        return float(row['Age'])
    except (ValueError, TypeError):
        return None


def preserve_fmt(original_str, new_value):
    """
    Format new_value to the same number of decimal places as original_str.

    Ensures the output CSV retains the same numeric precision as the input.
    """
    s = str(original_str).strip()
    if '.' in s:
        dec = len(s.split('.')[-1])
        return f"{new_value:.{dec}f}"
    return str(int(round(new_value)))


# ---------------------------------------------------------------------------
# VECTOR SINGULARITY-GUARDED MEAN
# ---------------------------------------------------------------------------

def safe_mean(lats, lons, canonical_hemisphere, context=''):
    """
    Compute Fisher vector mean with a singularity guard.

    If the resultant vector length R < SINGULARITY_R * N (indicating
    near-equal antipodal populations), falls back to nearest_sub_cluster
    and emits a RuntimeWarning.

    Parameters
    ----------
    lats                : list[float]
    lons                : list[float]
    canonical_hemisphere: str   — 'N' or 'S'
    context             : str   — label included in the warning message

    Returns
    -------
    (mean_lat, mean_lon) : tuple[float, float]
    """
    if not lats:
        return force_canonical(45.0, 0.0, canonical_hemisphere)

    mla, mlo, R = vector_mean(lats, lons)

    if mla is None or R < SINGULARITY_R * len(lats):
        warnings.warn(
            f"Vector sum near zero (R={R:.4f}) {context}. "
            "Falling back to canonical sub-cluster mean. "
            "This window may contain a true 50/50 mix of N and R poles.",
            RuntimeWarning,
            stacklevel=2,
        )
        mla, mlo = nearest_sub_cluster(lats, lons, canonical_hemisphere)

    return mla, mlo


# ---------------------------------------------------------------------------
# PHASE 1 — GLOBAL SEED WITH CANONICAL-HEMISPHERE ANCHORING
# ---------------------------------------------------------------------------

def phase1_global(raw_coords, valid_idx, canonical_hemisphere):
    """
    Iterative global Fisher-mean seed, anchored to the canonical hemisphere.

    The canonical anchor eliminates majority-vote bias: if the raw data
    has 80% reversed poles, the raw global mean lands in the reversed
    hemisphere.  We detect this and flip the seed mean to the canonical
    hemisphere before starting iteration, ensuring convergence to the
    canonical orientation regardless of the N/R ratio in the input.

    Parameters
    ----------
    raw_coords           : list[tuple[float, float] or None]
    valid_idx            : list[int]
    canonical_hemisphere : str — 'N' or 'S'

    Returns
    -------
    flip_set    : set[int]        — row indices of poles to flip
    global_mean : tuple[float, float] — stable mean in canonical hemisphere
    iterations  : int
    """
    all_lats = [raw_coords[i][0] for i in valid_idx]
    all_lons = [raw_coords[i][1] for i in valid_idx]

    mean_lat, mean_lon = safe_mean(
        all_lats, all_lons, canonical_hemisphere,
        context='Phase-1 initial mean',
    )

    # Canonical anchor: if the seed mean is in the wrong hemisphere, flip it.
    # This is the fix for majority-vote bias.  The seed now always points in
    # the canonical direction before iteration begins, so Phase 1 converges
    # toward the canonical hemisphere regardless of the raw N/R ratio.
    mean_lat, mean_lon = force_canonical(mean_lat, mean_lon, canonical_hemisphere)

    prev_flip_set = None
    iteration = 0
    for iteration in range(1, MAX_ITER + 1):
        flip_set = frozenset(
            i for i in valid_idx
            if angular_distance(
                raw_coords[i][0], raw_coords[i][1], mean_lat, mean_lon
            ) > 90.0
        )
        if flip_set == prev_flip_set:
            break
        prev_flip_set = flip_set

        c_lats, c_lons = [], []
        for i in valid_idx:
            la, lo = raw_coords[i]
            if i in flip_set:
                la, lo = flip(la, lo)
            c_lats.append(la)
            c_lons.append(lo)

        mean_lat, mean_lon = safe_mean(
            c_lats, c_lons, canonical_hemisphere,
            context=f'Phase-1 iter {iteration}',
        )

    # Final canonical check: if convergence drifted to wrong hemisphere, correct.
    if not is_canonical(mean_lat, canonical_hemisphere):
        mean_lat, mean_lon = flip(mean_lat, mean_lon)
        prev_flip_set = frozenset(
            (set(valid_idx) - prev_flip_set) if prev_flip_set else set(valid_idx)
        )

    return set(prev_flip_set), (mean_lat, mean_lon), iteration


# ---------------------------------------------------------------------------
# PHASE 2 — LOCAL COEVAL WINDOW REFINEMENT & CYCLE DETECTION
# ---------------------------------------------------------------------------

def phase2_local(
    raw_coords, valid_idx, initial_flip_set, global_mean,
    families, ages, canonical_hemisphere,
):
    """
    Iterative local coeval window refinement with full cycle detection.

    Each iteration:
      1. Build corrected coordinates from the current flip_set.
      2. Recompute the global mean from corrected coordinates.
      3. For each pole, build a local reference from coeval corrected peers
         (same craton family, ±WINDOW_MA Ma, >=MIN_PEERS); otherwise fall
         back to the global mean.
      4. Classify RAW coordinates against the local/global reference.
      5. Form a new flip_set.
      6. If new == current  -> converged.
         If new matches any earlier state -> cycle detected.
         Poles that changed on the last non-converging step are marked UNSTABLE.

    Parameters
    ----------
    raw_coords        : list[tuple[float, float] or None]
    valid_idx         : list[int]
    initial_flip_set  : set[int]
    global_mean       : tuple[float, float]  — from Phase 1 (unused after iter 1)
    families          : list[str or None]
    ages              : list[float or None]
    canonical_hemisphere : str — 'N' or 'S'

    Returns
    -------
    flip_set      : set[int]   — final row indices of poles to flip
    unstable_set  : set[int]   — row indices of oscillating poles
    iterations    : int
    converged     : bool
    """
    fam_idx = defaultdict(list)
    for i in valid_idx:
        if families[i]:
            fam_idx[families[i]].append(i)

    flip_set = set(initial_flip_set)
    history = [frozenset(flip_set)]   # full state history for cycle detection
    converged = False
    last_changed = set()

    for iteration in range(1, MAX_ITER + 1):

        # Build corrected coordinates from current flip_set
        corrected = {}
        for i in valid_idx:
            la, lo = raw_coords[i]
            if i in flip_set:
                la, lo = flip(la, lo)
            corrected[i] = (la, lo)

        # Global mean from corrected coordinates
        c_lats = [corrected[i][0] for i in valid_idx]
        c_lons = [corrected[i][1] for i in valid_idx]
        cur_global_lat, cur_global_lon = safe_mean(
            c_lats, c_lons, canonical_hemisphere,
            context=f'Phase-2 iter {iteration} global mean',
        )

        # Classify each pole
        new_flip_set = set()
        for i in valid_idx:
            raw_lat, raw_lon = raw_coords[i]
            fam = families[i]
            age = ages[i]

            ref_lat, ref_lon = cur_global_lat, cur_global_lon  # default reference

            if fam and age is not None:
                peers = [
                    j for j in fam_idx[fam]
                    if j != i
                    and ages[j] is not None
                    and abs(ages[j] - age) <= WINDOW_MA
                ]
                if len(peers) >= MIN_PEERS:
                    p_lats = [corrected[j][0] for j in peers]
                    p_lons = [corrected[j][1] for j in peers]
                    ref_lat, ref_lon = safe_mean(
                        p_lats, p_lons, canonical_hemisphere,
                        context=f'Phase-2 iter {iteration} local mean row {i}',
                    )

            # Always classify RAW coordinates — prevents double-flip
            if angular_distance(raw_lat, raw_lon, ref_lat, ref_lon) > 90.0:
                new_flip_set.add(i)

        new_fs = frozenset(new_flip_set)
        last_changed = set(new_flip_set).symmetric_difference(flip_set)

        # Convergence check
        if new_fs == frozenset(flip_set):
            converged = True
            break

        # Cycle detection: has this exact state appeared before?
        if new_fs in history:
            break   # converged=False; last_changed holds the oscillating poles

        history.append(new_fs)
        flip_set = new_flip_set

    return flip_set, last_changed, iteration, converged


# ---------------------------------------------------------------------------
# FILE PROCESSING
# ---------------------------------------------------------------------------

def process_file(filepath):
    """
    Run the full Polarity Standardizer pipeline on a single CSV file.

    Reads the file, runs Phase 1 and Phase 2, writes the output CSV with
    a PTest column added, and prints a progress summary to stdout.
    """
    filename = os.path.basename(filepath)
    print(f"\n{'-' * 70}")
    print(f"  Input : {filepath}")

    try:
        fieldnames, rows = read_csv_file(filepath)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    n = len(rows)
    raw_coords = []
    ages = []
    families = []
    n_skipped = 0

    for row in rows:
        fam = craton_family(row.get('Craton', ''))
        age = parse_age(row)
        try:
            raw_coords.append(parse_coords(row))
        except (ValueError, TypeError):
            raw_coords.append(None)
            n_skipped += 1
        ages.append(age)
        families.append(fam)

    valid_idx = [i for i, c in enumerate(raw_coords) if c is not None]

    if not valid_idx:
        print("  WARNING: no parseable poles found — nothing written.")
        return

    # Phase 1
    p1_flips, global_mean, p1_iters = phase1_global(
        raw_coords, valid_idx, CANONICAL_HEMISPHERE
    )
    print(
        f"  Phase-1 (global)  : {p1_iters} iter  |  "
        f"flipped: {len(p1_flips)}  |  "
        f"mean Plat={global_mean[0]:.2f}  Plon={global_mean[1]:.2f}"
    )

    # Phase 2
    final_flips, unstable_set, p2_iters, converged = phase2_local(
        raw_coords, valid_idx, p1_flips, global_mean,
        families, ages, CANONICAL_HEMISPHERE,
    )

    status = "converged" if converged else "CYCLE DETECTED - non-converged"
    print(
        f"  Phase-2 (local)   : {p2_iters} iter  |  "
        f"final reversed: {len(final_flips)}  |  "
        f"unstable: {len(unstable_set)}  |  {status}"
    )

    if not converged:
        warnings.warn(
            f"Phase-2 did not converge for '{filename}' after {p2_iters} iterations. "
            f"{len(unstable_set)} pole(s) are oscillating. "
            "They are flagged PTest='U' in the output.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Build output rows
    n_reversed = n_unstable = 0
    out_rows = []
    for i, row in enumerate(rows):
        out_row = dict(row)
        if raw_coords[i] is None:
            out_row['PTest'] = ''
        elif i in unstable_set:
            # Oscillating pole — output last assigned coords, flag as unstable
            if i in final_flips:
                lat, lon = flip(raw_coords[i][0], raw_coords[i][1])
                out_row['Plat'] = preserve_fmt(row['Plat'], lat)
                out_row['Plon'] = preserve_fmt(row['Plon'], lon)
            out_row['PTest'] = 'U'
            n_unstable += 1
        elif i in final_flips:
            lat, lon = flip(raw_coords[i][0], raw_coords[i][1])
            out_row['Plat'] = preserve_fmt(row['Plat'], lat)
            out_row['Plon'] = preserve_fmt(row['Plon'], lon)
            out_row['PTest'] = 'R'
            n_reversed += 1
        else:
            out_row['PTest'] = 'N'
        out_rows.append(out_row)

    out_fieldnames = fieldnames + (['PTest'] if 'PTest' not in fieldnames else [])

    # Write output to the same directory as the input file
    input_dir = os.path.dirname(os.path.abspath(filepath))
    out_dir = input_dir if os.path.isdir(input_dir) else os.getcwd()
    base = os.path.splitext(filename)[0]
    out_name = f"PolarityStd_{base}.csv"
    out_path = os.path.join(out_dir, out_name)
    write_csv_file(out_path, out_fieldnames, out_rows)

    n_normal = n - n_reversed - n_unstable - n_skipped
    print(f"  Output: {out_path}")
    print(
        f"  Total: {n}  |  Normal (N): {n_normal}  |  "
        f"Reversed (R): {n_reversed}  |  "
        f"Unstable (U): {n_unstable}  |  Skipped: {n_skipped}"
    )


# ---------------------------------------------------------------------------
# FILE DISCOVERY
# ---------------------------------------------------------------------------

def find_files(args):
    """
    Resolve the list of CSV files to process from command-line arguments.

    If the first argument is a directory, scans it for ConfidentInterval*.csv.
    If explicit file paths are given, uses those directly.
    Falls back to scanning the current working directory.

    Parameters
    ----------
    args : list[str] — typically sys.argv[1:]

    Returns
    -------
    list[str] — sorted, deduplicated list of absolute file paths
    """
    if args and os.path.isdir(args[0]):
        search_dir = args[0]
        explicit = []
    else:
        search_dir = os.getcwd()
        explicit = [a for a in args if os.path.isfile(a)]

    found = list(explicit)
    if not explicit:
        found += glob.glob(os.path.join(search_dir, 'ConfidentIntervals*.csv'))
        found += glob.glob(os.path.join(search_dir, 'ConfidentInterval_*.csv'))

    seen, unique = set(), []
    for f in found:
        key = os.path.abspath(f)
        if key not in seen and not os.path.basename(f).startswith('~$'):
            seen.add(key)
            unique.append(f)
    return sorted(unique)


# ---------------------------------------------------------------------------
# ENTRY POINT
# ---------------------------------------------------------------------------

def main():
    """Run the Polarity Standardizer on all discovered or specified CSV files."""
    warnings.filterwarnings('always', category=RuntimeWarning)

    if any(a in ('-h', '--help') for a in sys.argv[1:]):
        print(__doc__)
        sys.exit(0)

    script_name = os.path.basename(sys.argv[0])

    print("=" * 70)
    print("  POLARITY STANDARDIZER  v4")
    print("  (Coordinate hemisphere-unification for palaeomagnetic VGPs)")
    print(f"  Canonical hemisphere : {CANONICAL_HEMISPHERE}")
    print("  Phase-1 : global seed with canonical-hemisphere anchor")
    print("  Phase-2 : local coeval refinement with cycle detection")
    print("=" * 70)

    files = find_files(sys.argv[1:])
    if not files:
        cwd = os.getcwd()
        print(f"\nNo ConfidentIntervals*.csv files found in: {cwd}")
        print(f"Usage:  python {script_name}")
        print(f"        python {script_name} \"path/to/folder\"")
        print(f"        python {script_name} file1.csv file2.csv")
        sys.exit(1)

    print(f"\nFound {len(files)} file(s).")
    for fp in files:
        process_file(fp)

    print(f"\n{'=' * 70}  Done.\n")


if __name__ == '__main__':
    main()
