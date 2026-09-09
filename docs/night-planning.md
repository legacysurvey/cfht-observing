# CFHT night planning

`make_cfht_targets.py --night YYYY-MM-DD` produces a queue for the night
beginning on that calendar date in Hawaii (HST, UTC−10). It prints the
evening and following morning twilight times in HST and UTC, the local
**mean** sidereal time at both endpoints in hours and degrees, and the
number of complete exposure-plus-overhead slots that fit.

The existing Anaconda environment with Astropy 4.3.1 and NumPy 1.20.3 is
supported. In a new Python environment, install the astronomy dependencies:

```sh
python -m pip install -r requirements.txt
```

Populate priorities and design times, then generate a night:

```sh
python update_tile_priorities.py
python update_lmst_design.py
python make_cfht_targets.py --night 2026-09-09
```

The night mode writes three files in `plans/`, using the night date in
their names unless `--prefix` is supplied:

- `targets-YYYY-MM-DD.xml`: the existing CFHT target-import format.
- `targets-YYYY-MM-DD.fits`: the same queue and existing target-list columns.
- `targets-YYYY-MM-DD.schedule.ecsv`: every planned slot, including empty
  slots, exposure start/end and slot end in UTC, object name, midpoint LMST,
  design LMST, their signed difference, minimum endpoint altitude, priority,
  and status.

FITS and XML rows follow scheduled time order. The schedule is a planning
record; the target import format itself does not enforce start times or
set exposure durations in the CFHT queue. The queue exposure settings must
use the same exposure time assumed here. Existing outputs require `--overwrite`.

## Selection and timing

The default program is now **LBNL** in both ordinary and night modes.
`--program NAOC`, `--program LBNL NAOC`, or other names select an explicit
list; `--all-programs` includes all programs and blank labels, and
`--program ""` selects only blank labels. Matching ignores case.
The existing filter (default M4376), RA, DEC, `IN_IBIS=1`, HSC exclusion,
and `DONE=0` cuts apply before scheduling. RA bounds may wrap through zero.

Defaults are `--twilight 15`, `--exptime 130`, and `--overhead 44`.
Capacity is the integer floor of elapsed night seconds divided by the
174-second cadence. For a nine-hour night, this gives 186 slots and
36 seconds left over. Overhead follows each exposure, including the last;
no exposure or overhead extends past morning twilight.

Each slot is anchored to evening twilight. Candidate design times are
compared circularly to actual LMST at the **exposure midpoint**. The
default `--lmst-window 5` permits an absolute difference of up to five
degrees (twenty sidereal minutes). Among available targets in the window,
the highest numeric `PRIORITY` wins. Ties are broken by lowest declination,
then closest LMST, then object name. A target must be above the geometric horizon at the
start and end of its exposure. A target is used at most once per night.

`--non-overlapping` applies the existing projected rectangular separation
rule after each selection, so survey priority is honored in night mode.
It can reduce the number of slots filled. The window can also be changed:

```sh
python make_cfht_targets.py --night 2026-09-09 --twilight 18 \
    --exptime 150 --overhead 44 --lmst-window 3 \
    --ra 300 60 --dec -20 0 --program LBNL --non-overlapping
```

If no nearby, visible, unused target satisfies the cuts, that slot is
explicitly left empty. The summary reports the shortfall and total unused
time; the scheduler never silently selects outside the supplied cuts or
LMST window to claim a full night. There are currently no Moon, weather,
airmass, detailed telescope-limit, variable-slew, or filter-change models.
Ordinary (non-night) plans also use priority before greedy non-overlap
trimming, then restore RA order for output. Without trimming, ordinary
plans include every matching target in RA order.

## Survey-year priorities

`update_tile_priorities.py` populates the existing float64 `PRIORITY`
column in both tile files from the `YEAR` column of `desi2-stripes.ecsv`.
Y1 through Y5 correspond to 10.0, 9.0, 8.0, 7.0, and 6.0. All `IN_IBIS=1`
rows receive a priority, including completed targets and HSC areas;
planning still independently applies its normal completion and HSC cuts.

By default, each year's stripe uses its tabulated RA bounds and central
declination +/-2.1 degrees (1.6 for the DESI tile radius, plus 0.5 for the
CFHT half-footprint). The SGC RA bounds wrap through zero. Where stripes
overlap, the highest priority always wins, regardless of the order in
which the stripes are processed. In particular, a Y1 match cannot be
overwritten by a later-year match. IBIS rows outside these bounds get
priority 1.0; rows outside IBIS keep their existing priority.

The optional `--y1-only-padding` updater flag uses +/-2.1 degrees for Y1
and the table's DEC_MIN/DEC_MAX for later years. Recomputing priorities
does not change any other tile columns, including `LMST_DESIGN`.

## Rebalancing LMST_DESIGN

`update_lmst_design.py` updates both `obstatus/cfht-tiles.ecsv` and
`obstatus/cfht-tiles.fits`. It uses remaining rows with `IN_IBIS=1`,
`IN_HSC=0`, and `DONE=0`, independently in each filter, across all programs.
The ECSV table is authoritative for this calculation because it includes
five completion updates not yet present in FITS. Both formats get identical
designs; all pre-existing column values are preserved independently.

For each filter:

1. Find the largest empty RA gap and unwrap RA at its midpoint.
2. Sort by unwrapped RA, preserving input order for equal RAs.
3. Assign equally spaced circular ranks `(rank + 0.5) × 360/N` degrees.
4. Rotate the entire grid so the arithmetic mean of signed hour angles
   `wrap(LMST_DESIGN − catalog RA, −180°, 180°)` is zero.

The resulting LMST spacing is exactly 360/N degrees. This is a survey-wide
load-balancing design, not a request to observe each tile at transit.
The current M4376 mapping sends RA≈0° to LMST≈1.6° and RA≈59° to LMST≈86.5°.
The mapping is independent of any particular night or program cut, so a
restricted selection need not retain a uniform LMST distribution.

Design values are float64 degrees in [0, 360). Ineligible or completed rows
have NaN. Re-run the updater after completion or footprint changes; night
mode requires a finite design for every selected candidate. Designs are
planning angles relative to the catalog RA, rather than an epoch-dependent
astrometric hour-angle solution.

## Ephemeris and validation

The site is longitude −155°28′18.00″, latitude +19°49′41.86″, elevation
4204 m, from the [CFHT Observatory Manual](https://www.cfht.hawaii.edu/Instruments/ObservatoryManual/CFHT_ObservatoryManual_%28Sec_1%29.html).
Twilight is defined by the topocentric altitude of the Sun's center, using
[Astropy's built-in solar ephemeris](https://docs.astropy.org/en/stable/api/astropy.coordinates.get_sun.html)
and [AltAz with refraction disabled](https://docs.astropy.org/en/stable/api/astropy.coordinates.AltAz.html).
There is no horizon-dip or solar-radius correction to the specified solar
depression. Crossings are bisected to a 0.02-second time bracket.

[Astropy mean sidereal time](https://docs.astropy.org/en/stable/time/index.html)
uses the project's [offline IERS and leap-second tables](../data/iers/README.md),
including on older Astropy versions. Normal planning does not download
anything. Dates outside the local data coverage are rejected before any
coordinate transformations. New data can be copied into `data/iers/`
manually, or downloaded explicitly with the optional `update_iers_data.py`
utility. That utility is never invoked by the planner. The bisection
tolerance is numerical precision, not a claim of millisecond ephemeris
or operational accuracy.

Run the checks with:

```sh
python -m unittest -v test_make_cfht_targets test_update_tile_columns \
    test_lmst_design test_tile_priorities test_night_planning
```
