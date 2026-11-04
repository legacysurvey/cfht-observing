# cfht-observing
Observing logs for the DESI imaging on CFHT

[`desi2-stripes.ecsv`](desi2-stripes.ecsv) contains all 18 stripes of
Table 1 (`tab:stripes`) of the *DESI Run 2 High-z Imaging Plan*
(`desi_run2_imaging_plan/sample631.tex`, as of 2026-11-04), including the
proposed NGC-9 and SGC-9 extensions. It preserves the published coordinates,
areas, survey years, assignments, and notes, with angular units and source
details in the ECSV header, and adds `POLU_IBIS`: 1 where CFHT (POLU) and
DECam (IBIS) image the stripe outside the HSC-Niji+VST rectangles, 0 for
SGC-1 (HSC-Niji+VST only; it lies largely outside the Rubin footprint) and
for the extensions, which are not in any tile file yet.

[`footprint-regions.ecsv`](footprint-regions.ecsv) lists the areas inside
those stripes that CFHT and DECam do not image: the six HSC-Niji+VST
rectangles of Table 2 (`tab:niji_stripes`; NGC-5/6 at RA 128-226, SGC-1/2
at RA 330.2-39.7, and the HSC-SSP area south of SGC-2 at RA 28.7-39.7 as
parts of SGC-3/4) and the two dropped slivers (NGC-5 west of RA 128, NGC-6
east of RA 226). `footprint.py` turns the two files into a raster of
POLU+IBIS sky, shared by every camera.

`cfht-tiles.ecsv` is the sole tile catalog. It includes:

- `PRIORITY`: float64 survey-year priority from `desi2-stripes.ecsv`:
  Y1=10.0, Y2=9.0, Y3=8.0, Y4=7.0, and Y5=6.0. Every row with
  `IN_IBIS=1` participates, regardless of `DONE`, `IN_HSC`, filter, or
  program. Stripe membership uses the tabulated RA bounds (including SGC
  wraparound) and `DEC_CENTER +/- (1.6 + 0.5)` degrees for every year.
  Boundaries are inclusive, and the maximum priority wins at overlaps.
  IBIS rows outside all stripes get 1.0; rows outside IBIS retain their
  previous priority. Recompute with `python update_tile_priorities.py`.
  The optional `--y1-only-padding` flag uses tabulated DEC_MIN/DEC_MAX for
  later years instead, while keeping the +/-2.1-degree Y1 bounds.
- `PROGRAM`: `NAOC` for exact `OBJECT` matches to target `NAME` values in
  the FITS or XML files in `plans/`; otherwise `LBNL` where `IN_IBIS=1`.
  Other rows have an empty string. Plan membership takes precedence and is
  specific to the named target, including its filter.
- `IN_HSC`: an int16 flag, 1 for `IN_IBIS=1` tiles that CFHT should not
  observe. A tile is a MegaCam program tile (`IN_IBIS=1, IN_HSC=0`) when its
  center lies within one degree of POLU+IBIS sky (the same one-degree margin
  as `IN_IBIS`, so that fields reach the tabulated stripe bounds) unless its
  1 x 1 degree field lies entirely inside an HSC-Niji+VST rectangle. The whole
  of SGC-1 is therefore `IN_HSC=1` except where the margin of SGC-2 reaches
  into it. Despite its name the flag also covers sky nobody images (the
  SGC-1 ends and the dropped slivers). Recompute with
  `python update_tile_columns.py`.
- `LMST_DESIGN`: float64 local mean sidereal time in degrees, in [0, 360),
  balanced uniformly for remaining eligible rows separately in each filter.
  Rows outside `IN_IBIS=1, IN_HSC=0, DONE=0` have NaN.
  All programs participate in this mapping. Recompute with
  `python update_lmst_design.py` as completion or footprint flags change.
  See [night planning](../docs/night-planning.md) for the algorithm.

`PROGRAM` labels are sticky: `update_tile_columns.py` only labels rows that
have no label yet (`NAOC` for exact `OBJECT` matches to target `NAME` values
in the FITS or XML files in `plans/`, otherwise `LBNL` where `IN_IBIS=1`),
so tiles that enter or leave the program keep their existing label.

To regenerate the flags after a change to the plan tables (requires numpy
and astropy):

```sh
python update_tile_columns.py
python update_tile_priorities.py
python update_lmst_design.py
python -m unittest -v tests.test_update_tile_columns
```

The tile updaters read and write only `cfht-tiles.ecsv`. Each validates a
staged catalog before replacing the original and preserves all columns
other than those it updates. Observing plans still use FITS and XML.
