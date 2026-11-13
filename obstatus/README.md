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
- `IN_IBIS`: an int16 flag, 1 where the tile's 1 x 1 degree field overlaps
  the imaged DESI Run 2 footprint: POLU+IBIS sky or an HSC-Niji+VST
  rectangle (`footprint.py`). Tiles beyond the padded stripe bounds, in the
  dropped slivers, or at the unimaged ends of SGC-1 are 0. Recompute with
  `python update_tile_columns.py`.
- `IN_HSC`: an int16 flag, 1 for `IN_IBIS=1` tiles whose footprint sky is
  HSC-Niji+VST's alone, so CFHT does not observe them: the field lies
  entirely inside a Niji rectangle or overlaps the footprint only through
  one (the whole of SGC-1, for example). A tile is a MegaCam program tile
  when `IN_IBIS=1` and `IN_HSC=0`, i.e. its field overlaps POLU+IBIS sky
  (center within half a field of it, the same rule DECam uses with its own
  field size) and does not lie entirely inside a Niji rectangle.
- `DONE`, `DATE`, `EXPNUM`, `EXPTIME`, `SNR`: the observation that counts for
  the tile, set by `end-of-night.py`. An exposure counts when it is good
  (graded by CFHT, not grade 5, Elixir SNR > 4.2 = SNR goal 4.6 minus 0.4) or
  when CFHT validated it at lower SNR; the highest-SNR such exposure is
  recorded. `SNR <= 4.2` on a `DONE` tile marks one to be topped up later.
- `LMST_DESIGN`: float64 local mean sidereal time in degrees, in [0, 360),
  balanced uniformly for remaining eligible rows separately in each filter.
  Rows outside `IN_IBIS=1, IN_HSC=0, DONE=0` have NaN.
  All programs participate in this mapping. Recompute with
  `python update_lmst_design.py` as completion or footprint flags change.
  See [night planning](../docs/night-planning.md) for the algorithm.

`PROGRAM` labels are sticky for program tiles: `update_tile_columns.py`
only labels program rows (`IN_IBIS=1, IN_HSC=0`) that have no label yet
(`NAOC` for exact `OBJECT` matches to target `NAME` values in the FITS or XML
files in `plans/`, otherwise `LBNL`). Tiles outside the program carry no
label unless they have been observed: a tile that leaves the footprint
unobserved is cleared to the empty string, since it is not to be observed,
while an observed one keeps the label it was observed under.

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
other than those it updates (`update_tile_columns.py`: `IN_IBIS`, `IN_HSC`,
`PROGRAM`). Observing plans still use FITS and XML.
