# cfht-observing
Observing logs for the DESI imaging on CFHT

[`desi2-stripes.ecsv`](desi2-stripes.ecsv) contains
all 18 stripes from the table captioned Table 2 on page 3 of the September
9, 2026 imaging-plan draft, including the NGC-9 and SGC-9 extensions. It
preserves the published coordinates, areas, survey years, assignments, and
notes, with angular units and source details in the ECSV header. These are
the document's original bounds, before any CFHT overlap margins or
temporary HSC extensions. This file is also used by the footprint updater.

`cfht-tiles.ecsv` and `cfht-tiles.fits` include:

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
- `IN_HSC`: an int16 flag, set to 1 inside the HSC+VST footprint after a
  one-degree inset in both RA and declination, plus the temporary northern
  extensions described below, and 0 elsewhere. It is computed for all tile
  rows, independently of `PROGRAM`; `IN_IBIS` sets the SGC extension limit.
- `LMST_DESIGN`: float64 local mean sidereal time in degrees, in [0, 360),
  balanced uniformly for remaining eligible rows separately in each filter.
  Rows outside `IN_IBIS=1, IN_HSC=0, DONE=0` have NaN. The ECSV table is
  the source for eligibility, including its more recent `DONE` values;
  identical design values are written to FITS without changing its other
  columns. All programs participate in this mapping. Recompute with
  `python update_lmst_design.py` as completion or footprint flags change.
  See [night planning](../docs/night-planning.md) for the algorithm.

The HSC+VST assignments come from the stripe table on page 3 (captioned
Table 2) of *DESI Run 2 High-Z Imaging Plan*, draft September 9, 2026,
supplied as `DESI_Run_2_Imaging_Plan.pdf`. The four Niji+VST stripes are:

| Stripe | DEC center (deg) | RA min (deg) | RA max (deg) |
| --- | ---: | ---: | ---: |
| NGC-5 | 2.08 | 126.60 | 249.93 |
| NGC-6 | -0.72 | 128.16 | 228.30 |
| SGC-1 | 3.29 | 315.57 | 44.37 |
| SGC-2 | 0.49 | 313.82 | 47.15 |

Use a stripe half-height of 1.4 degrees about `DEC_CENTER`, rather than the wider
`DEC_MIN`/`DEC_MAX` bounds in `desi2-stripes.ecsv`. Merge adjoining stripes
before applying declination insets; shared edges have no inset, while
exposed steps from unequal RA limits retain the margin. RA ranges wrap
through 360 degrees in the SGC. Margins are in coordinate degrees, and
boundary values are included.

In the original footprint, the NGC interior extends from DEC=-1.12 to 2.48 degrees where
RA is between 129.16 and 227.30 degrees. The wider NGC-5-only portions,
within RA=127.60 to 248.93 degrees, extend from DEC=1.68 to 2.48 degrees.
The 2.48-degree northern edge follows the tabulated center of 2.08 degrees.

The current footprint temporarily gives HSC the opportunity to observe
farther north. NGC-5 extends to DEC=4.0 degrees, retaining its inset RA
range of 127.60 to 248.93 degrees. SGC-1 extends to the maximum DEC of any
SGC tile with `IN_IBIS=1` (currently 6.121 degrees), retaining its inset RA
range of 316.57 through 360 to 43.37 degrees. SGC tiles are identified by
RA >= 270 or RA < 90 degrees. These are final northern limits, with no
additional one-degree inset. Southern edges and RA boundaries stay as in
the original footprint. The extension can include rows with `IN_IBIS=0`.

To regenerate the current extended footprint and program labels from the
current plans (requires numpy and astropy):

```sh
python update_tile_columns.py --hsc-stripes NGC-5 NGC-6 SGC-1 SGC-2 --extend-hsc-north
python -m unittest -v tests.test_update_tile_columns
```

To restore the original Table 2 footprint, omit the extension switch:

```sh
python update_tile_columns.py --hsc-stripes NGC-5 NGC-6 SGC-1 SGC-2
```

The updater validates staged files before replacing the originals and
preserves each file's existing column values. In particular, five `DONE`
values already differed between ECSV and FITS before this update; those
differences are preserved.
