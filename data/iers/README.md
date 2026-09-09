# Offline Earth-orientation data

Night planning uses these local files with Astropy downloads disabled.
They work with the user's Astropy 4.3.1 / NumPy 1.20.3 Anaconda environment
as well as newer Astropy releases. No package upgrade or network access is
required to generate the September 9, 2026 plan.

The initial copies came from the already-installed
`astropy-iers-data` version `0.2026.8.31.0.57.9`, on September 9, 2026:

- `finals2000A.all`: IERS Earth orientation, including UT1−UTC and polar
  motion; valid rows cover MJD 41684 through 61645. Provider:
  [USNO/IERS Rapid Service](https://maia.usno.navy.mil/ser7/finals2000A.all).
- `Leap_Second.dat`: TAI−UTC values, updated through IERS Bulletin C 72
  (July 2026), expiring June 28, 2027. Provider:
  [IERS Earth Orientation Centre](https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat).

The planner checks the requested date against both files before computing
coordinates. It does not silently substitute zero UT1−UTC or extrapolate
beyond these data. Updating Astropy is optional; the planner explicitly
uses this project snapshot for reproducible offline results.

When newer dates need new data, the files can be replaced manually.
`python update_iers_data.py` is an optional, explicitly invoked network
utility that refreshes both files from the providers above. It is never
called by the planner and validates both downloads before replacing either
existing file. Normal planning continues to work offline afterward.
