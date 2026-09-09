"""Mauna Kea twilight and CFHT queue scheduling; requires numpy and astropy."""

from contextlib import contextmanager
from dataclasses import dataclass
import datetime as dt
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.table import Table
from astropy.time import Time
from astropy.utils import iers


# CFHT Observatory Manual, section 1 (west longitude is negative).
CFHT = EarthLocation.from_geodetic(-(155 + 28 / 60 + 18 / 3600) * u.deg,
                                 (19 + 49 / 60 + 41.86 / 3600) * u.deg,
                                 4204 * u.m)
HST = dt.timezone(dt.timedelta(hours=-10), name='HST')
IERS_DIRECTORY = Path(__file__).resolve().parent / 'data' / 'iers'


@lru_cache(maxsize=1)
def _iers_table():
    path = IERS_DIRECTORY / 'finals2000A.all'
    if not path.is_file():
        raise ValueError('Missing project IERS data; run python update_iers_data.py')
    return iers.IERS_A.open(str(path))


@contextmanager
def offline_iers():
    """Use project tables consistently on old and new Astropy versions."""
    leap_file = IERS_DIRECTORY / 'Leap_Second.dat'
    if not leap_file.is_file():
        raise ValueError('Missing project leap-second data; run python update_iers_data.py')
    # A finite auto_max_age avoids overflow in Astropy 4.x leap-second
    # initialization. We explicitly load IERS_A, so auto-refresh age is
    # irrelevant here; coverage is checked before any transformations.
    with iers.conf.set_temp('auto_download', False), \
            iers.conf.set_temp('auto_max_age', 30), \
            iers.conf.set_temp('system_leap_second_file', str(leap_file)):
        with iers.earth_orientation_table.set(_iers_table()):
            yield


def check_iers_coverage(times):
    """Reject unsupported dates before Astropy can substitute UT1-UTC=0."""
    _, status = times.get_delta_ut1_utc(return_status=True)
    if np.any(np.asarray(status) < 0):
        raise ValueError('Night is outside project IERS coverage; '
                         'run python update_iers_data.py')
    leap_seconds = iers.LeapSeconds.open(str(IERS_DIRECTORY / 'Leap_Second.dat'))
    # LeapSeconds expresses expiration in TAI. Compare in that scale to
    # avoid converting the future expiration date through old ERFA code.
    if np.any(times.tai > leap_seconds.expires):
        raise ValueError('Night is past the project leap-second table expiration; '
                         'run python update_iers_data.py')


@dataclass(frozen=True)
class Night:
    date: str
    twilight: float
    start: Time
    end: Time
    lmst_start: float
    lmst_end: float

    @property
    def seconds(self):
        return float((self.end - self.start).to_value(u.s))


def solar_altitude(time):
    return get_sun(time).transform_to(
        AltAz(obstime=time, location=CFHT, pressure=0 * u.hPa)).alt.deg


def night_bounds(date, twilight=15.0):
    """Sun-center crossings at -twilight, from local noon to next noon."""
    day = dt.date.fromisoformat(date)
    if not math.isfinite(twilight) or not 0 < twilight < 90:
        raise ValueError('Twilight must be between 0 and 90 degrees')
    noon = Time(dt.datetime.combine(day, dt.time(12), tzinfo=HST))
    with offline_iers():
        # Construct UTC bounds without time-scale arithmetic, so the
        # coverage check precedes any Earth-orientation transformations.
        following_noon = Time(dt.datetime.combine(day + dt.timedelta(days=1),
                                                  dt.time(12), tzinfo=HST))
        check_iers_coverage(Time([noon, following_noon]))
        offsets = np.arange(0, 86401, 600)
        times = noon + offsets * u.s
        values = solar_altitude(times) + twilight
        evening = np.flatnonzero((values[:-1] > 0) & (values[1:] <= 0))
        morning = np.flatnonzero((values[:-1] <= 0) & (values[1:] > 0))
        if not len(evening) or not len(morning):
            raise ValueError(f'No complete {twilight:g}-degree twilight night on {date}')
        left = int(evening[0])
        after_evening = morning[morning > left]
        if not len(after_evening):
            raise ValueError('Morning twilight did not follow evening twilight')

        def crossing(index):
            low, high = float(offsets[index]), float(offsets[index + 1])
            positive = values[index] > 0
            while high - low > 0.02:
                middle = (low + high) / 2
                if (solar_altitude(noon + middle * u.s) + twilight > 0) == positive:
                    low = middle
                else:
                    high = middle
            return noon + ((low + high) / 2) * u.s

        start, end = crossing(left), crossing(int(after_evening[0]))
        return Night(date, twilight, start, end,
                     float(start.sidereal_time('mean', longitude=CFHT.lon).deg),
                     float(end.sidereal_time('mean', longitude=CFHT.lon).deg))


def slot_count(seconds, exptime=130.0, overhead=44.0):
    if not all(math.isfinite(v) for v in (seconds, exptime, overhead)):
        raise ValueError('Durations must be finite')
    if seconds < 0 or exptime <= 0 or overhead < 0:
        raise ValueError('Night/overhead must be nonnegative and exposure positive')
    # Astropy time subtraction can put an exact slot boundary a few
    # picoseconds below its mathematical value.
    return math.floor((seconds + 1e-8) / (exptime + overhead))


def pick_target(targets, lmst, available, window):
    """Highest priority nearby design; tie-break by DEC, distance, name."""
    nearby = [i for i in available if abs(
        (targets[i].lmst_design_deg - lmst + 180) % 360 - 180) <= window]
    return sorted(nearby, key=lambda i: (-targets[i].priority, targets[i].dec_deg, abs(
        (targets[i].lmst_design_deg - lmst + 180) % 360 - 180), targets[i].object_name))


def schedule_night(targets, night, exptime=130.0, overhead=44.0,
                   lmst_window=5.0, non_overlapping=False, min_separation=1.0):
    """Return queue targets in time order plus an audit table of every slot.

    The exposure begins at slot start; overhead follows it. Selection uses
    LMST at exposure midpoint. A target must be above the geometric horizon
    at both exposure endpoints. Empty slots remain explicit in the audit.
    """
    if not math.isfinite(lmst_window) or not 0 < lmst_window <= 180:
        raise ValueError('LMST window must be in (0, 180] degrees')
    if len({t.object_name for t in targets}) != len(targets):
        raise ValueError('Candidate OBJECT names must be unique')
    if any(t.lmst_design_deg is None or not math.isfinite(t.lmst_design_deg) for t in targets):
        raise ValueError('Candidates require finite LMST_DESIGN; run update_lmst_design.py')
    if any(not math.isfinite(t.priority) for t in targets):
        raise ValueError('Candidates require finite PRIORITY; run update_tile_priorities.py')
    count = slot_count(night.seconds, exptime, overhead)
    cadence = exptime + overhead
    ra = np.asarray([t.ra_deg for t in targets])
    dec = np.asarray([t.dec_deg for t in targets])
    coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    available = set(range(len(targets)))
    selected, rows = [], []
    with offline_iers():
        check_iers_coverage(Time([night.start, night.end]))
        starts = night.start + np.arange(count) * cadence * u.s
        ends = starts + exptime * u.s
        midpoints = starts + (exptime / 2) * u.s
        lmsts = midpoints.sidereal_time('mean', longitude=CFHT.lon).deg
        for slot in range(count):
            nearby = pick_target(targets, float(lmsts[slot]), available, lmst_window)
            chosen, altitude = None, np.nan
            if nearby:
                alt_start = coords[nearby].transform_to(AltAz(
                    obstime=starts[slot], location=CFHT, pressure=0 * u.hPa)).alt.deg
                alt_end = coords[nearby].transform_to(AltAz(
                    obstime=ends[slot], location=CFHT, pressure=0 * u.hPa)).alt.deg
                visible = np.flatnonzero((alt_start > 0) & (alt_end > 0))
                if len(visible):
                    position = int(visible[0])
                    chosen = nearby[position]
                    altitude = float(min(alt_start[position], alt_end[position]))
            if chosen is None:
                name, design, offset, priority = '', np.nan, np.nan, np.nan
                status = 'BELOW_HORIZON' if nearby else 'NO_NEARBY_TARGET'
            else:
                target = targets[chosen]
                selected.append(target)
                available.remove(chosen)
                name, design = target.object_name, target.lmst_design_deg
                priority = target.priority
                offset = (design - lmsts[slot] + 180) % 360 - 180
                status = 'SCHEDULED'
                if non_overlapping:
                    ra_sep = np.abs((ra - target.ra_deg + 180) % 360 - 180)
                    ra_sep *= np.abs(np.cos(np.deg2rad((dec + target.dec_deg) / 2)))
                    blocked = (np.abs(dec - target.dec_deg) < min_separation) & (ra_sep < min_separation)
                    available.difference_update(np.flatnonzero(blocked))
            rows.append((slot + 1, starts[slot].utc.isot, ends[slot].utc.isot,
                         (starts[slot] + cadence * u.s).utc.isot,
                         name, float(lmsts[slot]), design, offset, altitude, priority, status))
    audit = Table(rows=rows or None, names=('SLOT', 'START_UTC', 'EXPOSURE_END_UTC', 'SLOT_END_UTC',
                                   'OBJECT', 'LMST', 'LMST_DESIGN', 'LMST_OFFSET', 'MIN_ALT', 'PRIORITY', 'STATUS'),
                  dtype=('i4', 'U23', 'U23', 'U23', 'U64', 'f8', 'f8', 'f8', 'f8', 'f8', 'U20'))
    for column in ('LMST', 'LMST_DESIGN', 'LMST_OFFSET', 'MIN_ALT'):
        audit[column].unit = 'deg'
    audit.meta.update(NIGHT=night.date, TWILIGHT=night.twilight,
                      START_UTC=night.start.utc.isot, END_UTC=night.end.utc.isot,
                      LMST_START=night.lmst_start, LMST_END=night.lmst_end,
                      EXPTIME=exptime, OVERHEAD=overhead, LMST_WINDOW=lmst_window,
                      AVAILABLE_SECONDS=night.seconds, CAPACITY=count,
                      SCHEDULED=len(selected), UNUSED_SECONDS=night.seconds-len(selected)*cadence)
    return selected, audit


def describe_night(night, exptime, overhead):
    def clock(degrees):
        seconds = round((degrees % 360) * 240) % 86400
        return f'{seconds//3600:02d}:{seconds//60%60:02d}:{seconds%60:02d}'
    for label, time, lmst in [('Evening', night.start, night.lmst_start),
                               ('Morning', night.end, night.lmst_end)]:
        local = time.to_datetime(timezone=HST).isoformat(timespec='seconds')
        print(f'{label} twilight: {local} HST; UTC {time.utc.isot}; '
              f'LMST {clock(lmst)} ({lmst:.6f} deg)')
    count = slot_count(night.seconds, exptime, overhead)
    print(f'Night: {night.seconds/3600:.6f} hours; capacity {count} targets '
          f'at {exptime:g}+{overhead:g} seconds; '
          f'{night.seconds-count*(exptime+overhead):.3f} seconds after final slot')
