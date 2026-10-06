"""Human-readable night outputs: a plain-text timeline and a sky-plot PDF.

Requires numpy; the plot additionally requires matplotlib (imported lazily).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


# Plot limits requested for a "looking up" view: RA decreases to the right.
PLOT_RA_LEFT = 300.0
PLOT_RA_RIGHT = -60.0
PLOT_DEC_RANGE = (-20.0, 20.0)


def lmst_hms(lmst_deg: float) -> str:
    """Format a sidereal time in degrees as hh:mm:ss."""
    total = int(round((lmst_deg % 360.0) / 15.0 * 3600.0)) % 86400
    hours, remainder = divmod(total, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f'{hours:02d}:{minutes:02d}:{seconds:02d}'


def write_schedule_text(path: Path, audit) -> None:
    """Write one line per slot: exposure-start UT, LMST at midpoint, object name."""
    meta = audit.meta
    lines = [
        f'# Night of {meta["NIGHT"]} (Hawaii date); {meta["TWILIGHT"]:g}-degree twilight',
        f'# Twilight UT: {meta["START_UTC"][:19]} to {meta["END_UTC"][:19]}',
        f'# LMST at twilight: {lmst_hms(meta["LMST_START"])} to {lmst_hms(meta["LMST_END"])}',
        f'# Exposure {meta["EXPTIME"]:g} s + overhead {meta["OVERHEAD"]:g} s; '
        f'{meta["SCHEDULED"]} of {meta["CAPACITY"]} slots filled',
        f'# Moon: {100 * meta["MOON_ILLUMINATION"]:.0f}% illuminated at mid-night'
        if 'MOON_ILLUMINATION' in meta else '# Moon: illumination not computed',
        '# UT is the exposure start; LMST is at the exposure midpoint.',
        f'# {"UT":19s}  {"LMST":8s}  {"LMST_DEG":8s}  OBJECT',
    ]
    for row in audit:
        name = row['OBJECT'] or f'({row["STATUS"]})'
        ut = str(row['START_UTC'])[:19].replace('T', ' ')
        lines.append(f'{ut}  {lmst_hms(row["LMST"])}  {row["LMST"]:8.3f}  {name}')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def read_tile_positions(ecsv_path: Path) -> dict[str, np.ndarray]:
    """Read RA, DEC, FILTER, and flag columns from the tile catalog without astropy."""
    wanted = ('RA', 'DEC', 'FILTER', 'IN_IBIS', 'IN_HSC', 'DONE')
    values = {name: [] for name in wanted}
    with ecsv_path.open(newline='', encoding='utf-8') as handle:
        rows = csv.reader((line for line in handle if not line.startswith('#')),
                          delimiter=' ', skipinitialspace=True)
        header = next(rows)
        index = {name: header.index(name) for name in wanted}
        for parts in rows:
            if not parts:
                continue
            for name in wanted:
                values[name].append(parts[index[name]])
    out = {
        'RA': np.asarray(values['RA'], dtype=float),
        'DEC': np.asarray(values['DEC'], dtype=float),
        'FILTER': np.asarray([v.upper() for v in values['FILTER']]),
    }
    for name in ('IN_IBIS', 'IN_HSC', 'DONE'):
        out[name] = np.asarray(values[name], dtype=float).astype(int)
    return out


def plot_ra(ra_deg):
    """Map RA so that the plotted x range runs continuously from 300 down to -60."""
    ra = np.asarray(ra_deg, dtype=float) % 360.0
    return np.where(ra > PLOT_RA_LEFT, ra - 360.0, ra)


BODY_STYLE = {
    'moon': dict(marker='o', color='0.3', size=55, label='Moon'),
    'venus': dict(marker='*', color='gold', size=110, label='Venus'),
    'mars': dict(marker='*', color='orangered', size=90, label='Mars'),
    'saturn': dict(marker='*', color='tan', size=90, label='Saturn'),
    'uranus': dict(marker='*', color='turquoise', size=70, label='Uranus'),
}


def plot_night(path: Path, audit, tiles: dict[str, np.ndarray], targets,
               filter_name: str, bodies: dict | None = None) -> None:
    """Draw tonight's targets over the tile footprint in a rectilinear RA/DEC view.

    ``bodies`` maps body name to (ra_deg, dec_deg) arrays sampled through the
    night, as returned by night_planning.bright_body_positions.
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit('The night sky plot requires matplotlib; '
                         'install requirements.txt in your Python environment') from exc

    meta = audit.meta
    in_filter = tiles['FILTER'] == filter_name.upper()
    footprint = in_filter & (tiles['IN_IBIS'] == 1)
    done = footprint & (tiles['DONE'] != 0)
    hsc = footprint & (tiles['IN_HSC'] == 1) & ~done
    remaining = footprint & ~done & ~hsc

    fig, ax = plt.subplots(figsize=(14, 3.8))
    x_all = plot_ra(tiles['RA'])
    for mask, color, label in ((remaining, '0.82', 'Footprint: remaining tiles'),
                               (hsc, '#c9d8ea', 'Footprint: IN_HSC tiles'),
                               (done, '0.45', 'Footprint: DONE tiles')):
        if np.any(mask):
            ax.scatter(x_all[mask], tiles['DEC'][mask], s=4, c=color, marker='s',
                       linewidths=0, label=f'{label} ({int(mask.sum())})', rasterized=True)

    # Tonight's targets, colored by the requested MAG_AB of each exposure.
    scheduled = audit[audit['STATUS'] == 'SCHEDULED']
    position = {t.object_name: (t.ra_deg, t.dec_deg) for t in targets}
    if len(scheduled):
        ra = np.asarray([position[name][0] for name in scheduled['OBJECT']])
        dec = np.asarray([position[name][1] for name in scheduled['OBJECT']])
        points = ax.scatter(plot_ra(ra), dec, s=22, c=np.asarray(scheduled['MAG_AB']),
                            cmap='viridis', edgecolors='black', linewidths=0.3, zorder=3,
                            label=f'Tonight ({len(scheduled)} exposures)')
        bar = fig.colorbar(points, ax=ax, pad=0.01, fraction=0.03, shrink=0.8)
        bar.set_label('Requested MAG_AB')

    for lmst, label in ((meta['LMST_START'], 'meridian at evening twilight'),
                        (meta['LMST_END'], 'meridian at morning twilight')):
        x = float(plot_ra(lmst))
        on_plot = PLOT_RA_RIGHT <= x <= PLOT_RA_LEFT
        if on_plot:
            ax.axvline(x, color='tab:red', linestyle='--', linewidth=0.8, zorder=2)
        # Keep the label inside the axes even when the meridian falls off the plot.
        x_text = min(max(x, PLOT_RA_RIGHT + 1.0), PLOT_RA_LEFT - 1.0)
        near_right = x_text < (PLOT_RA_LEFT + PLOT_RA_RIGHT) / 2
        ax.text(x_text, PLOT_DEC_RANGE[1] - 0.5,
                f'{label}{"" if on_plot else " (off plot)"}\nLMST {lmst_hms(lmst)}',
                color='tab:red', fontsize=7, va='top', ha='right' if near_right else 'left')

    # Planets are marked at the midpoint of the night (the sampled tracks are
    # centered on it). The Moon moves several degrees per night, so it is
    # drawn at evening and morning twilight with an arrow showing its motion.
    for body, (ra, dec) in (bodies or {}).items():
        style = BODY_STYLE.get(body, dict(marker='*', color='gray', size=70, label=body))
        x, y = plot_ra(ra), np.asarray(dec, dtype=float)
        middle = len(x) // 2
        in_view = ((PLOT_RA_RIGHT <= x) & (x <= PLOT_RA_LEFT)
                   & (PLOT_DEC_RANGE[0] <= y) & (y <= PLOT_DEC_RANGE[1]))
        if body == 'moon':
            if not (in_view[0] or in_view[-1]):
                continue
            # Break the track where it wraps across the plot edge.
            segments = np.split(np.arange(len(x)), np.flatnonzero(np.abs(np.diff(x)) > 180) + 1)
            for segment in segments:
                ax.plot(x[segment], y[segment], color=style['color'], linewidth=1.2, zorder=4)
            ax.scatter(x[0], y[0], marker='o', s=style['size'], facecolors='white',
                       edgecolors=style['color'], linewidths=1.2, zorder=5)
            ax.scatter(x[-1], y[-1], marker='o', s=style['size'], color=style['color'],
                       edgecolors='black', linewidths=0.5, zorder=5)
            if abs(x[-1] - x[-2]) < 180:
                ax.annotate('', xy=(x[-1], y[-1]), xytext=(x[-2], y[-2]), zorder=5,
                            arrowprops=dict(arrowstyle='-|>', color=style['color'], lw=1.2,
                                            shrinkA=0, shrinkB=6))
            label_x, label_y = (x[-1], y[-1]) if in_view[-1] else (x[0], y[0])
            label = 'Moon (evening \u2192 morning)'
            if 'MOON_ILLUMINATION' in meta:
                label = f'{100 * meta["MOON_ILLUMINATION"]:.0f}% ' + label
        else:
            if not in_view[middle]:
                continue
            ax.scatter(x[middle], y[middle], marker=style['marker'], s=style['size'],
                       color=style['color'], edgecolors='black', linewidths=0.5, zorder=5)
            label_x, label_y = x[middle], y[middle]
            label = style['label']
        near_top = label_y > PLOT_DEC_RANGE[1] - 3.0
        ax.annotate(label, (label_x, label_y), xytext=(6, -10 if near_top else 5),
                    textcoords='offset points', fontsize=8, zorder=5)

    ax.set_xlim(PLOT_RA_LEFT, PLOT_RA_RIGHT)
    ax.set_ylim(*PLOT_DEC_RANGE)
    ticks = np.arange(PLOT_RA_LEFT, PLOT_RA_RIGHT - 1, -30.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(t % 360)}' for t in ticks])
    ax.set_xlabel('RA (deg)')
    ax.set_ylabel('DEC (deg)')
    ax.set_aspect('equal')
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.set_title(f'CFHT {filter_name.upper()} plan for the night of {meta["NIGHT"]} (HST): '
                 f'{meta["SCHEDULED"]} of {meta["CAPACITY"]} slots, '
                 f'UT {meta["START_UTC"][11:16]} to {meta["END_UTC"][11:16]}; '
                 'planets at mid-night', fontsize=10)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.38), ncol=4, fontsize=7,
              markerscale=2, frameon=False)
    fig.savefig(path, format='pdf', dpi=200, bbox_inches='tight')
    plt.close(fig)
