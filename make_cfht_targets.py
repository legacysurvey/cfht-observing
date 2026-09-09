#!/usr/bin/env python3
"""Generate CFHT Kealohou MegaCam fixed-target lists.

The default input is ``obstatus/cfht-tiles.ecsv``.  The script writes XML in
the same ASTRO/CSV format as ``obstatus/megacam_fixed_target.xml`` and a FITS
binary table with the same target-list columns plus decimal-degree RA/DEC.
Only unfinished tiles with IN_IBIS=1 and IN_HSC!=1 are selected. PROGRAM
defaults to LBNL; --all-programs removes that restriction. --night builds
a twilight-to-twilight queue in observing order (requires numpy/astropy).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import struct
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


DEFAULT_INPUT = Path("obstatus/cfht-tiles.ecsv")
DEFAULT_TEMPLATE = Path("obstatus/megacam_fixed_target.xml")
DEFAULT_OUTDIR = Path("plans")
DEFAULT_FILTER = "M4376"
DEFAULT_MAG_AB = 24.25
DEFAULT_MIN_SEPARATION_DEG = 1.0

REQUIRED_COLUMNS = ("OBJECT", "RA", "DEC", "FILTER", "IN_IBIS", "IN_HSC", "DONE", "PRIORITY")

XML_TABLE_HEADER = [
    "NAME                                   |RA_J2000   |DEC_J2000   |MAG_AB|PM_RA |PM_DEC|POINT_RA|POINT_DEC|",
    "                                       |hh:mm:ss.ss|+dd:mm:ss.ss|      |      |      |+mm:ss.s|+mm:ss.s|",
    "123456789012345678901234567890123456789|12345678901|123456789012|123456|123456|123456| 1234567|12345678|",
    "---------------------------------------|-----------|------------|------|------|------|--------|--------|",
]


@dataclass(frozen=True)
class TileTarget:
    object_name: str
    ra_deg: float
    dec_deg: float
    lmst_design_deg: float | None = None
    priority: float = 1.0


@dataclass(frozen=True)
class OutputTarget:
    name: str
    ra_j2000: str
    dec_j2000: str
    ra_deg: float
    dec_deg: float
    mag_ab: float
    pm_ra: float = 0.0
    pm_dec: float = 0.0
    point_ra: str = "+00:00.0"
    point_dec: str = "+00:00.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate CFHT Kealohou MegaCam fixed-target XML and FITS target lists."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input CFHT tile ECSV file (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--template",
        type=Path,
        default=DEFAULT_TEMPLATE,
        help=f"XML template whose non-data header/footer should be retained (default: {DEFAULT_TEMPLATE}).",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help=f"Directory for generated plans (default: {DEFAULT_OUTDIR}).",
    )
    parser.add_argument(
        "--date",
        default=dt.date.today().isoformat(),
        help="Date string used in the default output prefix (default: today).",
    )
    parser.add_argument('--night', metavar='YYYY-MM-DD',
                        help='Schedule the night beginning on this Hawaii calendar date.')
    parser.add_argument('--twilight', type=float, default=15.0, metavar='DEG',
                        help='Solar depression defining night (default: 15 degrees).')
    parser.add_argument('--exptime', type=float, default=130.0, metavar='SECONDS',
                        help='Exposure time per night target (default: 130 seconds).')
    parser.add_argument('--overhead', type=float, default=44.0, metavar='SECONDS',
                        help='Overhead after each exposure (default: 44 seconds).')
    parser.add_argument('--lmst-window', type=float, default=5.0, metavar='DEG',
                        help='Maximum design LMST mismatch in night mode (default: 5 degrees).')
    parser.add_argument(
        "--prefix",
        help="Output filename prefix (default: targets-$NIGHT, or targets-$DATE).",
    )
    parser.add_argument(
        "--ra",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(0.0, 360.0),
        help="RA range in degrees. Use MIN > MAX for a wraparound range (default: 0 360).",
    )
    parser.add_argument(
        "--dec",
        nargs=2,
        type=float,
        metavar=("MIN_DEG", "MAX_DEG"),
        default=(-90.0, 90.0),
        help="DEC range in degrees (default: -90 90).",
    )
    parser.add_argument(
        "--filter",
        dest="filter_name",
        default=DEFAULT_FILTER,
        help=f"Filter name to select (default: {DEFAULT_FILTER}).",
    )
    program_group = parser.add_mutually_exclusive_group()
    program_group.add_argument(
        "--program",
        dest="programs",
        nargs="+",
        action="extend",
        metavar="PROGRAM",
        help=(
            "Select one or more programs, e.g. --program LBNL NAOC. "
            "Matches ignore case. "
            "Use --program '' for blank labels. "
            "Default: LBNL only."
        ),
    )
    program_group.add_argument('--all-programs', action='store_true',
                               help='Include every program, including blank labels.')
    parser.add_argument(
        "--mag-ab",
        type=float,
        default=DEFAULT_MAG_AB,
        help=f"MAG_AB value to write for every target (default: {DEFAULT_MAG_AB}).",
    )
    parser.add_argument(
        "--non-overlapping",
        action="store_true",
        help="Greedily trim to non-overlapping 1-degree footprints.",
    )
    parser.add_argument(
        "--min-separation",
        type=float,
        default=DEFAULT_MIN_SEPARATION_DEG,
        help=(
            "Footprint size in degrees for --non-overlapping; candidates are rejected "
            "when both DEC and cos(DEC)-projected RA separations are smaller than this "
            f"value (default: {DEFAULT_MIN_SEPARATION_DEG})."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )

    args = parser.parse_args()
    if args.programs is None and not args.all_programs:
        args.programs = ['LBNL']
    validate_args(parser, args)
    return args


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    dec_min, dec_max = args.dec
    if dec_min > dec_max:
        parser.error("--dec MIN_DEG must be less than or equal to MAX_DEG")
    if dec_min < -90.0 or dec_max > 90.0:
        parser.error("--dec values must be within [-90, 90] degrees")
    if args.min_separation <= 0.0:
        parser.error("--min-separation must be positive")
    if not args.filter_name:
        parser.error("--filter must not be empty")
    if not args.date:
        parser.error("--date must not be empty")
    if args.night:
        try:
            dt.date.fromisoformat(args.night)
        except ValueError:
            parser.error('--night must be a valid YYYY-MM-DD calendar date')
    if not math.isfinite(args.twilight) or not 0 < args.twilight < 90:
        parser.error('--twilight must be between 0 and 90 degrees')
    if not math.isfinite(args.exptime) or args.exptime <= 0:
        parser.error('--exptime must be finite and positive')
    if not math.isfinite(args.overhead) or args.overhead < 0:
        parser.error('--overhead must be finite and nonnegative')
    if not math.isfinite(args.lmst_window) or not 0 < args.lmst_window <= 180:
        parser.error('--lmst-window must be in (0, 180] degrees')


def normalize_program(program: str) -> str:
    """Normalize whitespace and case for program matching."""
    return program.strip().upper()


def iter_matching_targets(
    ecsv_path: Path,
    ra_range: Sequence[float],
    dec_range: Sequence[float],
    filter_name: str,
    programs: Sequence[str] | None = ("LBNL",),
    require_lmst: bool = False,
) -> Iterator[TileTarget]:
    """Stream targets from the CFHT ECSV file that match the observing cuts."""

    normalized_filter = filter_name.upper()
    normalized_programs = (
        {normalize_program(program) for program in programs}
        if programs is not None else None
    )
    required_columns = REQUIRED_COLUMNS + (
        ("PROGRAM",) if normalized_programs is not None else ()
    ) + (("LMST_DESIGN",) if require_lmst else ())
    columns = None
    column_index = None

    with ecsv_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            # ECSV uses CSV quoting for empty strings and strings with spaces.
            try:
                parts = next(csv.reader([stripped], delimiter=" ",
                                        skipinitialspace=True, strict=True))
            except csv.Error as exc:
                raise ValueError(f"Could not parse {ecsv_path}:{line_number}: {stripped}") from exc
            if columns is None:
                columns = parts
                column_index = {name: index for index, name in enumerate(columns)}
                missing = [name for name in required_columns if name not in column_index]
                if missing:
                    raise ValueError(
                        f"{ecsv_path} is missing required column(s): {', '.join(missing)}"
                    )
                continue

            if len(parts) < len(columns):
                raise ValueError(
                    f"{ecsv_path}:{line_number} has {len(parts)} fields, expected {len(columns)}"
                )

            assert column_index is not None
            try:
                object_name = parts[column_index["OBJECT"]]
                ra_deg = float(parts[column_index["RA"]])
                dec_deg = float(parts[column_index["DEC"]])
                row_filter = parts[column_index["FILTER"]].upper()
                in_ibis = int(parts[column_index["IN_IBIS"]])
                in_hsc = int(parts[column_index["IN_HSC"]])
                done = int(parts[column_index["DONE"]])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"Could not parse {ecsv_path}:{line_number}: {stripped}") from exc

            if row_filter != normalized_filter:
                continue
            if in_ibis != 1 or in_hsc == 1 or done != 0:
                continue
            if (normalized_programs is not None
                    and normalize_program(parts[column_index["PROGRAM"]])
                    not in normalized_programs):
                continue
            if not in_ra_range(ra_deg, ra_range):
                continue
            if not in_dec_range(dec_deg, dec_range):
                continue

            try:
                priority = float(parts[column_index['PRIORITY']])
            except ValueError as exc:
                raise ValueError(f'Invalid PRIORITY for {object_name}') from exc
            if not math.isfinite(priority):
                raise ValueError(f'{object_name} needs a finite PRIORITY; '
                                 'run update_tile_priorities.py')
            lmst = None
            if require_lmst:
                try:
                    lmst = float(parts[column_index['LMST_DESIGN']])
                except ValueError as exc:
                    raise ValueError(f'Invalid LMST_DESIGN for {object_name}') from exc
                if not math.isfinite(lmst) or not 0 <= lmst < 360:
                    raise ValueError(f'{object_name} needs a valid LMST_DESIGN; '
                                     'run update_lmst_design.py')
            yield TileTarget(object_name=object_name, ra_deg=ra_deg, dec_deg=dec_deg,
                             lmst_design_deg=lmst, priority=priority)

    if columns is None:
        raise ValueError(f"{ecsv_path} does not contain a data header")


def in_ra_range(ra_deg: float, ra_range: Sequence[float]) -> bool:
    ra_min, ra_max = ra_range
    if ra_max - ra_min >= 360.0:
        return True

    ra = ra_deg % 360.0
    lo = ra_min % 360.0
    hi = ra_max % 360.0
    if lo <= hi:
        return lo <= ra <= hi
    return ra >= lo or ra <= hi


def in_dec_range(dec_deg: float, dec_range: Sequence[float]) -> bool:
    dec_min, dec_max = dec_range
    return dec_min <= dec_deg <= dec_max


def trim_non_overlapping(
    targets: Iterable[TileTarget], min_separation_deg: float
) -> list[TileTarget]:
    """Greedily remove targets whose projected 1-degree footprints overlap."""

    selected: list[TileTarget] = []
    dec_bins: dict[int, list[TileTarget]] = defaultdict(list)

    for target in targets:
        dec_bin = math.floor(target.dec_deg / min_separation_deg)
        keep = True

        for neighbor_bin in range(dec_bin - 1, dec_bin + 2):
            for other in dec_bins.get(neighbor_bin, []):
                dec_sep = abs(target.dec_deg - other.dec_deg)
                if dec_sep >= min_separation_deg:
                    continue

                mean_dec_rad = math.radians(0.5 * (target.dec_deg + other.dec_deg))
                projected_ra_sep = angular_ra_separation(
                    target.ra_deg, other.ra_deg
                ) * abs(math.cos(mean_dec_rad))
                if projected_ra_sep < min_separation_deg:
                    keep = False
                    break
            if not keep:
                break

        if keep:
            selected.append(target)
            dec_bins[dec_bin].append(target)

    return selected


def angular_ra_separation(ra1_deg: float, ra2_deg: float) -> float:
    return abs((ra1_deg - ra2_deg + 180.0) % 360.0 - 180.0)


def make_output_targets(targets: Iterable[TileTarget], mag_ab: float) -> list[OutputTarget]:
    return [
        OutputTarget(
            name=target.object_name,
            ra_j2000=format_ra(target.ra_deg),
            dec_j2000=format_dec(target.dec_deg),
            ra_deg=target.ra_deg,
            dec_deg=target.dec_deg,
            mag_ab=mag_ab,
        )
        for target in targets
    ]


def format_ra(ra_deg: float) -> str:
    total_centiseconds = round(((ra_deg % 360.0) / 15.0) * 3600.0 * 100.0)
    total_centiseconds %= 24 * 3600 * 100

    hours, remainder = divmod(total_centiseconds, 3600 * 100)
    minutes, remainder = divmod(remainder, 60 * 100)
    seconds, centiseconds = divmod(remainder, 100)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}.{centiseconds:02d}"


def format_dec(dec_deg: float) -> str:
    sign = "-" if math.copysign(1.0, dec_deg) < 0.0 else "+"
    total_centiseconds = round(abs(dec_deg) * 3600.0 * 100.0)

    degrees, remainder = divmod(total_centiseconds, 3600 * 100)
    minutes, remainder = divmod(remainder, 60 * 100)
    seconds, centiseconds = divmod(remainder, 100)
    return f"{sign}{degrees:02d}:{minutes:02d}:{seconds:02d}.{centiseconds:02d}"


def write_xml(template_path: Path, output_path: Path, rows: Sequence[OutputTarget]) -> None:
    template = template_path.read_text(encoding="utf-8")
    cdata_marker = "<![CDATA[\n"
    data_end_marker = "]]></CSV></DATA>"
    try:
        prefix, remainder = template.split(cdata_marker, 1)
        _old_data, suffix = remainder.split(data_end_marker, 1)
    except ValueError as exc:
        raise ValueError(
            f"{template_path} does not look like the expected MegaCam target XML template"
        ) from exc

    data_lines = XML_TABLE_HEADER + [format_xml_row(row) for row in rows]
    output = prefix + cdata_marker + "\n".join(data_lines) + "\n" + data_end_marker + suffix
    output_path.write_text(output, encoding="utf-8")


def format_xml_row(row: OutputTarget) -> str:
    return (
        f"{row.name:<39}|"
        f"{row.ra_j2000:11}|"
        f"{row.dec_j2000:12}|"
        f"{row.mag_ab:6.2f}|"
        f"{row.pm_ra:6.1f}|"
        f"{row.pm_dec:6.1f}|"
        f"{row.point_ra:8}|"
        f"{row.point_dec:8}|"
    )


def write_fits(output_path: Path, rows: Sequence[OutputTarget]) -> None:
    name_width = max([39] + [len(row.name) for row in rows])
    columns = [
        ("NAME", f"{name_width}A", "str", name_width),
        ("RA_J2000", "11A", "str", 11),
        ("DEC_J2000", "12A", "str", 12),
        ("RA", "D", "float64", 8),
        ("DEC", "D", "float64", 8),
        ("MAG_AB", "E", "float32", 4),
        ("PM_RA", "E", "float32", 4),
        ("PM_DEC", "E", "float32", 4),
        ("POINT_RA", "8A", "str", 8),
        ("POINT_DEC", "8A", "str", 8),
    ]
    row_length = sum(column[3] for column in columns)

    with output_path.open("wb") as handle:
        write_fits_header(
            handle,
            [
                fits_card("SIMPLE", True),
                fits_card("BITPIX", 8),
                fits_card("NAXIS", 0),
                fits_card("EXTEND", True),
            ],
        )

        cards = [
            fits_card("XTENSION", "BINTABLE"),
            fits_card("BITPIX", 8),
            fits_card("NAXIS", 2),
            fits_card("NAXIS1", row_length),
            fits_card("NAXIS2", len(rows)),
            fits_card("PCOUNT", 0),
            fits_card("GCOUNT", 1),
            fits_card("TFIELDS", len(columns)),
            fits_card("EXTNAME", "TARGETS"),
        ]
        for index, (name, form, _kind, _width) in enumerate(columns, start=1):
            cards.append(fits_card(f"TTYPE{index}", name))
            cards.append(fits_card(f"TFORM{index}", form))
        write_fits_header(handle, cards)

        for row in rows:
            handle.write(pack_fits_row(row, columns))

        pad = (-len(rows) * row_length) % 2880
        if pad:
            handle.write(b"\0" * pad)


def pack_fits_row(row: OutputTarget, columns: Sequence[tuple[str, str, str, int]]) -> bytes:
    values = {
        "NAME": row.name,
        "RA_J2000": row.ra_j2000,
        "DEC_J2000": row.dec_j2000,
        "RA": row.ra_deg,
        "DEC": row.dec_deg,
        "MAG_AB": row.mag_ab,
        "PM_RA": row.pm_ra,
        "PM_DEC": row.pm_dec,
        "POINT_RA": row.point_ra,
        "POINT_DEC": row.point_dec,
    }

    packed = bytearray()
    for name, _form, kind, width in columns:
        value = values[name]
        if kind == "str":
            packed.extend(fits_string(str(value), width))
        elif kind == "float32":
            packed.extend(struct.pack(">f", float(value)))
        elif kind == "float64":
            packed.extend(struct.pack(">d", float(value)))
        else:
            raise ValueError(f"Unsupported FITS column kind: {kind}")
    return bytes(packed)


def fits_string(value: str, width: int) -> bytes:
    encoded = value.encode("ascii", errors="replace")[:width]
    return encoded.ljust(width, b" ")


def write_fits_header(handle, cards: Sequence[bytes]) -> None:
    header = b"".join(cards) + b"END".ljust(80)
    pad = (-len(header)) % 2880
    handle.write(header + b" " * pad)


def fits_card(keyword: str, value) -> bytes:
    if len(keyword) > 8:
        raise ValueError(f"FITS keyword is longer than 8 characters: {keyword}")

    if isinstance(value, bool):
        value_text = "T" if value else "F"
        text = f"{keyword:<8}= {value_text:>20}"
    elif isinstance(value, int):
        text = f"{keyword:<8}= {value:>20d}"
    elif isinstance(value, float):
        text = f"{keyword:<8}= {value:>20.10G}"
    else:
        escaped = str(value).replace("'", "''")
        text = f"{keyword:<8}= '{escaped}'"

    return text[:80].ljust(80).encode("ascii")


def output_paths(outdir: Path, prefix: str, overwrite: bool) -> tuple[Path, Path]:
    xml_path = outdir / f"{prefix}.xml"
    fits_path = outdir / f"{prefix}.fits"
    if not overwrite:
        existing = [str(path) for path in (xml_path, fits_path) if path.exists()]
        if existing:
            raise FileExistsError(
                "Output file(s) already exist; use --overwrite to replace them: "
                + ", ".join(existing)
            )
    return xml_path, fits_path


def main() -> None:
    args = parse_args()
    prefix = args.prefix or f"targets-{args.night or args.date}"
    xml_path, fits_path = output_paths(args.outdir, prefix, args.overwrite)
    audit_path = args.outdir / f'{prefix}.schedule.ecsv'
    if args.night and audit_path.exists() and not args.overwrite:
        raise FileExistsError(f'{audit_path} already exists; use --overwrite to replace it')

    targets = sorted(
        iter_matching_targets(
            args.input,
            ra_range=args.ra,
            dec_range=args.dec,
            filter_name=args.filter_name,
            programs=args.programs,
            require_lmst=bool(args.night),
        ),
        key=lambda target: target.ra_deg,
    )
    before_non_overlap = len(targets)

    audit = None
    if args.night:
        try:
            from night_planning import describe_night, night_bounds, schedule_night
        except ImportError as exc:
            raise SystemExit('Night planning requires numpy and astropy; '
                             'install requirements.txt in your Python environment') from exc
        night = night_bounds(args.night, args.twilight)
        describe_night(night, args.exptime, args.overhead)
        targets, audit = schedule_night(targets, night, args.exptime, args.overhead,
                                        args.lmst_window, args.non_overlapping,
                                        args.min_separation)
    elif args.non_overlapping:
        targets = trim_non_overlapping(
            sorted(targets, key=lambda t: (-t.priority, t.dec_deg, t.ra_deg, t.object_name)),
            args.min_separation)
        targets.sort(key=lambda t: t.ra_deg)

    output_rows = make_output_targets(targets, args.mag_ab)

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_xml(args.template, xml_path, output_rows)
    write_fits(fits_path, output_rows)
    if audit is not None:
        audit.write(audit_path, format='ascii.ecsv', overwrite=args.overwrite)

    print(f"Read: {args.input}")
    print(f"Selected targets: {before_non_overlap}")
    if audit is not None:
        missing = len(audit) - len(output_rows)
        print(f'Scheduled targets: {len(output_rows)} / {len(audit)} slots; '
              f'unfilled slots: {missing}; unused time: {audit.meta["UNUSED_SECONDS"]:.3f} seconds')
        if missing:
            print('Night is not full: inspect the schedule for slots without eligible nearby targets.')
        print(f'Wrote: {audit_path}')
    elif args.non_overlapping:
        print(f"After non-overlap trim: {len(output_rows)}")
    print(f"Wrote: {xml_path}")
    print(f"Wrote: {fits_path}")


if __name__ == "__main__":
    main()
