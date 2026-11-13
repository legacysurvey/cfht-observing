alias_run_ids = dict(
    NAOC = '26BZ01',
    LBNL = '26BZ50',
)

baseurl = 'https://api.cfht.hawaii.edu/'

# from astrometry.util.starutil_numpy.py
def timedeltatodays(dt):
    return dt.days + (dt.seconds + dt.microseconds/1e6)/86400.

def datetomjd(d):
    import datetime
    d0 = datetime.datetime(1858, 11, 17, 0, 0, 0, tzinfo=datetime.timezone.utc)
    dt = d - d0
    # dt is a timedelta object.
    return timedeltatodays(dt)

def mjdtodate(mjd):
    jd = mjdtojd(mjd)
    return jdtodate(jd)

def jdtodate(jd):
    import datetime
    unixtime = (jd - 2440587.5) * 86400. # in seconds
    return datetime.datetime.utcfromtimestamp(unixtime)

def mjdtojd(mjd):
    return mjd + 2400000.5


# Which exposures count as observing a tile.
SNR_GOAL = 4.6            # intended SNR goal of all templates so far
SNR_MARGIN = 0.4          # an exposure is good when SNR > SNR_GOAL - SNR_MARGIN
BAD_GRADE = 5             # CFHT grade that marks an exposure as bad
SNAP_COMMENT = 'MEGACAM SNR SNAP'


def exposure_quality(exposure_status, snr_goal=SNR_GOAL):
    """Classify one Kealahou exposure_status record.

    Returns (quality, snr): quality is
      'good'     graded, not marked bad (grade 5), and SNR > snr_goal - 0.4;
      'low-snr'  validated by CFHT but with lower SNR: counted as observed for
                 now (CFHT has charged it), to be topped up later;
      None       does not count: a snapshot, not graded yet, marked bad, or
                 neither good nor validated.
    SNR is the Elixir value in megacam_status (None if absent).
    """
    st = exposure_status or {}
    if st.get('obs_comment') == SNAP_COMMENT:
        return None, None
    snr = ((st.get('megacam_status') or {}).get('elixir_processing_result') or {}).get('snr')
    grade = st.get('grade')
    if grade is None or grade == BAD_GRADE:
        return None, snr
    if snr is not None and snr > round(snr_goal - SNR_MARGIN, 6):
        return 'good', snr
    if st.get('exp_status') == 'VALIDATED':
        return 'low-snr', snr
    return None, snr
