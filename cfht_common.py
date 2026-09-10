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
    d0 = datetime.datetime(1858, 11, 17, 0, 0, 0, tzinfo=datetime.UTC)
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
