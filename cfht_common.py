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
