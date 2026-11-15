"""Fetch the SEPA KiWIS snapshot the map is built from (run on a machine that can reach SEPA).

Writes data/sepa/sepa_snapshot.json with, per live gauge:
  v      30 days of 15-minute stage, sampled every 30 minutes  (ts_path */*/SG/15m.Cmd)
  dmean  five years of daily mean stage                        (*/*/SG/HDay.Mean)
  dmax   five years of daily maximum stage                     (*/*/SG/HDay.Max)
The published snapshot was produced by the same queries run in a browser.
"""
import json, sys, datetime as dt, requests
B = 'https://timeseries.sepa.org.uk/KiWIS/KiWIS'
BASE = dict(service='kisters', type='queryServices', datasource=0, format='json')
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'

def ts_list(path):
    r = requests.get(B, params=dict(BASE, request='getTimeseriesList', ts_path=path,
        returnfields='station_no,station_name,catchment_name,station_latitude,station_longitude,ts_id,coverage'))
    return r.json()[1:]
def values(ids, t_from, t_to):
    r = requests.get(B, params=dict(BASE, request='getTimeseriesValues', ts_id=','.join(ids), returnfields='Timestamp,Value',
        metadata='true', md_returnfields='ts_id', **{'from': t_from, 'to': t_to}))
    r.raise_for_status(); return r.json()
def parse(t): return dt.datetime.fromisoformat(t.replace('Z', '+00:00'))

now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0)
t0 = now - dt.timedelta(days=30); step = dt.timedelta(minutes=30); n = 30 * 48 + 1
l15, lmean, lmax = ts_list('*/*/SG/15m.Cmd'), ts_list('*/*/SG/HDay.Mean'), ts_list('*/*/SG/HDay.Max')
live = [x for x in l15 if x[7] and parse(x[7]) > now - dt.timedelta(days=3)]
st = {x[0]: dict(no=x[0], name=x[1], catchment=x[2], lat=float(x[3]), lon=float(x[4]), ts_id=x[5], v=[None] * n) for x in live}
by_ts = {s['ts_id']: s for s in st.values()}
ids = list(by_ts)
for i in range(0, len(ids), 12):
    for o in values(ids[i:i + 12], t0.isoformat(), now.isoformat()):
        s = by_ts[o['ts_id']]
        for t, v in o['data']:
            k, rem = divmod(parse(t) - t0, step)
            if v is not None and not rem and 0 <= k < n:
                s['v'][k] = round(float(v), 3)
d1 = now.date(); d0 = d1 - dt.timedelta(days=1826)
for key, lst in (('dmean', lmean), ('dmax', lmax)):
    m = {x[5]: st[x[0]] for x in lst if x[0] in st}; ids = list(m)
    for i in range(0, len(ids), 10):
        for o in values(ids[i:i + 10], d0.isoformat(), d1.isoformat()):
            a = [None] * 1827
            for t, v in o['data']:
                k = (parse(t).date() - d0).days
                if v is not None and 0 <= k < 1827:
                    a[k] = round(float(v), 3)
            m[o['ts_id']][key] = a
json.dump(dict(generated=dt.datetime.now(dt.timezone.utc).isoformat(), t0=t0.isoformat(), step_min=30, n=n,
               d0=d0.isoformat(), stations=list(st.values())), open(f'{ROOT}/data/sepa/sepa_snapshot.json', 'w'))
