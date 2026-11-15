"""Pack telemetry + baseline + network into the single data.json the page loads.

Normalisation knots per gauge (metres -> 0..1), from five years of SEPA daily series:
  0.00  1st percentile of daily mean stage
  0.32  median daily mean
  0.62  90th percentile daily mean
  0.84  99th percentile daily mean
  1.00  highest daily maximum on record in the window (and never below the 30-day peak)
"""
import json, sys, base64, numpy as np
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
snap = json.load(open(f'{ROOT}/data/sepa/sepa_snapshot.json'))
net = json.load(open(f'{ROOT}/build/network.json'))
S = snap['stations']; step = 60 // snap['step_min']

def arr(a): return np.array([np.nan if x is None else x for x in a], float)
def fill(v):
    ok = np.isfinite(v); i = np.arange(len(v))
    return np.interp(i, i[ok], v[ok])

stations = []; rows = []; lrows = []; llo = []; lhi = []
for s, n in zip(S, net['stations']):
    v = fill(arr(s['v']))[::step]
    dm = arr(s['dmean']); dx = arr(s['dmax'])
    dm = dm[np.isfinite(dm)]; dx = dx[np.isfinite(dx)]
    p1, p50, p90, p99 = np.percentile(dm, [1, 50, 90, 99])
    top = max(dx.max() if len(dx) else p99, v.max(), p99 + 0.01)
    k = [min(p1, v.min()), p50, p90, p99, top]
    for i in range(1, 5):                        # strictly increasing
        k[i] = max(k[i], k[i - 1] + 0.002)
    lo, hi = v.min(), max(v.max(), v.min() + 0.002)
    rows.append(np.round(254 * (v - lo) / (hi - lo)).astype(np.uint8))
    # long record: five years of daily means; 255 marks a gap longer than a week
    d = arr(s['dmean']); ok = np.isfinite(d); i = np.arange(len(d))
    near = np.abs(i - i[ok][np.clip(np.searchsorted(i[ok], i), 0, ok.sum() - 1)])
    prev = np.abs(i - i[ok][np.clip(np.searchsorted(i[ok], i) - 1, 0, ok.sum() - 1)])
    gap = np.minimum(near, prev) > 7
    df = np.interp(i, i[ok], d[ok]); a, b = d[ok].min(), max(d[ok].max(), d[ok].min() + 0.002)
    q = np.round(254 * (df - a) / (b - a)).astype(np.uint8); q[gap] = 255
    lrows.append(q); llo.append(round(float(a), 3)); lhi.append(round(float(b), 3))
    stations.append({'no': s['no'], 'n': s['name'], 'c': s['catchment'], 'r': n['river'],
                     'x': n['x'], 'y': n['y'], 's': int(n['snapped']), 'u': n['up_km'],
                     'lo': round(lo, 3), 'hi': round(hi, 3), 'k': [round(x, 3) for x in k],
                     'yrs': round(len(dm) / 365.25, 1)})
V = np.stack(rows); L = np.stack(lrows)
out = {'generated': snap['generated'], 't0': snap['t0'], 'stepMin': 60, 'frames': int(V.shape[1]),
       'extent': net['extent'], 'q': net['q'], 'stations': stations,
       'series': base64.b64encode(V.tobytes()).decode(), 'groups': net['groups'], 'labels': net['labels'],
       'long': {'t0': snap['d0'] + 'T00:00:00Z', 'stepMin': 1440, 'frames': int(L.shape[1]), 'lo': llo, 'hi': lhi,
                'series': base64.b64encode(L.tobytes()).decode()}}
json.dump(out, open(f'{ROOT}/build/data.json', 'w'), separators=(',', ':'))
print('stations', len(stations), 'frames', V.shape[1], 'bytes', len(json.dumps(out)))

# diagnostics: does the pulse travel west -> east?
def norm(s, v): return np.interp(v, s['k'], [0, .32, .62, .84, 1])
N = np.stack([norm(s, s['lo'] + r / 254 * (s['hi'] - s['lo'])) for s, r in zip(stations, V)])
x = np.array([s['x'] for s in stations]) * net['q'] / 1000
pk = N.argmax(1); strong = N.max(1) > 0.62
print('gauges peaking above 0.62:', strong.sum(), ' >0.84:', (N.max(1) > .84).sum(), ' =1.0:', (N.max(1) >= .999).sum())
print('corr(easting, peak hour):', round(float(np.corrcoef(x[strong], pk[strong])[0, 1]), 2))
for a, b in ((0, 150), (150, 220), (220, 290), (290, 360), (360, 440)):
    m = strong & (x >= a) & (x < b)
    print(f'  easting {a}-{b} km: n={m.sum()} median peak hour {np.median(pk[m]) if m.any() else None}')
print('national mean by day:', np.round(N.mean(0)[::24], 2).tolist())
