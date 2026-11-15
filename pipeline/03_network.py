"""Snap SEPA gauges to OS Open Rivers and assign every drawn watercourse link
a gauge (or a blend of two gauges) so whole rivers can be coloured.

Inputs
  data/rivers/Data/oprvrs_gb.gpkg   OS Open Rivers (watercourse_link layer)
  data/ne/uk_units.gpkg             Natural Earth map units (England polygon, BNG)
  data/sepa/sepa_snapshot.json      output of 01_fetch_sepa
Output
  build/network.json                grouped, simplified river lines + gauge snapping
"""
import json, sys, os, collections
import numpy as np, geopandas as gpd, shapely
from shapely.ops import linemerge
from pyproj import Transformer

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
os.makedirs(f'{ROOT}/build', exist_ok=True)
X0, X1, Y0, Y1 = 40000, 480000, 525000, 1225000   # map extent, BNG metres
KEEP_UPSTREAM_M = 7000      # draw links draining at least this much network length
SNAP_M = 350                # max gauge-to-river snapping distance
SIMPLIFY_M = 90
Q = 50                      # coordinate quantum, metres

g = gpd.read_file(f'{ROOT}/data/rivers/Data/oprvrs_gb.gpkg', layer='watercourse_link',
                  bbox=(X0, Y0, X1, Y1))
eng = gpd.read_file(f'{ROOT}/data/ne/uk_units.gpkg')
eng = eng[eng.GEOUNIT == 'England'].geometry.iloc[0]
cent = g.geometry.representative_point()
g = g[~cent.within(eng)].reset_index(drop=True)
n = len(g)
print('links in Scotland:', n, 'km:', round(g.length.sum() / 1000))

start = g.start_node.values; end = g.end_node.values; length = g.length.values
by_start = collections.defaultdict(list); by_end = collections.defaultdict(list)
for i in range(n):
    by_start[start[i]].append(i); by_end[end[i]].append(i)

# --- upstream network length (Kahn topological pass, source -> sea) -----------
indeg = np.array([len(by_end[start[i]]) for i in range(n)])
up_len = length.astype(float).copy()
pending = indeg.copy()
queue = [i for i in range(n) if pending[i] == 0]
order = []
while queue:
    i = queue.pop()
    order.append(i)
    outs = by_start[end[i]]
    for j in outs:
        up_len[j] += up_len[i] / len(outs)       # split at braids / distributaries
        pending[j] -= 1
        if pending[j] == 0:
            queue.append(j)
print('ordered', len(order), 'of', n, '(remainder sit on loops)')
done = np.zeros(n, bool); done[order] = True

# main downstream link and main upstream link for each link
down = np.full(n, -1); upm = np.full(n, -1)
for i in range(n):
    outs = by_start[end[i]]
    if outs:
        down[i] = max(outs, key=lambda j: up_len[j])
    ins = by_end[start[i]]
    if ins:
        upm[i] = max(ins, key=lambda j: up_len[j])

# --- gauges -------------------------------------------------------------------
snap = json.load(open(f'{ROOT}/data/sepa/sepa_snapshot.json'))
S = snap['stations']
tr = Transformer.from_crs(4326, 27700, always_xy=True)
sx, sy = tr.transform([s['lon'] for s in S], [s['lat'] for s in S])
tree = shapely.STRtree(g.geometry.values)
st_link = np.full(len(S), -1); st_dist = np.full(len(S), np.nan)
for k in range(len(S)):
    p = shapely.Point(sx[k], sy[k])
    cand = tree.query(p.buffer(SNAP_M))
    if len(cand) == 0:
        continue
    d = np.array([g.geometry.values[c].distance(p) for c in cand])
    ok = d <= SNAP_M
    if not ok.any():
        continue
    cand, d = cand[ok], d[ok]
    # a gauge sits on the river it is named for: among close candidates prefer the larger channel
    near = d <= max(120, d.min() * 1.5)
    best = cand[near][np.argmax(up_len[cand[near]])]
    st_link[k] = best; st_dist[k] = g.geometry.values[best].distance(p)
print('gauges snapped:', (st_link >= 0).sum(), 'of', len(S),
      'median snap distance m:', round(float(np.nanmedian(st_dist)), 1))

gauge_of = {}
for k, l in enumerate(st_link):
    if l >= 0 and (l not in gauge_of or st_dist[k] < st_dist[gauge_of[l]]):
        gauge_of[l] = k

# --- which links to draw ------------------------------------------------------
keep = up_len >= KEEP_UPSTREAM_M
for l in gauge_of:                       # every gauged burn, and its path to the sea
    j = l; hops = 0
    while j >= 0 and not keep[j] and hops < 5000:
        keep[j] = True; j = down[j]; hops += 1
print('links drawn:', keep.sum(), 'km:', round(length[keep].sum() / 1000))

# --- nearest gauge downstream / upstream along the main stem ------------------
NONE = (-1, np.inf)
dn = [NONE] * n; upg = [NONE] * n
for i in reversed(order):                # sea -> source
    if i in gauge_of:
        dn[i] = (gauge_of[i], 0.0)
    elif down[i] >= 0 and done[down[i]]:
        s, d = dn[down[i]]
        if s >= 0:
            dn[i] = (s, d + length[i])
for i in order:                          # source -> sea
    if i in gauge_of:
        upg[i] = (gauge_of[i], 0.0)
    elif upm[i] >= 0:
        s, d = upg[upm[i]]
        if s >= 0:
            upg[i] = (s, d + length[i])

a = np.full(n, -1); b = np.full(n, -1); wa = np.ones(n)
for i in range(n):
    (sd, dd), (su, du) = dn[i], upg[i]
    if sd >= 0 and su >= 0 and sd != su:
        a[i], b[i] = sd, su
        wa[i] = du / (du + dd) if du + dd > 0 else 0.5     # closer gauge weighs more
    elif sd >= 0:
        a[i] = sd
    elif su >= 0:
        a[i] = su
wq = np.round(wa * 4) / 4
swap = (b >= 0) & (wq == 0)
a[swap] = b[swap]; b[swap] = -1; wq[swap] = 1
b[wq == 1] = -1
cls = np.digitize(up_len, [20000, 80000, 300000])           # 0..3 stroke-width class
print('drawn links by state: gauged-blend', int(((b >= 0) & keep).sum()),
      'single gauge', int(((a >= 0) & (b < 0) & keep).sum()),
      'ungauged', int(((a < 0) & keep).sum()))

# --- group, merge, simplify, quantise -----------------------------------------
groups = collections.defaultdict(list)
for i in np.nonzero(keep)[0]:
    groups[(int(a[i]), int(b[i]), float(wq[i]), int(cls[i]))].append(g.geometry.values[i])
out = []; npts = 0
for (ga, gb, w, c), geoms in sorted(groups.items()):
    m = linemerge(shapely.MultiLineString(geoms)) if len(geoms) > 1 else geoms[0]
    m = m.simplify(SIMPLIFY_M, preserve_topology=False)
    parts = list(m.geoms) if m.geom_type == 'MultiLineString' else [m]
    lines = []
    for p in parts:
        xy = np.asarray(p.coords)[:, :2]
        qx = np.round((xy[:, 0] - X0) / Q).astype(int); qy = np.round((Y1 - xy[:, 1]) / Q).astype(int)
        keepm = np.r_[True, (np.diff(qx) != 0) | (np.diff(qy) != 0)]
        qx, qy = qx[keepm], qy[keepm]
        if len(qx) < 2:
            continue
        flat = np.empty(len(qx) * 2, int)
        flat[0], flat[1] = qx[0], qy[0]
        flat[2::2] = np.diff(qx); flat[3::2] = np.diff(qy)       # delta-encoded
        lines.append(flat.tolist()); npts += len(qx)
    if lines:
        out.append({'a': ga, 'b': gb, 'w': w, 'c': c, 'l': lines})
print('groups', len(out), 'points', npts)

# --- label anchors for the principal rivers -----------------------------------
labels = []
wn = g.watercourse_name.fillna('').values
for nm, cond in (('Tay', None), ('Spey', None), ('Tweed', None), ('Clyde', None), ('Forth', None),
                 ('Dee', lambda y: y > 700000), ('Don', lambda y: y > 700000), ('Nith', None),
                 ('Findhorn', None), ('Deveron', None), ('Annan', None), ('Earn', None)):
    idx = [i for i in np.nonzero(wn == 'River ' + nm)[0]
           if cond is None or cond(g.geometry.values[i].centroid.y)]
    if not idx:
        continue
    idx.sort(key=lambda i: up_len[i])
    m = g.geometry.values[idx[int(len(idx) * 0.55)]].interpolate(0.5, normalized=True)
    labels.append({'n': nm, 'x': int(round((m.x - X0) / Q)), 'y': int(round((Y1 - m.y) / Q))})
print('labels', [l['n'] for l in labels])

names = g.watercourse_name.fillna('').values
stations = []
for k, s in enumerate(S):
    l = int(st_link[k])
    stations.append({'x': int(round((sx[k] - X0) / Q)), 'y': int(round((Y1 - sy[k]) / Q)),
                     'river': names[l] if l >= 0 else '', 'snapped': l >= 0,
                     'up_km': round(float(up_len[l]) / 1000, 1) if l >= 0 else None})
json.dump({'extent': [X0, Y0, X1, Y1], 'q': Q, 'groups': out, 'stations': stations, 'labels': labels},
          open(f'{ROOT}/build/network.json', 'w'), separators=(',', ':'))
