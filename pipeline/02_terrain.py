"""Render the watercolour base map: multi-directional hillshade of OS Terrain 50
tinted in parchment / ochre, cerulean sea with a shoreline wash, paper grain.

Inputs
  work/dem100.npy            100 m mosaic of OS Terrain 50 (see 02a_dem_mosaic.py)
  data/ne/mu.geojson         Natural Earth 10m map units (England / Ireland outlines)
Output
  build/basemap.jpg          1320 x 2100 px, 333.3 m per pixel
"""
import sys, os, numpy as np, geopandas as gpd
from scipy import ndimage as ndi
from PIL import Image, ImageDraw

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
X0, X1, Y0, Y1 = 40000, 480000, 525000, 1225000
R = 100; OUT_W, OUT_H = 1320, 2100
rng = np.random.default_rng(7)
os.makedirs(f'{ROOT}/build', exist_ok=True)

dem = np.load(f'{ROOT}/work/dem100.npy')
H, W = dem.shape
nodata = np.isnan(dem)
z = np.where(nodata, 0, dem)

# --- water masks: Terrain 50 flattens the sea (to low-water level) and lochs --
def eq(a, axis):                         # equals a neighbour along an axis (within a tile)
    d = np.abs(np.diff(a, axis=axis)) < 0.011
    pad = [(0, 0), (0, 0)]; pad[axis] = (1, 1)
    d = np.pad(d, pad)
    return (d[1:, :] | d[:-1, :]) if axis == 0 else (d[:, 1:] | d[:, :-1])
flat = eq(z, 0) & eq(z, 1)
sea = nodata | (flat & (z < 0.6))
sea = ndi.binary_closing(sea, iterations=1) | sea        # heal tile seams
lab, nl = ndi.label(sea)                      # keep only sea connected to open water
sizes = ndi.sum(sea, lab, np.arange(1, nl + 1))
sea = np.isin(lab, 1 + np.nonzero(sizes > 400)[0])
loch = flat & ~sea & (z >= 0.6)
lab, nl = ndi.label(loch)
sizes = ndi.sum(loch, lab, np.arange(1, nl + 1))
loch = np.isin(lab, 1 + np.nonzero(sizes > 40)[0])      # > 0.4 km2
land = ~sea
print('land %.3f loch km2 %.0f' % (land.mean(), loch.sum() / 100))

# --- neighbouring land from Natural Earth (flat tint, no relief) --------------
mu = gpd.read_file(f'{ROOT}/data/ne/mu.geojson')
def raster(names):
    geo = mu[mu.GEOUNIT.isin(names)].to_crs(27700).geometry
    im = Image.new('L', (W, H), 0); d = ImageDraw.Draw(im)
    for geom in geo:
        for poly in getattr(geom, 'geoms', [geom]):
            xy = [((x - X0) / R, (Y1 - y) / R) for x, y in poly.exterior.coords]
            d.polygon(xy, fill=255)
    return np.asarray(im) > 0
england = raster(['England'])
ireland = raster(['Northern Ireland', 'Ireland']) & nodata
sea &= ~ireland
scot = land & ~england & ~ireland
neighbour = (land & england) | ireland

# --- multi-directional hillshade ----------------------------------------------
def hillshade(zz, zf):
    gy, gx = np.gradient(zz, R)               # gy: +row = south
    nx, ny, nz = -zf * gx, zf * gy, 1.0       # surface normal (east, north, up)
    norm = np.sqrt(nx * nx + ny * ny + 1.0)
    out = np.zeros_like(zz, dtype=np.float32)
    for az, wt in ((315, .50), (270, .20), (0, .18), (225, .12)):
        a, alt = np.radians(az), np.radians(45)
        l = (np.sin(a) * np.cos(alt), np.cos(a) * np.cos(alt), np.sin(alt))
        out += wt * np.clip((nx * l[0] + ny * l[1] + nz * l[2]) / norm, 0, 1)
    return out
zs = np.where(land, z, 0).astype(np.float32)
hs = 0.6 * hillshade(ndi.gaussian_filter(zs, 2.5), 4.0) + 0.4 * hillshade(ndi.gaussian_filter(zs, 0.8), 3.0)
flat_ref = np.sin(np.radians(45))
shade = np.clip((hs - flat_ref) / flat_ref, -1, 1)          # -1 deep shadow .. +0.4 lit

# --- colour --------------------------------------------------------------------
def hexrgb(h): return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], np.float32) / 255
ramp_e = np.array([0, 60, 250, 550, 900, 1300], np.float32)
ramp_c = np.stack([hexrgb(c) for c in ('#F7F0DC', '#F1E4C3', '#EAD8B1', '#D3B47E', '#A98755', '#8C6D46')])
e = ndi.gaussian_filter(zs, 1.2)
tint = np.stack([np.interp(e, ramp_e, ramp_c[:, k]) for k in range(3)], -1)
umber = hexrgb('#5E4326'); cream = hexrgb('#FFFBEF')
sh = shade[..., None]
rgb = np.where(sh < 0, tint + (umber - tint) * np.clip(-sh * 0.85, 0, .8), tint + (cream - tint) * np.clip(sh * 1.1, 0, .6))

# neighbouring land: same relief where we have it, washed out
pale = hexrgb('#EFE9DA')
rgb[neighbour] = pale + (rgb[neighbour] - pale) * 0.28
rgb[ireland] = pale * 0.985

# sea with shoreline wash
sea_c, sea_deep, sea_edge = hexrgb('#D9E8F2'), hexrgb('#C3D9E8'), hexrgb('#A9C6DB')
dist = ndi.distance_transform_edt(~land & ~ireland) * R / 1000.0          # km from shore
seargb = sea_c + (sea_deep - sea_c) * (np.exp(-dist / 14.0) * 0.9)[..., None]
seargb = seargb + (sea_edge - seargb) * (np.exp(-dist / 1.6) * 0.75)[..., None]
water = sea & ~ireland
rgb[water] = seargb[water]
loch_c = hexrgb('#B4CDE0')
rgb[loch & scot] = loch_c
rgb[loch & neighbour] = pale + (loch_c - pale) * 0.5

# --- resample to output, then ink line + paper texture at output resolution ----
def down(a):
    im = Image.fromarray((np.clip(a, 0, 1) * 255).astype(np.uint8))
    return np.asarray(im.resize((OUT_W, OUT_H), Image.LANCZOS)).astype(np.float32) / 255
img = down(rgb)
landf = down(land.astype(np.float32) + ireland); scotf = down(scot.astype(np.float32))
edge = np.clip(1 - np.abs(landf - 0.5) * 2.4, 0, 1) ** 1.5                # coastline ink
ink = hexrgb('#6B5236')
img = img + (ink - img) * (edge * 0.55)[..., None]

def noise(sigma):
    nz = ndi.gaussian_filter(rng.standard_normal((OUT_H, OUT_W)).astype(np.float32), sigma)
    return nz / nz.std()
wash = 0.55 * noise(60) + 0.3 * noise(18) + 0.15 * noise(6)               # pigment pooling
grain = 0.6 * noise(0.7) + 0.4 * noise(1.6)                               # paper tooth
fib = ndi.gaussian_filter(rng.standard_normal((OUT_H, OUT_W)).astype(np.float32), (0.6, 5)); fib /= fib.std()
tex = 1 + 0.030 * wash + 0.022 * grain + 0.010 * fib
img = np.clip(img * tex[..., None], 0, 1)
Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8)).save(f'{ROOT}/build/basemap.jpg', quality=86, subsampling=0, optimize=True)
Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8)).resize((528, 840), Image.LANCZOS).save(f'{ROOT}/build/basemap_preview.jpg', quality=80)
Image.fromarray((scotf > 0.5).astype(np.uint8) * 255).save(f'{ROOT}/build/scotmask.png')
