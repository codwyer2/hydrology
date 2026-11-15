"""Mosaic OS Terrain 50 ASCII grid tiles (50 m) into a 100 m grid over the map extent.

Usage: python pipeline/02a_dem_mosaic.py [ROOT] [TILE_DIR]
  TILE_DIR  folder of OS Terrain 50 .asc tiles (default ROOT/data/terr50)
Output: ROOT/work/dem100.npy
"""
import numpy as np, glob, os, sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
TILES = sys.argv[2] if len(sys.argv) > 2 else f'{ROOT}/data/terr50'
os.makedirs(f'{ROOT}/work', exist_ok=True)
X0,X1,Y0,Y1=40000,480000,525000,1225000
R=100
W=(X1-X0)//R; H=(Y1-Y0)//R
dem=np.full((H,W),np.nan,np.float32)
n=0
for f in glob.glob(f'{TILES}/**/*.asc', recursive=True):
    with open(f) as fh:
        hdr=[fh.readline().split() for _ in range(5)]
        a=np.loadtxt(fh,dtype=np.float32)
    h={k.lower():float(v) for k,v in hdr}
    if 'xllcorner' not in h or a.shape!=(200,200): print('skip',f,a.shape); continue
    x=int(round(h['xllcorner']/10000)*10000); y=int(round(h['yllcorner']/10000)*10000)
    if x+10000<=X0 or x>=X1 or y+10000<=Y0 or y>=Y1: continue
    assert a.shape==(200,200) and h['cellsize']==50, (f,a.shape,h)
    a=a.reshape(100,2,100,2).mean(axis=(1,3))
    c0=(x-X0)//R; r0=(Y1-(y+10000))//R
    rs=slice(max(r0,0),min(r0+100,H)); cs=slice(max(c0,0),min(c0+100,W))
    dem[rs,cs]=a[rs.start-r0:rs.stop-r0, cs.start-c0:cs.stop-c0]
    n+=1
print(n,'tiles',np.isnan(dem).mean(),np.nanmin(dem),np.nanmax(dem))
np.save(f'{ROOT}/work/dem100.npy',dem)
