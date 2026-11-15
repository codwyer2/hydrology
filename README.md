# Scotland's Hydrological Pulse

Animated map of Scottish river levels. SEPA gauge telemetry is carried along the OS Open Rivers
network over a watercolour relief of OS Terrain 50, with an optional generative soundtrack driven
by the same levels.

`web/` is the published page and is self-contained: open `web/index.html` through any static
server (it loads `data.json` and `basemap.jpg` beside it).

## Page

- **Records:** last 30 days at hourly resolution, or five years at daily resolution.
- **Normalisation:** each gauge against its own five-year record (default), or the low-to-high of the current view.
- **Time:** the west-to-east strip and the selected gauge's trace share one date axis; click or drag either to seek.
- **Rendering:** river and gauge colours ease towards the data (fast rise, slow recession). In fast playback each frame averages the frames it skips. Detail steps down automatically if frames run long.
- **Sound:** Web Audio, no libraries. A filtered saw pad, a pink-noise river bed and an FM-bell arpeggio on D minor pentatonic. Intensity is the mean gauge level from 0.20 to 0.70. The arpeggio's eight steps are eight west-to-east groups of gauges, panned left to right.

## Pipeline

Run from the repository root; each script takes the working root as an optional first argument
(default `.`). Inputs go in `data/`, intermediates in `work/` and `build/` (all git-ignored).

| Script | Does | Needs |
|---|---|---|
| `pipeline/01_fetch_sepa.py` | SEPA KiWIS: 30 days of 15-min stage, five years of daily mean and max per live gauge -> `data/sepa/sepa_snapshot.json` | network access to `timeseries.sepa.org.uk` |
| `pipeline/02a_dem_mosaic.py` | Mosaics OS Terrain 50 ASCII tiles to a 100 m grid -> `work/dem100.npy` | OS Terrain 50 `.asc` tiles in `data/terr50/` (or pass a folder) |
| `pipeline/02_terrain.py` | Multi-directional hillshade, height tint, sea wash, paper grain -> `build/basemap.jpg` | `data/ne/mu.geojson` (Natural Earth 10m admin-0 map units) |
| `pipeline/03_network.py` | Snaps gauges to OS Open Rivers; assigns each drawn link a gauge or a distance blend of two -> `build/network.json` | `data/rivers/Data/oprvrs_gb.gpkg`, `data/ne/uk_units.gpkg` |
| `pipeline/04_pack.py` | Normalisation knots, 30-day and five-year series, packing -> `build/data.json` | outputs above |

Then copy `build/basemap.jpg` and `build/data.json` into `web/`.

`pip install -r requirements.txt`

`data/ne/uk_units.gpkg` is the England, Wales, Scotland and Northern Ireland rows of the Natural
Earth map units, reprojected to EPSG:27700.

The committed snapshot was taken 2026-11-15 19:18 UTC. `01_fetch_sepa.py` is a port of the
browser queries that produced it and has not itself been run end to end.

## Method

- **Normalisation knots** (per gauge, from five years of SEPA daily series): 1st percentile of daily mean = 0.00, median = 0.32, 90th percentile = 0.62, 99th percentile = 0.84, highest daily maximum = 1.00.
- **River lines:** OS Open Rivers links draining at least 7 km of network (about 18,400 km). Gauges snap to the nearest link within 350 m (342 of 348; median 13 m). A link takes the first gauge downstream; between two gauges on one stem it blends them by distance. This is display interpolation, not a hydraulic model.

## Licences

Contains SEPA data © Scottish Environment Protection Agency and database right 2026.
Contains OS data © Crown copyright and database right 2026. Made with Natural Earth.
