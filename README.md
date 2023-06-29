# DeltaDEM.jl
Julia package used to create DeltaDEM

In `src` there's the source code of the package. More useful are the files in `scripts`, which make use of the `src` files.
They can be called like so:

- `julia --project scripts/bias_call` to process the bias correction
- `julia --project scripts/dtm_call` to generate DeltaDEM tiles
- `julia --project scripts/make_release` to generate DeltaDEM zips

## Prerequisites
To run this as is, one needs at least CopernicusDEM tiles (used both with and without water masked out), the original CopernicusDEM error and mask files. One also need to download ESA-WorldCover 2021 and subtile the dataset into CopernicusDEM tiles. Last but not least, one needs to download all ICESat-2 and GEDI data for their area of interest (can be done with `SpaceLiDAR.jl`) and save these as GeoParquet (`.pq`), also tiled to the CopernicusDEM tiles, by saving multiple GeoParquet files (one for each granule) in a folder per tile.

For all these `.tif` tiles, `gdalbuildvrt` is run, to include all direct neighbours of such tiles. Each 1 by 1 degree `.tif` tile thus as a `.tif.vrt` tile that's 3 by 3 degrees. For the GeoParquet files `.pq` a similar buffer is made, by making a text (`.txt`) file containing the paths of all files, including all filepaths in the neighboring tilefolders.

## Plots
The plots are made in Python (with the exception of the bias plot in `src/bias.jl`), and can be found in the `plots` folder.

## Tables
The validation tables can be generated as `.csv` by calling `scripts/validation.jl`, once DeltaDEM has been generated for areas with a validation tile, also tiled to the CopernicusDEM specification.
