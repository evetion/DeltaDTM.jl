import numpy as np
import rasterio as rio
import seaborn as sns
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cmocean.cm as cm
import matplotlib
from rasterio.plot import show
from rasterio.windows import transform
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from os.path import isfile, splitext
import matplotlib.patheffects as pe

sns.set_context("paper")
sns.set_style("ticks")


matplotlib.rcParams["font.size"] = "10"


datadir = "data"

ocean_10m = cfeature.NaturalEarthFeature(
    "physical", "ocean", "10m", facecolor="lightgrey"
)

shp = shapereader.Reader("data/water-polygons-split-4326")
esa = "data/esa-worldcover-2021/copernicus2/"

esahcolors = [
    ("#ffffff", 0),
    ("#006400", 10),
    ("#ffbb22", 20),
    ("#ffff4c", 30),
    ("#f096ff", 40),
    ("#fa0000", 50),
    ("#b4b4b4", 60),
    ("#f0f0f0", 70),
    ("#0064c8", 80),
    ("#0096a0", 90),
    ("#00cf75", 95),
    ("#fae6a0", 100),
]
esacolors = [(v / 100, colors.to_rgb(c)) for (c, v) in esahcolors]
esacmap = colors.LinearSegmentedColormap.from_list("esa", esacolors)

hlimit = 10
llimit = 0
scalebar_size = 0.1
compare = {
    # "NASADEM": "data/nasadem",
    # "CopernicusDEM": "data/copernicus/copernicus2",
    "MERIT": "data/merit",
    "CoastalDEM": "data/coastaldem",
    "FABDEM": "data/FABDEM",
    "DeltaDEM": "code/DeltaDEM/data/deltadem/v1",
}


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def scalebar_size(latitude, factor=0.1, width_of_degree=111_320):
    w = abs(np.cos(np.deg2rad(latitude))) * width_of_degree
    exp = np.log10(w * factor)
    exp = np.ceil(exp)
    s = 10**exp / w * factor
    ss = round(s * abs(np.cos(np.deg2rad(latitude))) * width_of_degree)
    return s, ss


diffs = {}

def plot_compare(tile, window, x, y, vmin=0, vmax=10, vs=10):
    fig = plt.figure(figsize=(8.27, 3.31), dpi=300)

    spec = gridspec.GridSpec(
        ncols=len(compare) + 2,
        nrows=2,
        figure=fig,
        wspace=0.05,
        hspace=0.05,
        width_ratios=(100,) * (len(compare) + 1) + (5,),
    )

    # Always make it square!
    window = rio.windows.Window(*w)

    # deep = plt.get_cmap(cm.deep, 10)
    deep = discrete_cmap(vs, cm.deep_r)
    curl = discrete_cmap(20, cm.curl)

    src = rio.open(f"{datadir}/validation/{tile}")
    val = src.read(1, masked=True, window=window)
    val = val.filled(np.nan)
    mask = val > hlimit
    val[mask] = np.nan
    m = np.isnan(val)
    ax = fig.add_subplot(spec[0, 0], projection=ccrs.PlateCarree())

    show(
        val,
        ax=ax,
        cmap=deep,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    geometries = []
    for row in shp.records():
        if row.attributes["x"] == x and row.attributes["y"] == y:
            geometries.append(row.geometry)
    ax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    ax.set_title("Reference")
    ax.set_ylabel("Latitude [°]")
    gl = ax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    ax.set_aspect("auto")
    for i, (key, fn) in enumerate(compare.items()):
        ii = i + 1
        fn = f"{fn}/{tile}"
        if not isfile(fn):
            print(f"Skipping {key} because {fn} does not exist.")
            continue
        src = rio.open(fn)
        cd = src.read(1, masked=True, window=window)
        cd = cd.astype("float32")
        cd = cd.filled(np.nan)
        cd[m] = np.nan

        hax = fig.add_subplot(spec[0, ii], projection=ccrs.PlateCarree())
        hax.add_geometries(
            geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1
        )

        imh = show(
            cd,
            ax=hax,
            cmap=deep,
            vmin=vmin,
            vmax=vmax,
            interpolation="none",
            transform=transform(window, src.transform),
        )
        hax.tick_params(direction="in")
        gl = hax.gridlines(
            xlocs=[x + 0.33, x + 0.66],
            ylocs=[y + 0.33, y + 0.66],
            rotate_labels=False,
            draw_labels=True,
            dms=False,
            x_inline=False,
            y_inline=False,
            linewidth=0.5,
        )
        gl.top_labels = False
        gl.left_labels = False
        gl.bottom_labels = False
        gl.right_labels = False
        hax.set_title(key)
        diff = cd - val

        if key in diffs:
            diffs[key] = np.concatenate((diffs[key], cd[~np.isnan(diff)].ravel()))
        else:
            diffs[key] = cd[~np.isnan(diff)].ravel()

        if key == "FABDEM":
            k = "Reference"
            if k in diffs:
                diffs[k] = np.concatenate((diffs[k], val[~np.isnan(diff)].ravel()))
            else:
                diffs[k] = val[~np.isnan(diff)].ravel()

        hax.set_xticklabels([])
        hax.set_yticklabels([])
        hax.set_aspect("auto")

        dax = fig.add_subplot(spec[1, ii], projection=ccrs.PlateCarree())
        dax.add_geometries(
            geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1
        )
        imd = show(
            diff,
            ax=dax,
            cmap=curl,
            vmin=-5,
            vmax=5,
            transform=transform(window, src.transform),
        )
        dax.tick_params(direction="in")
        dax.set_xlabel("Longitude [°]")
        dax.set_yticklabels([])
        dax.set_aspect("auto")

        gl = dax.gridlines(
            xlocs=[x + 0.33, x + 0.66],
            ylocs=[y + 0.33, y + 0.66],
            rotate_labels=False,
            draw_labels=True,
            dms=False,
            x_inline=False,
            y_inline=False,
            linewidth=0.5,
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.left_labels = False
        gl.bottom_labels = True

        if i == len(compare) - 1:
            cax = fig.add_subplot(spec[0, ii + 1])
            plt.colorbar(
                imh.get_images()[0],
                cax=cax,
                extend="both",
                label="Elevation [m +EGM08]",
            )
            cax = fig.add_subplot(spec[1, ii + 1])
            plt.colorbar(
                imd.get_images()[0],
                cax=cax,
                extend="both",
                label="Difference [m]",
            )

    _, lat = src.xy(0, 0)
    f = 0.1
    s, ss = scalebar_size(lat, f)
    scalebar = AnchoredSizeBar(
        ax.transData,
        s,
        f"{int(ss/1000)} km",
        "lower center",
        color="black",
        # sep=2,
        pad=0.15,
        # borderpad=2,
        size_vertical=0.001,
        frameon=True,
    )
    scalebar.set(path_effects=[pe.withStroke(linewidth=1.5, foreground="red")])
    scalebar.set_path_effects([pe.withStroke(linewidth=1.5, foreground="red")])
    scalebar.get_child()._children[0].set_path_effects(
        [pe.withStroke(linewidth=1.5, foreground="red")]
    )

    ax.add_artist(scalebar)

    # ESA
    gax = fig.add_subplot(spec[1, 0], projection=ccrs.PlateCarree())
    gax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    fn = f"{esa}/{tile}"
    src = rio.open(fn)
    cd = src.read(1, masked=True, window=window).astype("float")
    cd[m] = np.nan
    imh = show(
        cd,
        ax=gax,
        cmap=esacmap,
        vmin=0,
        vmax=100,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    gax.set_aspect("auto")

    gl = gax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.right_labels = False

    fig.savefig(f"v1/{tile}_compare_all.png", bbox_inches="tight", dpi=600)
    plt.close(fig)
    return None

orig_folder = "data/copernicus/copernicus_orig/elevation/"
dd_folder = compare["DeltaDEM"]

def plot_explain(tile, window, x, y, vmin=0, vmax=10, vs=10):
    fig = plt.figure(figsize=(3.31, 5.1), dpi=300)

    spec = gridspec.GridSpec(
        ncols=3,
        nrows=3,
        figure=fig,
        wspace=0.05,
        hspace=0.25,
        width_ratios=(100,) * (2) + (5,),
        height_ratios=(100,) * (3),
    )

    # Always make it square!
    window = rio.windows.Window(*w)

    deep = discrete_cmap(vs, cm.deep_r)
    tempo = discrete_cmap(vs, cm.tempo)

    src = rio.open(f"{datadir}/validation/{tile}")
    val = src.read(1, masked=True, window=window)
    val = val.filled(np.nan)
    mask = val > hlimit
    val[mask] = np.nan
    m = np.isnan(val)
    ax = fig.add_subplot(spec[0, 1], projection=ccrs.PlateCarree())

    show(
        val,
        ax=ax,
        cmap=deep,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    geometries = []
    for row in shp.records():
        if row.attributes["x"] == x and row.attributes["y"] == y:
            geometries.append(row.geometry)
    ax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    ax.set_title("Reference")
    ax.set_ylabel("Latitude [°]")
    gl = ax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        ylabel_style={"rotation": 90},  # xlabel_style={"fontsize": 5},
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    ax.set_aspect("auto")

    # CopernicusDEM
    fn = f"{orig_folder}/{tile}"
    src = rio.open(fn)
    ori = src.read(1, masked=True, window=window)
    ori = ori.astype("float32")
    ori = ori.filled(np.nan)
    ori[m] = np.nan

    hax = fig.add_subplot(spec[0, 0], projection=ccrs.PlateCarree())
    hax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    imh = show(
        ori,
        ax=hax,
        cmap=deep,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    hax.tick_params(direction="in")
    gl = hax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        ylabel_style={"rotation": 90},  # xlabel_style={"fontsize": 5},
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = False
    gl.right_labels = False
    hax.set_title("CopernicusDEM")
    hax.set_xticklabels([])
    hax.set_yticklabels([])
    hax.set_aspect("auto")

    # Terrain
    tilen = splitext(tile)[0] + "_f.tif"
    fn = f"{dd_folder}/{tilen}"
    src = rio.open(fn)
    cd = src.read(1, masked=True, window=window)
    cd = cd.astype("float32")
    cd = cd.filled(np.nan)
    cd[m] = np.nan

    hax = fig.add_subplot(spec[1, 0], projection=ccrs.PlateCarree())
    hax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    imh = show(
        cd,
        ax=hax,
        cmap=deep,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    hax.set_title("Terrain")
    hax.set_aspect("auto")

    gl = hax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        ylabel_style={"rotation": 90},  # xlabel_style={"fontsize": 5},
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = False
    gl.right_labels = False

    # DeltaDEM
    fn = f"{dd_folder}/{tile}"
    src = rio.open(fn)
    dd = src.read(1, masked=True, window=window)
    dd = dd.astype("float32")
    dd = dd.filled(np.nan)
    dd[m] = np.nan

    hax = fig.add_subplot(spec[2, 0], projection=ccrs.PlateCarree())
    hax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    imh = show(
        dd,
        ax=hax,
        cmap=deep,
        vmin=vmin,
        vmax=vmax,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    hax.set_title("DeltaDEM")
    hax.set_aspect("auto")

    gl = hax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        ylabel_style={"rotation": 90},
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = True
    gl.bottom_labels = True
    gl.right_labels = False

    # nDSM (by DeltaDEM)
    hax = fig.add_subplot(spec[2, 1], projection=ccrs.PlateCarree())
    hax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    imndsm = show(
        ori - dd,
        ax=hax,
        cmap=tempo,
        vmin=0,
        vmax=25,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    hax.set_title("normalised DSM")
    hax.set_aspect("auto")

    gl = hax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        ylabel_style={"rotation": 90},
      
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = False
    gl.bottom_labels = True
    gl.right_labels = False

    # Colorbar
    cax = fig.add_subplot(spec[0, 2])
    plt.colorbar(
        imh.get_images()[0], cax=cax, extend="both", label="Elevation [m +EGM08]"
    )
    cax = fig.add_subplot(spec[2, 2])
    plt.colorbar(
        imndsm.get_images()[0], cax=cax, extend="both", label="Elevation [m +terrain]"
    )

    _, lat = src.xy(0, 0)
    f = 0.1
    s, ss = scalebar_size(lat, f)
    scalebar = AnchoredSizeBar(
        ax.transData,
        s,
        f"{int(ss/1000)} km",
        "lower center",
        color="black",
        # sep=2,
        pad=0.15,
        # borderpad=2,
        size_vertical=0.001,
        frameon=True,
    )
    scalebar.set(path_effects=[pe.withStroke(linewidth=1.5, foreground="red")])
    scalebar.set_path_effects([pe.withStroke(linewidth=1.5, foreground="red")])
    scalebar.get_child()._children[0].set_path_effects(
        [pe.withStroke(linewidth=1.5, foreground="red")]
    )

    ax.add_artist(scalebar)

    # ESA
    gax = fig.add_subplot(spec[1, 1], projection=ccrs.PlateCarree())
    gax.add_geometries(geometries, ccrs.PlateCarree(), facecolor="lightgray", zorder=-1)
    fn = f"{esa}/{tile}"
    src = rio.open(fn)
    cd = src.read(1, masked=True, window=window).astype("float")
    cd[m] = np.nan
    imh = show(
        cd,
        ax=gax,
        cmap=esacmap,
        vmin=0,
        vmax=100,
        interpolation="none",
        transform=transform(window, src.transform),
    )
    gax.set_aspect("auto")
    gax.set_title("ESA WorldCover")

    gl = gax.gridlines(
        xlocs=[x + 0.33, x + 0.66],
        ylocs=[y + 0.33, y + 0.66],
        # ylabel_style={"fontsize": 5},
        # xlabel_style={"fontsize": 5},
        rotate_labels=False,
        draw_labels=True,
        dms=False,
        x_inline=False,
        y_inline=False,
        linewidth=0.5,
    )
    gl.top_labels = False
    gl.left_labels = False
    gl.bottom_labels = False
    gl.right_labels = False


    fig.savefig(f"v1/{tile}_explain.png", bbox_inches="tight", dpi=600)
    plt.close(fig)
    return None

tiles = [
    (
        "Copernicus_DSM_10_N25_00_W082_00_DEM.tif",
        25,
        -82,
        (500, 0, 3100, 3100),
        (-2, 8, 10),
    ),  # florida
    (
        "Copernicus_DSM_10_N25_00_W081_00_DEM.tif",
        25,
        -81,
        (1600, 0, 2000, 2000),
        (0, 10, 10),
    ),  # miami south
    (
        "Copernicus_DSM_10_N26_00_W082_00_DEM.tif",
        26,
        -82,
        (0, 500, 2000, 2000),
        (0, 10, 10),
    ),  # florida north
    (
        "Copernicus_DSM_10_N07_00_E171_00_DEM.tif",
        7,
        171,
        (0, 2200, 1400, 1400),
        (0, 10, 10),
    ),  # majuro
    (
        "Copernicus_DSM_10_N52_00_E004_00_DEM.tif",
        52,
        4,
        (800, 1200, 1600, 1600),
        (-5, 10, 15),
    ),  # nl
    (
        "Copernicus_DSM_10_N52_00_E005_00_DEM.tif",
        52,
        5,
        (0, 0, 3600, 3600),
        (-5, 10, 15),
    ),  # nl mid
    (
        "Copernicus_DSM_10_N51_00_E003_00_DEM.tif",
        51,
        3,
        (1000, 400, 1400, 1400),
        (-5, 10, 15),
    ),  # zeeland
    (
        "Copernicus_DSM_10_N53_00_E005_00_DEM.tif",
        53,
        5,
        (600, 1800, 1800, 1800),
        (-5, 10, 15),
    ),  # friesland
    (
        "Copernicus_DSM_10_S03_00_E104_00_DEM.tif",
        -3,
        104,
        (1800, 0, 1800, 1800),
        (0, 10, 10),
    ),  # sumatra
    (
        "Copernicus_DSM_10_S02_00_E104_00_DEM.tif",
        -2,
        104,
        (0, 0, 1800, 1800),
        (0, 10, 10),
    ),  # sumatra 2
    (
        "Copernicus_DSM_10_S13_00_E096_00_DEM.tif",
        -13,
        96,
        (2700, 100, 800, 800),
        (0, 10, 10),
    ),  # cocos
    (
        "Copernicus_DSM_10_S13_00_E130_00_DEM.tif",
        -13,
        130,
        (1800, 1000, 1800, 1800),
        (0, 10, 10),
    ),  # australia
    (
        "Copernicus_DSM_10_S13_00_E132_00_DEM.tif",
        -13,
        132,
        (0, 0, 3600, 3600),
        (0, 10, 10),
    ),  # australia 2
    (
        "Copernicus_DSM_10_N18_00_W094_00_DEM.tif",
        18,
        -94,
        (1800, 1800, 1800, 1800),
        (0, 10, 10),
    ),  # mexico
    (
        "Copernicus_DSM_10_N18_00_W093_00_DEM.tif",
        18,
        -93,
        (1800, 500, 2000, 2000),
        (-2, 8, 10),
    ),  # mexico 2
    (
        "Copernicus_DSM_10_N52_00_E000_00_DEM.tif",
        52,
        0,
        (0, 0, 1800, 1800),
        (-2, 10, 12),
    ),  # uk
    (
        "Copernicus_DSM_10_N52_00_W001_00_DEM.tif",
        52,
        -1,
        (1200, 200, 1800, 1800),
        (-2, 10, 12),
    ),  # uk 2
    (
        "Copernicus_DSM_10_N56_00_E023_00_DEM.tif",
        56,
        23,
        (600, 0, 1800, 1800),
        (0, 10, 10),
    ),  # latvia
    (
        "Copernicus_DSM_10_N54_00_E018_00_DEM.tif",
        54,
        18,
        (800, 2000, 1600, 1600),
        (-2, 10, 12),
    ),  # poland
    (
        "Copernicus_DSM_10_N54_00_E019_00_DEM.tif",
        54,
        19,
        (0, 2000, 1600, 1600),
        (-2, 10, 12),
    ),  # poland 2
    (
        "Copernicus_DSM_10_S03_00_E114_00_DEM.tif",
        -3,
        114,
        (800, 800, 2000, 2000),
        (0, 10, 10),
    ),  # kalimantan 1
    (
        "Copernicus_DSM_10_S04_00_E114_00_DEM.tif",
        -4,
        114,
        (300, 0, 1000, 1000),
        (0, 10, 10),
    ),  # kalimantan 2
]
for tile, y, x, w, (vmin, vmax, vs) in tiles:
    plot_compare(tile, w, x, y, vmin=vmin, vmax=vmax, vs=vs)
    plot_explain(tile, w, x, y, vmin=vmin, vmax=vmax, vs=vs)


# Locations
fig = plt.figure(figsize=(6, 3), dpi=600)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.add_feature(cfeature.OCEAN, facecolor="lightgrey", zorder=-1)
for tile, y, x, w, (vmin, vmax, vs) in tiles:
    rect = patches.Rectangle(
        (x, y),
        1,
        1,
        linewidth=0.5,
        edgecolor="r",
        facecolor="none",
        transform=ccrs.Geodetic(),
    )
    ax.add_patch(rect)

ax.set_global()
ax.spines["geo"].set_visible(False)
plt.savefig("locations.pdf", dpi=600, bbox_inches="tight")


pltk = {
    "Reference": {"label": "Reference", "c": "black", "lw": 3, "zorder": -10.0},
    "MERIT": {"label": "MERIT", "c": "purple", "lw": 1, "linestyle": "dashed"},
    "CoastalDEM": {"label": "CoastalDEM", "c": "blue", "lw": 1, "linestyle": "dashdot"},
    "FABDEM": {"label": "FABDEM", "c": "green", "lw": 1, "linestyle": "dotted"},
    "DeltaDEM": {"label": "DeltaDEM", "c": "red", "lw": 1, "zorder": -5},
}


sns.set_context("paper")
sns.set_style("ticks")

fig = plt.figure(figsize=(6, 4), dpi=600)
ax = fig.add_subplot(1, 1, 1)

for k in pltk.keys():
    v = diffs[k]
    v = v[~np.isnan(v)]
    x = np.flip(np.sort(v))
    y = 100 * np.arange(len(v)) / (len(v) - 1)
    ax.plot(y, x, **pltk[k])

ax.set_ylim(-6, 12)
ax.set_xlabel("Cumulative percentage [%]")
ax.set_ylabel("Height [m +EGM08]")
ax.legend()
sns.despine(trim=True)
plt.savefig("cdf.pdf", dpi=600)
plt.close(fig)
