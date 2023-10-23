import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import cmocean.cm as cm
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd

sns.set_context("paper")
sns.set_style("ticks")

df = gpd.read_file(
    "DeltaDTM/crossval_tiles.gpkg"
)
df.longitude = df.geometry.x
df.latitude = df.geometry.x


def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


cmap = discrete_cmap(20, cm.curl)
fig = plt.figure(figsize=(8.27, 4.02), dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
cb = plt.scatter(
    df.geometry.x,
    df.geometry.y,
    transform=ccrs.Geodetic(),
    s=(np.log10(df.n) ** 2) / 10,
    c=df["mean"],
    cmap=cmap,
    vmin=-5,
    vmax=5,
    edgecolors="none",
    linewidths=0,
    zorder=10,
)
ax.spines["geo"].set_visible(False)
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="2%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)
c = plt.colorbar(cb, cax=ax_cb, extend="both", label="Mean error [m]")
plt.savefig("crossval_globe_mean.png", dpi=300, bbox_inches="tight")


df = gpd.read_file(
    "data/biasv1.gpkg"
)
df.longitude = df.geometry.x
df.latitude = df.geometry.x

cmap = discrete_cmap(8, cm.curl)
fig = plt.figure(figsize=(6, 3), dpi=600)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
cb = plt.scatter(
    df.geometry.x,
    df.geometry.y,
    transform=ccrs.Geodetic(),
    s=(np.log10(df.n) ** 2) / 250,
    c=df["bias"],
    cmap=cmap,
    vmin=-1,
    vmax=1,
    edgecolors="none",
    linewidths=0,
    zorder=10,
)
ax.spines["geo"].set_visible(False)
divider = make_axes_locatable(ax)
ax_cb = divider.new_horizontal(size="2%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)
c = plt.colorbar(cb, cax=ax_cb, extend="both", label="Bias [m]")
plt.savefig("bias_globe.png", dpi=600, bbox_inches="tight")
