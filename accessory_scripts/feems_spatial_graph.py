import numpy as np
import pkg_resources
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Read in the coordinate and outer files
coord = np.loadtxt(snakemake.input.coord)
outer = np.loadtxt(snakemake.input.outer)
grid_path = snakemake.input.grid

# Prepare the graph inputs
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)

# Setup the spatial graph
sp_graph = SpatialGraph(genotypes=None, coord=coord, grid=grid, edges=edges, scale_snps=True)

# Create a plot
projection = ccrs.AlbersEqualArea(central_longitude=-119.624, central_latitude=35.602,
                                  standard_parallels=(34, 42))
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)

# Save the plot
fig.savefig(snakemake.output.graph_img)
