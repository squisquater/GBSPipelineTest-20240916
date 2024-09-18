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

# Fit FEEMS
lambda_val = float(snakemake.params.lambda_val)
sp_graph.fit(lamb=lambda_val)

# Create a plot of the fitted graph
projection = ccrs.AlbersEqualArea(central_longitude=-119.624, central_latitude=35.602,
                                  standard_parallels=(34, 42))
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=5, cbar_ticklabelsize=5)

# Adjust colorbar properties
v.cbar_width = "100%"
v.cbar_height = "10%"
v.cbar_bbox_to_anchor = (0.03, 0.05, 0.2, 0.3)

v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()

# Save the plot
fig.savefig(snakemake.output.feems_img)
