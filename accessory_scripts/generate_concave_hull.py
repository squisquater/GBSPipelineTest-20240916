import pandas as pd
from shapely.geometry import Polygon, LineString
import alphashape
import matplotlib.pyplot as plt
import sys

# Get the input parameters from Snakemake
alpha = float(snakemake.params.alpha)
buffer_distance = float(snakemake.params.buffer_distance)

# Read the 'individual_centroids.csv' file
df = pd.read_csv(snakemake.input.centroids)

# Convert to a list of tuples (longitude, latitude)
points = [(lon, lat) for lon, lat in zip(df['longitude'], df['latitude'])]

# Generate the alpha shape (concave hull)
concave_hull = alphashape.alphashape(points, alpha)

# Buffer the polygon (this will create a buffered polygon around the existing polygon)
buffered_polygon = concave_hull.buffer(buffer_distance)

# Optionally, smooth the polygon by interpolating points on the buffered polygon
def smooth_polygon(polygon, num_points=100):
    # Create a LineString from the polygon's exterior
    line = LineString(polygon.exterior.coords)
    
    # Interpolate points evenly spaced along the length of the LineString
    interpolated_line = [line.interpolate(float(i) / num_points, normalized=True) for i in range(num_points)]
    
    # Create a smoothed polygon from the interpolated points
    smoothed_polygon = Polygon([(p.x, p.y) for p in interpolated_line])
    
    return smoothed_polygon

# Apply the smoothing function to the buffered polygon
smoothed_polygon = smooth_polygon(buffered_polygon, num_points=200)

# Plot the result
fig, ax = plt.subplots()

# Plot the original data points
ax.scatter(df['longitude'], df['latitude'], color='blue', label='Data Points', zorder=5)

# Plot the original concave hull
x_concave, y_concave = concave_hull.exterior.xy
ax.plot(x_concave, y_concave, color='red', label='Concave Hull', zorder=3)

# Plot the smoothed, buffered polygon
x_smoothed, y_smoothed = smoothed_polygon.exterior.xy
ax.plot(x_smoothed, y_smoothed, color='green', label='Smoothed Buffered Polygon', zorder=4)

# Set axis labels
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_title(f"Smoothed Buffered Polygon (alpha={alpha}, buffer={buffer_distance})")

# Save the plot to an external .png file
plt.legend()
plt.savefig(snakemake.output.polygon_img, dpi=300)
plt.close()

# Save the smoothed buffered polygon coordinates to a tab-delimited file without headers
def save_polygon_coordinates(polygon, filename):
    with open(filename, 'w') as f:
        for coord in polygon.exterior.coords:
            f.write(f"{coord[0]}\t{coord[1]}\n")

# Save the coordinates of the smoothed buffered polygon
save_polygon_coordinates(smoothed_polygon, snakemake.output.outer_coords)
