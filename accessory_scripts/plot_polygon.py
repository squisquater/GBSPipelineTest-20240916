import matplotlib.pyplot as plt

def plot_data_and_polygon(df, smoothed_polygon, output_file="smoothed_polygon_plot.png"):
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the original points
    ax.scatter(df['longitude'], df['latitude'], color='blue', label='Data Points', zorder=5)

    # Plot the smoothed polygon
    polygon_x, polygon_y = smoothed_polygon.exterior.xy
    ax.plot(polygon_x, polygon_y, color='red', label='Smoothed Polygon', zorder=3)

    # Add labels to the data points
    for i, row in df.iterrows():
        ax.text(row['longitude'], row['latitude'], row['FID'], fontsize=9,
                verticalalignment='bottom', horizontalalignment='right')

    # Add a title and legend
    ax.set_title("Smoothed Polygon and Data Points")
    ax.legend()

    # Set axis labels
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # Save the plot to a file
    plt.savefig(output_file, dpi=300)

    # Optionally, show the plot (this step can be skipped if you only want to save the file)
    # plt.show()

