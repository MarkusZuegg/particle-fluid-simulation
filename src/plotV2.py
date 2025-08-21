import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import csv
import sys
import pandas as pd
import json
import matplotlib.patches as patches

def visualize_particles(csv_filepath, meta_filepath):
    """
    Reads particle simulation data from a CSV file and its corresponding JSON
    metadata file to create a particle animation.

    Args:
        base_filepath (str): The base path for the simulation files (e.g., 'logs/x_sim1').
                             The script will look for '.csv' and '_meta.json' extensions.
    """

    try:
        # --- 1. Read Metadata from JSON file ---
        with open(meta_filepath, 'r') as f:
            metadata = json.load(f)

        num_particles = metadata['numParticles']
        radius = metadata['particleRadius']
        bounds = metadata['bounds']
        xmin = metadata['xmin']
        xmax = metadata['xmax']
        ymax = metadata['ymax']
        ymin = bounds['ymin']

        # --- 2. Read Particle Position Data using Pandas ---
        df = pd.read_csv(csv_filepath)

    except FileNotFoundError:
        print(f"Error: Could not find required files.")
        print(f"Please make sure '{csv_filepath}' and '{meta_filepath}' exist.")
        return
    except (KeyError, json.JSONDecodeError):
        print(f"Error: The metadata file '{meta_filepath}' is corrupted or missing required keys.")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: The CSV file '{csv_filepath}' is empty.")
        return

    # --- 3. Reshape the Data for Animation ---
    # The dataframe has columns like 'step', 'time', 'x0', 'y0', 'vx0', 'vy0', ...
    # We only need the position columns for visualization.
    position_columns = []
    for i in range(num_particles):
        position_columns.extend([f'x{i}', f'y{i}'])
    
    # Select only the position columns and convert to a NumPy array
    try:
        positions_flat = df[position_columns].values
    except KeyError:
        print("Error: CSV file is missing expected position columns (e.g., 'x0', 'y0', ...).")
        return

    # Reshape the (num_frames, num_particles * 2) array into (num_frames, num_particles, 2)
    num_frames = len(positions_flat)
    if num_frames == 0:
        print("Error: No data frames found in the CSV file.")
        return
    positions = positions_flat.reshape(num_frames, num_particles, 2)

    # --- 4. Set up the Plot and Animation ---
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(xmin - radius, xmax + radius)
    ax.set_ylim(ymin - radius, ymax + radius)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('Particle Simulation')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.grid(True)
    
    # Draw the bounding box for context
    rect = patches.Rectangle(
        (xmin, ymin), xmax - xmin, ymax - ymin,
        linewidth=2, edgecolor='r', facecolor='none', label='Bounds'
    )
    ax.add_patch(rect)

    # Calculate an appropriate marker size based on the particle radius
    # This scales the visual size with the plot dimensions
    marker_size = (ax.get_window_extent().width / (xmax - xmin) * radius) ** 2
    
    # Initialize the scatter plot with the first frame of data
    scatter = ax.scatter(positions[0, :, 0], positions[0, :, 1], s=marker_size)

    # Add a text element to display the current time and step
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=12)

    def update(frame):
        """Function called for each frame of the animation."""
        # Update particle positions
        scatter.set_offsets(positions[frame])
        
        # Update the time text from the 'time' column in the DataFrame
        current_time = df['time'].iloc[frame]
        time_text.set_text(f'Time: {current_time:.2f}s')
        
        return scatter, time_text

    # Create and run the animation
    ani = FuncAnimation(fig, update, frames=num_frames, blit=True, interval=50, repeat=True)

    plt.show()

def main():
    if len(sys.argv) < 2:
        print("Usage: python your_script_name.py <base_filepath>")
        print("Example: python your_script_name.py logs/x_sim1")
        sys.exit(1)

    csv_file = sys.argv[1]
    json_file = sys.argv[2]

    visualize_particles(csv_file, json_file)

if __name__ == '__main__':
    main()