# All made by chatGPT for plotting purposes only
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import numpy as np
import json


def animate_simulation(csv_file, json_file, output_file=None):
    # ==================== Load metadata ====================
    with open(json_file, "r") as f:
        metadata = json.load(f)

    num_particles = metadata["numParticles"]
    radius = metadata["radius"]
    xmin = metadata["xmin"]
    xmax = metadata["xmax"]
    ymin = metadata["ymin"]
    ymax = metadata["ymax"]

    # ==================== Load CSV ====================
    df = pd.read_csv(csv_file)
    steps = df["step"].values
    num_frames = len(steps)

    # Extract x,y columns
    x_cols = [f"x{i}" for i in range(num_particles)]
    y_cols = [f"y{i}" for i in range(num_particles)]

    # ==================== Setup plot ====================
    fig, ax = plt.subplots()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")

    # The size of the points in the scatter plot needs to be calculated based on the display units
    # This is a more robust way to set particle size to match the radius
    # We calculate the radius in points (72 points per inch)
    fig.set_size_inches(6, 6) # Assuming a 6x6 inch figure size
    x_range_inches = (xmax - xmin) / fig.get_figwidth()
    radius_in_points = (radius / x_range_inches) * 72
    
    scat = ax.scatter([], [], s=radius_in_points**2, c="blue")

    title = ax.set_title("Step 0")

    # ==================== Animation function ====================
    def init():
        scat.set_offsets(np.empty((0, 2)))
        title.set_text("Step 0")
        return scat, title

    def update(frame):
        x = df.loc[frame, x_cols].values
        y = df.loc[frame, y_cols].values
        positions = list(zip(x, y))
        scat.set_offsets(positions)
        title.set_text(f"Step {steps[frame]}")
        # Print progress to the console
        print(f"Processing frame {frame + 1}/{num_frames}", end='\r')
        return scat, title

    # The interval determines the speed of the live animation (in milliseconds)
    # The fps in ani.save() determines the speed of the saved file
    ani = FuncAnimation(fig, update, frames=num_frames,
                        init_func=init, blit=True, interval=60, repeat=False)

    # ==================== MODIFICATION START ====================
    if output_file:
        # Determine the writer and fps for saving
        # 1000ms / 60ms_interval = 16.67 fps. We'll use 15.
        fps = 15
        print(f"Saving animation to {output_file} at {fps} fps...")
        
        # The writer depends on the file extension and what you have installed.
        # 'imagemagick' is common for .gif
        # 'ffmpeg' is common for .mp4
        if output_file.endswith('.gif'):
            ani.save(output_file, writer='imagemagick', fps=fps)
        elif output_file.endswith('.mp4'):
            ani.save(output_file, writer='ffmpeg', fps=fps)
        else:
            print("Unsupported file format. Please use .gif or .mp4")
            sys.exit(1)
            
        print(f"\nAnimation successfully saved to {output_file}")
    else:
        # If no output file is specified, show the plot interactively
        plt.show()
    # ==================== MODIFICATION END ====================


if __name__ == "__main__":
    # MODIFICATION: Updated argument checking to allow for an optional output file
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python plot.py <csv_file> <json_file> [output_file]")
        sys.exit(1)

    csv_file = sys.argv[1]
    json_file = sys.argv[2]
    
    # Check if an output file was provided
    output_file = sys.argv[3] if len(sys.argv) == 4 else None

    animate_simulation(csv_file, json_file, output_file)