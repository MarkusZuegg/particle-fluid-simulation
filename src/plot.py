import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import csv
import sys

def visualize_particles(csv_filepath):
    """
    Reads particle simulation data from a CSV file and creates an animation.

    Args:
        csv_filepath (str): The path to the input CSV file.
    """
    try:
        with open(csv_filepath, 'r') as f:
            # --- 1. Read Header Information ---
            num_particles = int(f.readline().strip())
            radius = float(f.readline().strip())
            xmin = float(f.readline().strip())
            xmax = float(f.readline().strip())
            ymin = float(f.readline().strip())
            ymax = float(f.readline().strip())

            # --- 2. Read Particle Position Data ---
            reader = csv.reader(f)
            frames_data = [list(map(float, row)) for row in reader]

    except FileNotFoundError:
        print(f"Error: The file '{csv_filepath}' was not found.")
        print("Please make sure you have run the C++ simulation first and the file path is correct.")
        return
    except (ValueError, IndexError):
        print(f"Error: The file '{csv_filepath}' is empty or has an invalid format.")
        return

    # Convert the flat list of positions into a structured numpy array
    # Shape: (num_frames, num_particles, 2) where the last dimension is [x, y]
    num_frames = len(frames_data)
    positions = np.array(frames_data).reshape(num_frames, num_particles, 2)

    fig, ax = plt.subplots()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title('Particle Simulation')
    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    fig.set_facecolor('lightgray') # Set background color for the figure
    ax.set_facecolor('white') # Set background color for the plot area

    # The `get_window_extent` is used to convert data units (radius) to point size.
    # This is an approximation but works well for visualization.
    marker_size = (ax.get_window_extent().width / (xmax - xmin) * radius) ** 2
    scatter = ax.scatter(positions[0, :, 0], positions[0, :, 1], s=marker_size)

    # Add a text element to display the current frame number
    frame_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    def update(frame):
        scatter.set_offsets(positions[frame])
        frame_text.set_text(f'Frame: {frame + 1}/{num_frames}')
        return scatter, frame_text

    ani = FuncAnimation(fig, update, frames=num_frames, blit=True, interval=50, repeat=True)

    plt.show()

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <csv_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    visualize_particles(csv_file)

if __name__ == '__main__': main()
