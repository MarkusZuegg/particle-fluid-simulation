import json
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import numpy as np

def animate_simulation(csv_file, json_file):
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

    scat = ax.scatter([], [], s=(radius * 20)**2, c="blue")  # marker size scaled

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
        return scat, title

    ani = FuncAnimation(fig, update, frames=num_frames,
                        init_func=init, blit=True, interval=50, repeat=False)

    plt.show()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot.py <csv_file> <json_file>")
        sys.exit(1)

    csv_file = sys.argv[1]
    json_file = sys.argv[2]

    animate_simulation(csv_file, json_file)