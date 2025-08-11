import csv
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def main():
    # === Handle command-line argument ===
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <csv_file>")
        sys.exit(1)

    csv_file = sys.argv[1]

    # === Read CSV ===
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        header = [next(reader) for _ in range(6)]
        
        num_particles = int(header[0][0])
        radius = float(header[1][0])
        xmin = float(header[2][0])
        xmax = float(header[3][0])
        ymin = float(header[4][0])
        ymax = float(header[5][0])
        
        frames = []
        for row in reader:
            coords = list(map(float, row))
            xs = coords[0::2]
            ys = coords[1::2]
            frames.append((xs, ys))

    # === Setup Figure ===
    fig, ax = plt.subplots()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect('equal', adjustable='box')

    # Marker size in points^2 — scaled from radius in data coords
    marker_size = (radius * 72 / fig.dpi) ** 2 * 3.14 * 100

    scat = ax.scatter([], [], s=marker_size, c='blue', edgecolors='black')

    # === Animation update ===
    def update(frame_index):
        xs, ys = frames[frame_index]
        scat.set_offsets(list(zip(xs, ys)))
        return scat

    ani = animation.FuncAnimation(
        fig, update, frames=len(frames), interval=50, blit=True, repeat=True
    )

    plt.show()

if __name__ == '__main__': main()
