import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
import matplotlib.animation as animation

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="filename of the state file")
    parser.add_argument("side_length", type=float, help="side length of the simulation box")
    args = parser.parse_args()

    def read_sample():
        samples = np.loadtxt(args.filename)

        for sample in samples:
            reshaped = sample.reshape(-1, 3)
            yield [x.flatten() for x in np.hsplit(reshaped, 3)]

    fig, ax = plt.subplots()
    ax.set_xlim(-0.5 * args.side_length, 0.5 * args.side_length)
    ax.set_ylim(-0.5 * args.side_length, 0.5 * args.side_length)
    ax.set_aspect('equal')

    radius = 1.0

    ellipse_pos = EllipseCollection(widths=radius, heights=radius, angles=0, units='xy', facecolors='r', offsets=[],
                                    transOffset=ax.transData)
    ellipse_neg = EllipseCollection(widths=radius, heights=radius, angles=0, units='xy', facecolors='b', offsets=[],
                                    transOffset=ax.transData)
    ax.add_collection(ellipse_pos)
    ax.add_collection(ellipse_neg)

    def update(data):
        q, x, y = data
        idx_pos = q > 0.0
        idx_neg = np.logical_not(idx_pos)

        global ellipse_pos, ellipse_neg
        ellipse_pos.set_offsets(np.column_stack((x[idx_pos], y[idx_pos])))
        ellipse_neg.set_offsets(np.column_stack((x[idx_neg], y[idx_neg])))

    ani = animation.FuncAnimation(fig, update, read_sample, interval=0)
    plt.show()










