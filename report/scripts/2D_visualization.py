import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
import matplotlib.animation as animation

plt.style.use("publication")

samples = np.loadtxt("../data/visualization/out.tsv", delimiter='\t')

q = samples[:,0]
x = samples[:,1]
y = samples[:,2]


fig, ax = plt.subplots()
ax.set_xlim(-0.5 * 15, 0.5 * 15)
ax.set_ylim(-0.5 * 15, 0.5 * 15)
ax.set_aspect('equal')

plt.ylabel('position~$y$')
plt.xlabel('position~$x$')


radius = 1.0

ellipse_pos = EllipseCollection(widths=radius, heights=radius, angles=0, units='xy', facecolors='r', offsets=[],
								transOffset=ax.transData)
ellipse_neg = EllipseCollection(widths=radius, heights=radius, angles=0, units='xy', facecolors='b', offsets=[],
								transOffset=ax.transData)
ax.add_collection(ellipse_pos)
ax.add_collection(ellipse_neg)

idx_pos = q > 0.0
idx_neg = np.logical_not(idx_pos)

ellipse_pos.set_offsets(np.column_stack((x[idx_pos], y[idx_pos])))
ellipse_neg.set_offsets(np.column_stack((x[idx_neg], y[idx_neg])))

#plt.scatter(x, y, s = area, c = ["r","b"])
plt.savefig('../figures/Kristall.pdf', bbox_inches='tight')
