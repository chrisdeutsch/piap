import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

samples = np.loadtxt("../data/visualization/out_3d.tsv", delimiter='\t')
#plt.style.use("publication")
q = samples[:,0]
x = samples[:,1]
y = samples[:,2]
z = samples[:,3]

qr = q[0::2]
xr = x[0::2]
yr = y[0::2]
zr = z[0::2]

ql = q[1::2]
xl = x[1::2]
yl = y[1::2]
zl = z[1::2]



def drawSphere(xCenter, yCenter, zCenter, r):
    #draw sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    # shift and scale sphere
    x = r*x + xCenter
    y = r*y + yCenter
    z = r*z + zCenter
    return (x,y,z)


r = 20 * [0.5]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


plt.setp( ax.get_xticklabels(), visible=False)
plt.setp( ax.get_yticklabels(), visible=False)
plt.setp( ax.get_zticklabels(), visible=False)
ax.set_xlim(-0.5 * 15, 0.5 * 15)
ax.set_ylim(-0.5 * 15, 0.5 * 15)
ax.set_zlim(-0.5 * 15, 0.5 * 15)
ax.set_aspect('equal')

ax.view_init(elev=30, azim=-50)

for ii in range(0,180):
	ax.view_init(elev=30., azim=2*ii)

	# draw a sphere for each data point
	for (xi,yi,zi,ri) in zip(xr,yr,zr,r):
		(xs,ys,zs) = drawSphere(xi,yi,zi,ri)
		ax.plot_surface(xs, ys, zs, color="r", linewidth = 0)
		
	for (xi,yi,zi,ri) in zip(xl,yl,zl,r):
		(xs,ys,zs) = drawSphere(xi,yi,zi,ri)
		ax.plot_surface(xs, ys, zs, color="b", linewidth = 0)
	
	anistring = "3D_"
	anistring += str(ii).zfill(3)
	anistring += ".png"
	plt.savefig(anistring, bbox_inches="tight")
