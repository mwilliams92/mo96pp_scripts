import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches

fig, ax = plt.subplots()

x = np.linspace(0,9*np.pi,151)
y = np.sin(x)
ax.plot(x,y, color="green", lw="3")

verts = np.array([[0,1],[0,-1],[2,0],[0,1]]).astype(float)*1.3
verts[:,0] += 9*np.pi
path = mpath.Path(verts)
patch = mpatches.PathPatch(path, fc='green', ec="green")
ax.add_patch(patch)

ax.axis("off")
ax.set_aspect("equal",'datalim')
ax.relim()
ax.autoscale_view()
plt.show()