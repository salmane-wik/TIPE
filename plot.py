#!/usr/bin/env python3
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import animation as anim
import numpy as np
import sys
import random
import epidemic_model
np.set_printoptions(linewidth=np.inf)
length = 50

fig, ax = plt.subplots(figsize=(15,5))
plt.axes().set_aspect('equal')
ax.set_xlim(0, 3*epidemic_model.N)
ax.set_ylim(0, epidemic_model.N)


arr = np.zeros((epidemic_model.N, 3*epidemic_model.N))

im = plt.imshow(arr, interpolation='nearest',
        cmap = "hot", vmin=0, vmax=1.,
              origin='lower', animated=True) # small changes here

T = 0
def animate(i):
    global x,y,T
    suspected = epidemic_model.model.suspected.reshape(epidemic_model.N, epidemic_model.N)
    infected = epidemic_model.model.infected.reshape(epidemic_model.N, epidemic_model.N)
    recovered = epidemic_model.model.recovered.reshape(epidemic_model.N, epidemic_model.N)
    arr = np.hstack((suspected, infected, recovered))
    arr /= np.max(arr)
    im.set_array(arr)
    old_T = T
    print()
    print(np.sum(suspected))
    print(np.sum(infected))
    print(np.sum(recovered))
    while int(T*60) == int(old_T*60):
        epidemic_model.model.iterate()
        T += epidemic_model.dt
    return [im]

figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
anim = anim.FuncAnimation(fig, animate, frames=200, interval=20, blit=True)
plt.show()