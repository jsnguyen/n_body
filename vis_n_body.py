import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

f = open('particles.dat','rb')
t=0
x=[]
y=[]
z=[]

for i in range(10000):
    x.append([])
    y.append([])
    z.append([])

for line in f:
    line = line.split()

    if len(line) == 1:
        t+=1
        continue
    x[t].append(float(line[0]))
    y[t].append(float(line[1]))
    z[t].append(float(line[2]))
f.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
box_size = 1000000.0
for t in range(0,len(x),10):
    ax.scatter(x[t], y[t], z[t])
    ax.set_xlim(-box_size/2.0, box_size/2.0)
    ax.set_ylim(-box_size/2.0, box_size/2.0)
    ax.set_zlim(-box_size/2.0, box_size/2.0)
    plt.savefig('./frames/n_body_'+str(t)+'.png', bbox_inches='tight')
    ax.clear()
