import matplotlib.pyplot as plt
import numpy as np
import sys
def vis(dot,x,y,name,color):
    plt.figure(figsize=(10, 10))
    length=dot.shape[0]
    plt.plot(dot[0:length:2],dot[1:length:2],'o',color=color) 
    i = 0
    while (i<x.shape[0]):
        plt.plot(x[i:i+2],y[i:i+2],color=color)
        i+=2
    i = 0
    for i in range(dot.shape[0]//2):
      plt.text(dot[i*2],dot[i*2+1]+0.5,i)
    plt.axis('off')
    plt.savefig(name)
file=sys.argv[1]
color=sys.argv[2]
name=sys.argv[3]
f=open(file,"r")
verticesNum=int(f.readline())
edgesNum=int(f.readline())
vertices=np.zeros(verticesNum*2)
x=np.zeros(edgesNum*2)
y=np.zeros(edgesNum*2)
for i in range(verticesNum):
  vertArg=f.readline().split()
  vertices[i*2]=float(vertArg[0])
  vertices[i*2+1]=float(vertArg[1])
for i in range(edgesNum):
  edgeArg=f.readline().split()
  x[i*2]=vertices[2*int(edgeArg[0])]
  x[i*2+1]=vertices[2*int(edgeArg[1])]
  y[i*2]=vertices[2*int(edgeArg[0])+1]
  y[i*2+1]=vertices[2*int(edgeArg[1])+1]
vis(vertices,x,y,name,color)