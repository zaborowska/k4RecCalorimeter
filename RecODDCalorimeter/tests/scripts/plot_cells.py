import ROOT as root
import pandas as pd
from numpy import sqrt,unique,cos,sin
from matplotlib import pyplot as plt

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

plt.rcParams["figure.figsize"] = [7.00, 7.0]
plt.rcParams["figure.autolayout"] = True
f = root.TFile("ODD_cells.root")
myTree = f.Get("cells")
columns = ["x", "y", "z", "xA","yA","zA", "xB", "yB", "zB", "xC", "yC", "zC", "xD", "yD", "zD","phi","dphi","R","dR","dz"]
df = pd.read_csv("ECalBarrel_cellPositions.txt", skipinitialspace=True, usecols=columns)
print("Contents in csv file:", df)
fig1, (ax1, ax2) = plt.subplots(1, 2,  width_ratios= (1,3), figsize=(24,7), sharey = True)
ax1.set_xlabel("x (mm)")
ax1.set_ylabel("y (mm)")
ax2.set_xlabel("z (mm)")
errorboxes_xyCart = []
errorboxes_xyCylindr = []
errorboxes_xyAC = []
for x, y, xa, ya, xc, yc in zip(df.x, df.y, df.xA, df.yA, df.xC, df.yC):
        errorboxes_xyCart.append(Polygon(((x-xc+xa,y-yc+ya),(x+xa+xc,y+ya+yc), (x+xc-xa,y+yc-ya), (x-xa-xc,y-ya-yc))))
for event in myTree:
        r, phi, dr, dphi = event.r, event.phi, event.dr, event.dphi
        errorboxes_xyCylindr.append(Polygon((
                ((r-0.5*dr)*cos(phi-0.5*dphi),(r-0.5*dr)*sin(phi-0.5*dphi)), ((r-0.5*dr)*cos(phi+0.5*dphi),(r-0.5*dr)*sin(phi+0.5*dphi)),
                ((r+0.5*dr)*cos(phi+0.5*dphi),(r+0.5*dr)*sin(phi+0.5*dphi)), ((r+0.5*dr)*cos(phi-0.5*dphi),(r+0.5*dr)*sin(phi-0.5*dphi))
                                            )))

# Create patch collection with specified colour/alpha
pc_xyCart = PatchCollection(errorboxes_xyCart, facecolor='b', alpha=0.2, edgecolor='white', linestyle="-", linewidth=1)
pc_xyCylindr = PatchCollection(errorboxes_xyCylindr, facecolor='r', alpha=0.2, edgecolor='r', linestyle="-", linewidth=1)
ax1.add_collection(pc_xyCart)
ax1.add_collection(pc_xyCylindr)
errorboxes_zyCart = []
errorboxes_zyCylindr = []
for x, y, xa, ya, xc, yc in zip(df.z, df.y, df.zD, df.yD, df.zB, df.yB):
        errorboxes_zyCart.append(Polygon(((x-xc+xa,y-yc+ya),(x+xa+xc,y+ya+yc), (x+xc-xa,y+yc-ya), (x-xa-xc,y-ya-yc))))
for entry in myTree:
        r, phi, dr, dphi, z, dz = event.r, event.phi, event.dr, event.dphi, event.z, event.dz
        ymax = max ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        ymin = min ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        xmax = max ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        xmin = min ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        errorboxes_zyCylindr.append(Polygon((
                (z+0.5*dz,ymax), (z-0.5*dz,ymax), (z-0.5*dz,ymin), (z+0.5*dz,ymin)
                                            )))
pc_zyCart = PatchCollection(errorboxes_zyCart, facecolor='b', alpha=0.2, edgecolor='white')
pc_zyCylindr = PatchCollection(errorboxes_zyCylindr, facecolor='r', alpha=0.1, edgecolor='r')
ax1.scatter([ev.r*cos(ev.phi) for ev in myTree], [ev.r*sin(ev.phi) for ev in myTree],color='r', marker='+')
ax2.add_collection(pc_zyCart)
ax2.add_collection(pc_zyCylindr)
ax2.scatter(df.z, df.y)
plt.show()
