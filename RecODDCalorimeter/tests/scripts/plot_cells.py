import ROOT as root
import pandas as pd
from numpy import sqrt,unique,cos,sin
from matplotlib import pyplot as plt

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

plt.rcParams["figure.figsize"] = [7.00, 7.0]
plt.rcParams["figure.autolayout"] = True
f = root.TFile("ODD_cells_barrel.root")
myTree = f.Get("cells")
columns = ["x", "y", "z", "xA","yA","zA", "xB", "yB", "zB", "xC", "yC", "zC", "xD", "yD", "zD","phi","dphi","R","dR","dz"]
df = pd.read_csv("ECalBarrel_cellPositions.txt", skipinitialspace=True, usecols=columns)
print("Contents in csv file:", df)
fig1, (ax1, ax2) = plt.subplots(1, 2,  width_ratios= (1,3), figsize=(24,7), sharey = True)
fig2, (ax3, ax4) = plt.subplots(1, 2,  width_ratios= (1,3), figsize=(24,7), sharey = True)
ax1.set_xlabel("x (mm)")
ax1.set_ylabel("y (mm)")
ax2.set_xlabel("z (mm)")
ax3.set_xlabel("y (mm)")
ax3.set_ylabel("x (mm)")
ax4.set_xlabel("z (mm)")
errorboxes_xyCart = []
errorboxes_yxCart = []
errorboxes_xyCylindr = []
errorboxes_yxCylindr = []
errorboxes_xyAC = []
for x, y, xa, ya, xc, yc in zip(df.x, df.y, df.xA, df.yA, df.xC, df.yC):
        errorboxes_xyCart.append(Polygon(((x-xc+xa,y-yc+ya),(x+xa+xc,y+ya+yc), (x+xc-xa,y+yc-ya), (x-xa-xc,y-ya-yc))))
        errorboxes_yxCart.append(Polygon(((y-yc+ya,x-xc+xa),(y+ya+yc,x+xa+xc), (y+yc-ya,x+xc-xa), (y-ya-yc,x-xa-xc))))
for event in myTree:
        r, phi, dr, dphi = event.r, event.phi, event.dr, event.dphi
        errorboxes_xyCylindr.append(Polygon((
                ((r-0.5*dr)*cos(phi-0.5*dphi),(r-0.5*dr)*sin(phi-0.5*dphi)), ((r-0.5*dr)*cos(phi+0.5*dphi),(r-0.5*dr)*sin(phi+0.5*dphi)),
                ((r+0.5*dr)*cos(phi+0.5*dphi),(r+0.5*dr)*sin(phi+0.5*dphi)), ((r+0.5*dr)*cos(phi-0.5*dphi),(r+0.5*dr)*sin(phi-0.5*dphi))
                                            )))
        errorboxes_yxCylindr.append(Polygon((
                ((r-0.5*dr)*sin(phi-0.5*dphi),(r-0.5*dr)*cos(phi-0.5*dphi)), ((r-0.5*dr)*sin(phi+0.5*dphi),(r-0.5*dr)*cos(phi+0.5*dphi)),
                ((r+0.5*dr)*sin(phi+0.5*dphi),(r+0.5*dr)*cos(phi+0.5*dphi)), ((r+0.5*dr)*sin(phi-0.5*dphi),(r+0.5*dr)*cos(phi-0.5*dphi))
                                            )))

# Create patch collection with specified colour/alpha
pc_xyCart = PatchCollection(errorboxes_xyCart, facecolor='b', alpha=0.2, edgecolor='white', linestyle="-", linewidth=1)
pc_yxCart = PatchCollection(errorboxes_yxCart, facecolor='b', alpha=0.2, edgecolor='white', linestyle="-", linewidth=1)
pc_xyCylindr = PatchCollection(errorboxes_xyCylindr, facecolor='r', alpha=0.2, edgecolor='r', linestyle="-", linewidth=1)
pc_yxCylindr = PatchCollection(errorboxes_yxCylindr, facecolor='r', alpha=0.2, edgecolor='r', linestyle="-", linewidth=1)
ax1.add_collection(pc_xyCart)
ax1.add_collection(pc_xyCylindr)
ax3.add_collection(pc_yxCart)
ax3.add_collection(pc_yxCylindr)
errorboxes_zyCart = []
errorboxes_zxCart = []
errorboxes_zyCylindr = []
errorboxes_zxCylindr = []
for x, y, xa, ya, xc, yc in zip(df.z, df.y, df.zD, df.yD, df.zB, df.yB):
        errorboxes_zyCart.append(Polygon(((x-xc+xa,y-yc+ya),(x+xa+xc,y+ya+yc), (x+xc-xa,y+yc-ya), (x-xa-xc,y-ya-yc))))
for x, y, xa, ya, xc, yc in zip(df.z, df.x, df.zD, df.xD, df.zB, df.xB):
        errorboxes_zxCart.append(Polygon(((x-xc+xa,y-yc+ya),(x+xa+xc,y+ya+yc), (x+xc-xa,y+yc-ya), (x-xa-xc,y-ya-yc))))
for entry in myTree:
        r, phi, dr, dphi, z, dz = event.r, event.phi, event.dr, event.dphi, event.z, event.dz
        ymax = max ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        ymin = min ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        xmax = max ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        xmin = min ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        errorboxes_zyCylindr.append(Polygon((
                (z+0.5*dz,ymax), (z-0.5*dz,ymax), (z-0.5*dz,ymin), (z+0.5*dz,ymin)
                                            )))
        errorboxes_zxCylindr.append(Polygon((
                (z+0.5*dz,xmax), (z-0.5*dz,xmax), (z-0.5*dz,xmin), (z+0.5*dz,xmin)
                                            )))
pc_zyCart = PatchCollection(errorboxes_zyCart, facecolor='b', alpha=0.2, edgecolor='white')
pc_zxCart = PatchCollection(errorboxes_zxCart, facecolor='b', alpha=0.2, edgecolor='white')
pc_zyCylindr = PatchCollection(errorboxes_zyCylindr, facecolor='r', alpha=0.1, edgecolor='r')
pc_zxCylindr = PatchCollection(errorboxes_zxCylindr, facecolor='r', alpha=0.1, edgecolor='r')
ax1.scatter(df.x, df.y)
ax1.scatter([ev.r*cos(ev.phi) for ev in myTree], [ev.r*sin(ev.phi) for ev in myTree],color='r', marker='+')
ax3.scatter(df.y, df.x)
ax3.scatter([ev.r*sin(ev.phi) for ev in myTree], [ev.r*cos(ev.phi) for ev in myTree],color='r', marker='+')
ax2.add_collection(pc_zyCart)
ax2.add_collection(pc_zyCylindr)
ax2.scatter(df.z, df.y)
ax4.add_collection(pc_zxCart)
ax4.add_collection(pc_zxCylindr)
ax4.scatter(df.z, df.x)
plt.show()
