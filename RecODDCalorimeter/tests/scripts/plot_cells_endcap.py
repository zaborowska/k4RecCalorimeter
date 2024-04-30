import ROOT as root
import pandas as pd
from numpy import sqrt,unique,cos,sin
from matplotlib import pyplot as plt

from matplotlib.patches import Rectangle, Polygon
from matplotlib.collections import PatchCollection

plt.rcParams["figure.figsize"] = [7.00, 7.0]
plt.rcParams["figure.autolayout"] = True
dfr = root.RDataFrame("cells", "ODD_cells.root")
dfr = dfr.Filter("isBarrel==0").AsNumpy(["r","phi","dr","dphi","z","dz"])
fig1, (ax1, ax2) = plt.subplots(1, 2,  width_ratios= (1,3), figsize=(24,7), sharey = True)

ax1.set_xlabel("x (mm)")
ax1.set_ylabel("y (mm)")
ax2.set_xlabel("z (mm)")
errorboxes_xyCylindr = []
errorboxes_zyCylindr = []
xBarrel = []
yBarrel = []

for r, phi, dr, dphi, z, dz in zip(dfr["r"], dfr["phi"], dfr["dr"], dfr["dphi"], dfr["z"], dfr["dz"]):
        errorboxes_xyCylindr.append(Polygon((
                ((r-0.5*dr)*cos(phi-0.5*dphi),(r-0.5*dr)*sin(phi-0.5*dphi)), ((r-0.5*dr)*cos(phi+0.5*dphi),(r-0.5*dr)*sin(phi+0.5*dphi)),
                ((r+0.5*dr)*cos(phi+0.5*dphi),(r+0.5*dr)*sin(phi+0.5*dphi)), ((r+0.5*dr)*cos(phi-0.5*dphi),(r+0.5*dr)*sin(phi-0.5*dphi))
        )))
        xBarrel.append(r*cos(phi))
        yBarrel.append(r*sin(phi))
        ymax = max ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        ymin = min ((r+0.5*dr)*sin(phi-0.5*dphi), (r-0.5*dr)*sin(phi-0.5*dphi), (r+0.5*dr)*sin(phi+0.5*dphi), (r-0.5*dr)*sin(phi+0.5*dphi))
        xmax = max ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        xmin = min ((r+0.5*dr)*cos(phi-0.5*dphi), (r-0.5*dr)*cos(phi-0.5*dphi), (r+0.5*dr)*cos(phi+0.5*dphi), (r-0.5*dr)*cos(phi+0.5*dphi))
        errorboxes_zyCylindr.append(Polygon((
                (z+0.5*dz,ymax), (z-0.5*dz,ymax), (z-0.5*dz,ymin), (z+0.5*dz,ymin)
        )))
print(len(xBarrel))
# Create patch collection with specified colour/alpha
pc_xyCylindr = PatchCollection(errorboxes_xyCylindr, facecolor='r', alpha=0.2, edgecolor='r', linestyle="-", linewidth=1)
ax1.add_collection(pc_xyCylindr)
pc_zyCylindr = PatchCollection(errorboxes_zyCylindr, facecolor='r', alpha=0.1, edgecolor='r')
ax1.scatter(xBarrel, yBarrel, color='r', marker='+')
ax2.add_collection(pc_zyCylindr)
plt.show()
